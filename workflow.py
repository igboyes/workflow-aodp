from logging import getLogger
from pathlib import Path
from typing import List, Any

import aiofiles
from pyfixtures import fixture
from virtool_core.utils import file_length
from virtool_workflow import step, hooks
from virtool_workflow.data_model.indexes import WFIndex

import utils

AODP_MAX_HOMOLOGY = 0
AODP_OLIGO_SIZE = 8

logger = getLogger("workflow")


@fixture
def index(indexes: List[WFIndex]):
    return indexes[0]


@hooks.on_result
async def upload_results(results, analysis_provider):
    await analysis_provider.upload_result(results)


@step
async def join_reads(
    joined_path: Path, proc: int, run_subprocess, results: dict, sample, work_path: Path
):
    """
    Join overlapping paired reads into single reads.

    """
    max_overlap = round(0.65 * sample.read_length)

    command = [
        "flash",
        "--max-overlap",
        str(max_overlap),
        "-d",
        work_path,
        "-o",
        "flash",
        "-t",
        str(proc - 1),
        sample.read_paths,
    ]

    await run_subprocess(command)

    hist_path = work_path / "flash.hist"
    remainder_path = work_path / "flash.notCombined_1.fastq"

    results.update(
        {
            "join_histogram": await utils.parse_flash_histogram(hist_path),
            "joined_pair_count": await file_length(joined_path) / 4,
            "remainder_pair_count": await file_length(remainder_path) / 4,
        }
    )


@step
async def deduplicate_reads(
    joined_path: Path,
    results: dict,
    run_in_executor,
    unique_path: Path,
):
    """Remove duplicate reads. Store the counts for unique reads."""
    counts = await run_in_executor(utils.run_deduplication, joined_path, unique_path)

    results["sequence_counts"] = counts


@step
async def aodp(
    index: WFIndex,
    proc: int,
    results: dict[str, Any],
    run_subprocess,
    unique_path: Path,
    work_path: Path,
):
    """
    Run AODP, parse the output, and update the.

    TODO: Upload data and finalize analysis

    """
    aodp_path = work_path / "aodp.out"
    base_name = work_path / "aodp"

    await run_subprocess(
        [
            "aodp",
            f"--basename={base_name}",
            f"--threads={proc}",
            f"--oligo-size={AODP_OLIGO_SIZE}",
            f"--match={unique_path}",
            f"--match-output={aodp_path}",
            f"--max-homolo={AODP_MAX_HOMOLOGY}",
            index.fasta_path,
        ],
        cwd=work_path,
    )

    parsed = list()

    async with aiofiles.open(aodp_path, "r") as f:
        async for line in f:
            split = line.rstrip().split("\t")
            assert len(split) == 7

            sequence_id = split[1]

            if sequence_id == "-":
                continue

            identity = split[2]

            if identity[0] == "<":
                continue
            else:
                identity = float(identity.replace("%", ""))

            read_id = split[0]

            sequence_id = split[1]

            otu_id = index.get_otu_id_by_sequence_id(sequence_id)
            otu_version = index.get_otu_version_by_sequence_id(sequence_id)

            parsed.append(
                {
                    "id": read_id,
                    "sequence_id": sequence_id,
                    "identity": identity,
                    "matched_length": int(split[3]),
                    "read_length": int(split[4]),
                    "min_cluster": int(split[5]),
                    "max_cluster": int(split[6]),
                    "count": results["sequence_counts"][read_id],
                    "otu": {"version": otu_version, "id": otu_id},
                }
            )

    results["hits"] = parsed
