FROM python:3.10-buster as pip
RUN pip install --upgrade pip
RUN pip install --user biopython

FROM ghcr.io/virtool/workflow:5.3.0
WORKDIR /workflow
COPY --from=pip /root/.local /root/.local
COPY workflow.py fixtures.py utils.py ./
