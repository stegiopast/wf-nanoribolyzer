FROM docker.io/mambaorg/micromamba:1.5-jammy

USER root

COPY remora.yaml /opt/conda.yml

RUN micromamba install -y -n base -f /opt/conda.yml \
    && micromamba clean -afy \
    && micromamba env export -n base --explicit > /opt/conda.lock

ENV PATH "${MAMBA_ROOT_PREFIX}/bin:${PATH}"

USER $MAMBA_USER