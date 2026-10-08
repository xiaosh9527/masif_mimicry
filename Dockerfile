ARG MASIF_PLATFORM=linux/amd64
FROM --platform=${MASIF_PLATFORM} pymesh/pymesh

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG DEBIAN_FRONTEND=noninteractive
ARG MINICONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# The base image supplies the legacy Python required by TensorFlow 1.12 and
# PyMesh. Keep it off PATH aliases owned by Conda and give it a stable name.
RUN sed -i 's|deb.debian.org|archive.debian.org|g' /etc/apt/sources.list && \
    sed -i 's|security.debian.org|archive.debian.org|g' /etc/apt/sources.list && \
    sed -i '/buster-updates/d' /etc/apt/sources.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
      build-essential cmake curl git libgl1-mesa-glx unzip vim wget dssp && \
    rm -rf /var/lib/apt/lists/* && \
    native_python="$(command -v python)" && \
    "$native_python" --version 2>&1 | grep -E '^Python 3\.6\.' && \
    ln -sf "$native_python" /usr/local/bin/masif-python

# APBS, MSMS, and PDB2PQR.
RUN mkdir -p /install && \
    git clone https://github.com/Electrostatics/apbs-pdb2pqr /install/apbs-pdb2pqr && \
    cd /install/apbs-pdb2pqr && \
    git checkout b3bfeec && \
    git submodule update --init --recursive && \
    cmake -DGET_MSMS=ON apbs && \
    make -j"$(nproc)" && \
    make install && \
    test -x /usr/local/bin/msms && \
    curl -fsSL https://bootstrap.pypa.io/pip/3.6/get-pip.py -o /tmp/get-pip.py && \
    /usr/local/bin/masif-python /tmp/get-pip.py && \
    cd /install/apbs-pdb2pqr/pdb2pqr && \
    python2.7 scons/scons.py install && \
    rm -f /tmp/get-pip.py

ENV MSMS_BIN=/usr/local/bin/msms \
    APBS_BIN=/usr/local/bin/apbs \
    MULTIVALUE_BIN=/usr/local/share/apbs/tools/bin/multivalue \
    PDB2PQR_BIN=/root/pdb2pqr/pdb2pqr.py

# Reduce and USalign.
RUN git clone https://github.com/rlabduke/reduce.git /install/reduce && \
    cd /install/reduce && \
    make install && \
    test -x /usr/local/bin/reduce && \
    git clone https://github.com/pylelab/USalign.git /install/USalign && \
    make -C /install/USalign -j"$(nproc)" && \
    test -x /install/USalign/USalign

ENV USALIGN_BIN=/install/USalign/USalign

# Legacy MaSIF runtime: this is deliberately isolated from Conda.
RUN /usr/local/bin/masif-python -m pip install \
      matplotlib ipython Biopython scikit-learn tensorflow==1.12 \
      networkx open3d==0.8.0.0 dask==1.2.2 packaging \
      python-igraph==0.9.6 requests==2.27.1 pandas==1.1.5

# Modern bonsai runtime. Miniconda is not added to PATH, so `python` remains
# the native MaSIF interpreter everywhere in the image.
RUN wget -q "$MINICONDA_URL" -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm -f /tmp/miniconda.sh

COPY scripts/bonsai_scripts/environment.yml /tmp/bonsai-environment.yml
RUN /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main || true && \
    /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r || true && \
    /opt/conda/bin/conda env create --file /tmp/bonsai-environment.yml && \
    /opt/conda/bin/conda clean --all --yes && \
    rm -f /tmp/bonsai-environment.yml

# Install the user-facing CLI into the same Python 3.11 environment used by
# bonsai. This is the second and final Python environment in the image.
COPY pyproject.toml README.md /opt/masif-mimicry-cli/
COPY src /opt/masif-mimicry-cli/src
RUN /opt/conda/envs/bonsai/bin/python -m pip install '/opt/masif-mimicry-cli[ui]' && \
    ln -sf /opt/conda/envs/bonsai/bin/masif-mimicry /usr/local/bin/masif-mimicry

# Bonsai native executables and precomputed DeepTMHMM data.
COPY scripts/bonsai_scripts/src/stride.tar.gz /tmp/bonsai-src/stride.tar.gz
COPY scripts/bonsai_scripts/src/deepTMHMM.tar.gz /tmp/bonsai-src/deepTMHMM.tar.gz
RUN git clone https://github.com/tommyhuangthu/EvoEF2.git /opt/EvoEF2 && \
    chmod +x /opt/EvoEF2/build.sh && \
    cd /opt/EvoEF2 && \
    /bin/bash ./build.sh && \
    mkdir -p /opt/stride && \
    tar -xzf /tmp/bonsai-src/stride.tar.gz -C /opt/stride && \
    make -C /opt/stride -j"$(nproc)" && \
    tar -xzf /tmp/bonsai-src/deepTMHMM.tar.gz -C /opt && \
    rm -rf /tmp/bonsai-src

ENV MASIF_PYTHON=/usr/local/bin/masif-python \
    BONSAI_PYTHON=/opt/conda/envs/bonsai/bin/python \
    DEEPTMHMM_DIR=/opt/deepTMHMM/out \
    STRIDE_EXEC=/opt/stride/stride \
    EVOEF2_DIR=/opt/EvoEF2

# Fail the image build if either side of the compatibility boundary is wrong.
RUN /usr/local/bin/masif-python -c \
      'import sys; assert sys.version_info[:2] == (3, 6), sys.version' && \
    /opt/conda/envs/bonsai/bin/python -c \
      'import sys; assert sys.version_info[:2] == (3, 11), sys.version; import Bio, numpy, pandas, rdkit' && \
    /usr/local/bin/masif-mimicry --version

WORKDIR /workspace
CMD ["bash"]
