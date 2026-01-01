#!/bin/bash -l

masif_root=/scratch/shxiao/masif_seed/
IMAGE=$masif_root/masif_mimicry/masif_mimicry.sif

apptainer exec \
  --bind $masif_root:/workspace \
  --pwd /workspace \
  $IMAGE /bin/bash -c '  
    cd /workspace
    curl -fsSL -o Miniconda3.sh "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" && \
    bash Miniconda3.sh -p /workspace/miniconda -b && \
    rm -f Miniconda3.sh

    cd /workspace/masif_mimicry/scripts/bonsai_scripts/
    source /workspace/miniconda/etc/profile.d/conda.sh
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
    conda env create -y -n bonsai --file environment.yml
    
    cd /workspace
    conda activate bonsai
    if [ ! -d "/workspace/EvoEF2" ]; then
        git clone https://github.com/tommyhuangthu/EvoEF2.git
    fi
    cd /workspace/EvoEF2
    chmod +x ./build.sh
    bash ./build.sh

    cd /workspace
    mkdir -p stride
    cd stride
    tar -xzf /workspace/masif_mimicry/scripts/bonsai_scripts/src/stride.tar.gz
    make

    cd /workspace
    tar -xzf /workspace/masif_mimicry/scripts/bonsai_scripts/src/deepTMHMM.tar.gz
'