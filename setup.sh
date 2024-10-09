#!/bin/bash -e

# Notes and todos:
# - Try linking to the normal cuda runtime, instead of the HPC SDK one.
# - Try to directly download the SDK instead of using apt-get.
# - Try FLang instead of NVHPC.

if [ ! -f "/usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg" ]; then
    sudo apt-get update
    sudo apt-get install -y autoconf automake autotools-dev make git m4 libopenblas-dev
    curl https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg
    echo 'deb [signed-by=/usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' | sudo tee /etc/apt/sources.list.d/nvhpc.list
    sudo apt-get update -y
    sudo apt-get install -y nvhpc-24-9
fi
NVARCH=$(uname -s)_$(uname -m)
export NVARCH
NVCOMPILERS=/opt/nvidia/hpc_sdk
export NVCOMPILERS
PATH=$NVCOMPILERS/$NVARCH/24.9/compilers/bin:$PATH
export PATH
export PATH=$NVCOMPILERS/$NVARCH/24.9/comm_libs/mpi/bin:$PATH
printenv
echo "os-version=$(lsb_release -ds | tr " " -)"

export JSON_FORTRAN_DIR="external/json-fortran"

if [ ! -d $JSON_FORTRAN_DIR ]; then
    git clone https://github.com/jacobwilliams/json-fortran.git external/json-fortran
    JSON_FORTRAN_DIR="$(realpath $JSON_FORTRAN_DIR)"
    FC=nvfortran cmake -B$JSON_FORTRAN_DIR/build -S$JSON_FORTRAN_DIR \
        --install-prefix $JSON_FORTRAN_DIR \
        -DCMAKE_BUILD_TYPE=Release \
        -Wno-dev \
        -DUSE_GNU_INSTALL_CONVENTION=ON \
        -DSKIP_DOC_GEN=ON

    cmake --build $JSON_FORTRAN_DIR/build
    cmake --install $JSON_FORTRAN_DIR/build
    rm -rf $JSON_FORTRAN_DIR/build
fi
JSON_FORTRAN_DIR="$(realpath $JSON_FORTRAN_DIR)"

JSON_FORTRAN_LIB=$(find $JSON_FORTRAN_DIR -type d \
    -exec test -f '{}'/libjsonfortran.so \; -print)

export JSON_FORTRAN_LIB="$(realpath $JSON_FORTRAN_LIB)"
export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:$JSON_FORTRAN_LIB/pkgconfig/"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$JSON_FORTRAN_LIB/"

echo "PKG_CONFIG_PATH=$PKG_CONFIG_PATH"
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

[ $(git apply patches/nvhpc_bge.patch 2>/dev/null) ] || echo "Already applied"
./regen.sh
make clean

# Direct version from the HPC SDK
# ./configure FC=nvfortran FCFLAGS="-O3" \
#     --enable-real=dp \
#     --with-cuda=/opt/nvidia/hpc_sdk/Linux_x86_64/24.9/cuda

# Version use the normal cuda runtime
./configure FC=nvfortran FCFLAGS="-O3" \
    --enable-real=dp \
    --with-cuda=/usr/local/cuda

make
