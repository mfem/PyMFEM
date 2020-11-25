#!/bin/bash
set -ev

UBUNTU_VERSION=$(lsb_release -sr)
UBUNTU_VERSION=ubuntu"${UBUNTU_VERSION//.}"
CUDA=10.1.105-1
CUDA_SHORT=10.1

INSTALLER=cuda-repo-${UBUNTU_VERSION}_${CUDA}_amd64.deb
wget http://developer.download.nvidia.com/compute/cuda/repos/${UBUNTU_VERSION}/x86_64/${INSTALLER}
sudo dpkg -i ${INSTALLER}
wget https://developer.download.nvidia.com/compute/cuda/repos/${UBUNTU_VERSION}/x86_64/7fa2af80.pub
sudo apt-key add 7fa2af80.pub
sudo apt update -qq
sudo apt install -y cuda-core-${CUDA_SHORT/./-} cuda-cudart-dev-${CUDA_SHORT/./-} cuda-cufft-dev-${CUDA_SHORT/./-}
sudo apt clean



