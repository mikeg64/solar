#!/bin/bash


sudo yum update -y
sudo yum groupinstall -y "Development tools"
sudo yum install -y kernel-devel-`uname -r`
wget http://us.download.nvidia.com/XFree86/Linux-x86_64/352.99/NVIDIA-Linux-x86_64-352.99.run
wget http://developer.download.nvidia.com/compute/cuda/7.5/Prod/local_installers/cuda_7.5.18_linux.run
chmod +x NVIDIA-Linux-x86_64-352.99.run
sudo ./NVIDIA-Linux-x86_64-352.99.run
chmod +x cuda_7.5.18_linux.run
sudo ./cuda_7.5.18_linux.run   # Don't install driver, just install CUDA and sample
sudo nvidia-smi -pm 1
sudo nvidia-smi -acp 0
sudo nvidia-smi --auto-boost-permission=0
sudo nvidia-smi -ac 2505,875