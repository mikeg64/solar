#!/bin/bash
# Use the latest ubuntu image as base for the new image
# ubuntu is the image name and latest is a tag that 
# references a particular version of the image.
# In this case latest is always the latest LTS image
# at the time of writing this, it's 16.04.
#FROM ubuntu:latest



mkdir proj
cd proj

git clone https://github.com/littlelogking/virtual-astro
git clone https://github.com/mikeg64/solar
git clone https://github.com/mikeg64/hermes
git clone https://github.com/mikeg64/solar
git clone https://github.com/mikeg2105/datascience-py

cd ..

#install spack
git clone https://github.com/spack/spack
cd spack
git checkout releases/v0.12



#install anaconda
curl -O https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh
bash Anaconda3-2019.03-Linux-x86_64.sh



