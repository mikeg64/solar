#!/bin/bash
# Use the latest ubuntu image as base for the new image
# ubuntu is the image name and latest is a tag that 
# references a particular version of the image.
# In this case latest is always the latest LTS image
# at the time of writing this, it's 16.04.
#FROM ubuntu:latest


#
apt -y update
apt -y upgrade
apt install -y git-core 
apt install -y curl 
apt install -y wget
apt install -y build-essential




#RUN apt-get install -y firefox
# Run a system update to get it up to speed
# Then install python3 and pip
apt -y update 

#install codeblocks
apt -y upgrade
apt -y install g++
apt -y install codeblocks

#mkdir proj
#cd proj

#git clone https://github.com/littlelogking/virtual-astro
#git clone https://github.com/mikeg64/solar
#git clone https://github.com/mikeg64/hermes
#git clone https://github.com/mikeg64/solar
#git clone https://github.com/mikeg2105/datascience-py

#cd ..

#install spack
#git clone https://github.com/spack/spack
#cd spack
#git checkout releases/v0.12

#install x2go
apt -y install python-software-properties
apt -y install software-properties-common
add-apt-repository ppa:x2go/stable
apt -y update
apt -y install x2goserver x2goserver-xsession

#install docker
apt -y remove docker docker-engine docker.io containerd runc
apt -y install \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg-agent \
    software-properties-common
	
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -


add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"
   
apt -y update
apt -y install docker-ce docker-ce-cli containerd.io


#install google chrome
wget https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb
dpkg -i google-chrome-stable_current_amd64.deb

#install anaconda
curl -O https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh
bash Anaconda3-2019.03-Linux-x86_64.sh



