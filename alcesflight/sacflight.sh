#!/bin/bash

#for install instructions see
# https://github.com/alces-software/clusterware



#script to be run as root
# sudo ./sacflight.sh


yum -y install update
yum -y install svn
yum -y install git
yum -y install gcc-gfortran

#first install alces clusterware
cd /usr/src
git clone https://github.com/alces-software/clusterware

export cw_BUILD_source_dir=/usr/src/clusterware
export cw_DIST=el7

/usr/src/clusterware/scripts/bootstrap
   
source /etc/profile.d/alces-clusterware.sh


#setup alces flight
alces handler enable cluster-sge
alces handler enable task-session
alces handler enable task-session-gnome
alces handler enable clusterable
alces handler enable clusterable-aws-compat
alces handler enable cluster-appliances
alces handler enable cluster-customizer
alces handler enable cluster-firewall
alces handler enable gridware
alces handler enable cluster-gridware
alces handler enable cluster-nfs
alces handler enable flight
alces handler enable taskable

#enable and install services
alces service enable s3cmd
alces service install s3cmd
alces service install aws
alces service install clusterware-dropbox-cli
alces service install alces-flight-www

#set up packages using alces gridware
alces gridware install compilers/gcc
alces gridware install ffmpeg
alces gridware install imagemagick
alces gridware install graphicsmagick
alces gridware install grace
alces gridware install paraview
alces gridware install nvidia-cuda
alces gridware install anaconda3

#spack software for package management
git clone https://github.com/LLNL/spack
cd /usr/src/spack/bin
./spack install libelf

cd ~


