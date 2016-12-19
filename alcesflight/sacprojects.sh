#!/bin/bash

#script to upload an install projects in an alces flight installation

cd ~
mkdir proj


svn checkout --username mikeg64 https://ccpforge.cse.rl.ac.uk/svn/sac/branches/sac_working
svn checkout --username mikeg64 https://ccpforge.cse.rl.ac.uk/svn/sac/dev/smaug


git clone https://github.com/mikeg64/athena_pmode
git clone https://github.com/mikeg64/hermes
git clone https://github.com/mikeg64/solar



