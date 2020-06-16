from ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

run apt-get update -y
run apt-get -y install cmake curl git zip unzip bzip2 wget ftp vim
run apt-get -y install build-essential
run apt-get -y install libblas-dev liblapack-dev libboost-dev libatlas-base-dev
run apt-get -y install gfortran
run apt-get -y install libopenblas-dev
run apt-get -y install python3 python3-numpy python3-matplotlib python3-boto3
run curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh |  bash
run apt-get -y install git-lfs
run git lfs install
