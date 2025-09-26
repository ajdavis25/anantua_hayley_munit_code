#!/usr/bin/env bash
# install system deps
sudo apt-get update
sudo apt-get install -y git cmake build-essential

# clone ipole into runner’s search_dir (~Desktop/ipole_versions) and build it
mkdir -p ~/Desktop/ipole_versions
cd ~/Desktop/ipole_versions
git clone https://github.com/AFD-Illinois/ipole.git
cd ipole
mkdir build && cd build
cmake ..
make -j$(nproc)

# python env + deps
cd "$(dirname "$0")"      # back to project root
python3 -m venv venv      # create venv
source venv/bin/activate
pip install -r requirements.txt
