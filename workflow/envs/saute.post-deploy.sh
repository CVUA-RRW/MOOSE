#!/usr/bin/env bash
set -Eeu

export LIBRARY_PATH=${CONDA_PREFIX}/lib
export CPP_INCLUDE_PATH=${CONDA_PREFIX}/include
export CPLUS_INCLUDE_PATH=${CONDA_PREFIX}/include
export CXX_INCLUDE_PATH=${CONDA_PREFIX}/include


# Need the boost static libraries, re-building form source

cd ${CONDA_PREFIX}
mkdir -p src_cache
cd src_cache
wget https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.bz2
tar --bzip2 -xf boost_1_81_0.tar.bz2
cd boost_1_81_0
./bootstrap.sh --prefix=${CONDA_PREFIX} --exec-prefix=${CONDA_PREFIX}
./b2 install

# Fetch last skesa saute release and build

cd ${CONDA_PREFIX}/src_cache
wget https://github.com/ncbi/SKESA/archive/refs/tags/skesa.2.4.0_saute.1.3.0_2.tar.gz
tar -xzf skesa.2.4.0_saute.1.3.0_2.tar.gz
mv SKESA-skesa.2.4.0_saute.1.3.0_2 ${CONDA_PREFIX}/include/skesa
cd ${CONDA_PREFIX}/include/skesa
make -f Makefile.nongs BOOST_PATH=${CONDA_PREFIX}
mkdir -p ${CONDA_PREFIX}/bin
mv skesa ${CONDA_PREFIX}/bin/
mv saute ${CONDA_PREFIX}/bin/
mv saute_prot ${CONDA_PREFIX}/bin
mv gfa_connector ${CONDA_PREFIX}/bin
mv kmercounter ${CONDA_PREFIX}/bin