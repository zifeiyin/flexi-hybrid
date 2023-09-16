#!/bin/sh

if [ ! -d "build" ]; then
  echo "build folder non-exist, create build folder"
  mkdir build
fi

cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=~/flexi-hybrid/platform \
    -DFLEXI_EDDYVISCOSITY=ON \
    -DFLEXI_EQNSYSNAME=rans_sa \
    -DFLEXI_EXACT_MASSMATRIX=OFF \
    -DFLEXI_FV=SWITCH \
    -DFLEXI_FV_RECONSTRUCTION=ON \
    -DFLEXI_LIFTING=br1 \
    -DFLEXI_NODETYPE=GAUSS \
    -DFLEXI_PARABOLIC=ON \
    -DFLEXI_PERFORMANCE=OFF \
    -DFLEXI_POLYNOMIAL_DEGREE=4 \
    -DFLEXI_PP_LIMITER=OFF \
    -DFLEXI_SPLIT_DG=OFF \
    -DFLEXI_TESTCASE=default \
    -DFLEXI_VISCOSITY=sutherland \
    -DLIBS_BUILD_HDF5=ON \
    -DLIBS_BUILD_MATH_LIB=OFF \
    -DLIBS_USE_MKL=OFF \
    -DLIBS_USE_MPI=ON \
    -DLIBS_USE_OPENMP=OFF \
    -DLIBS_USE_PAPI=OFF \
    -DFLEXI_UNITTESTS=OFF \
    -DPOSTI=OFF
make -j 16
make install
