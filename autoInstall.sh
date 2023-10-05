#!/bin/sh

if [ ! -d "build" ]; then
  echo "build folder non-exist, create build folder"
  mkdir build
fi

cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=~/flexi-hybrid/platform \
    -DFLEXI_EDDYVISCOSITY=OFF \
    -DFLEXI_EQNSYSNAME=rans_komega \
    -DFLEXI_EXACT_MASSMATRIX=OFF \
    -DFLEXI_FV=SWITCH \
    -DFLEXI_FV_RECONSTRUCTION=ON \
    -DFLEXI_LIFTING=br1 \
    -DFLEXI_NODETYPE=GAUSS \
    -DFLEXI_PARABOLIC=ON \
    -DFLEXI_PERFORMANCE=OFF \
    -DFLEXI_POLYNOMIAL_DEGREE=N \
    -DFLEXI_PP_LIMITER=OFF \
    -DFLEXI_SPLIT_DG=OFF \
    -DFLEXI_TESTCASE=default \
    -DFLEXI_VISCOSITY=sutherland \
    -DLIBS_BUILD_HDF5=OFF \
    -DHDF5_ROOT=/Users/saayzf/softwares/hdf5 \
    -DHDF5_LIBRARIES=/Users/saayzf/softwares/hdf5/lib \
    -DLIBS_BUILD_MATH_LIB=OFF \
    -DLIBS_USE_MKL=OFF \
    -DLIBS_USE_MPI=ON \
    -DLIBS_USE_OPENMP=OFF \
    -DLIBS_USE_PAPI=OFF \
    -DFLEXI_UNITTESTS=OFF \
    -DPOSTI=ON \
    -DPOSTI_VISU=ON \
    -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc \
    -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++ \
    -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran
make -j 16
make install
