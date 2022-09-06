#!/usr/bin/env bash

if [ -n "$MHM2_BUILD_ENV" ]; then
    source $MHM2_BUILD_ENV
fi

upcxx_exec=`which upcxx`

if [ -z "$upcxx_exec" ]; then
    echo "upcxx not found. Please install or set path."
    exit 1
fi

upcxx_exec_canonical=$(readlink -f $upcxx_exec)
if [ "$upcxx_exec_canonical" != "$upcxx_exec" ]; then
    echo "Found symlink for upcxx - using target at $upcxx_exec_canonical"
    export PATH=`dirname $upcxx_exec_canonical`:$PATH
fi

set -e

SECONDS=0

rootdir=`pwd`

INSTALL_PATH=${MHM2_INSTALL_PATH:=$rootdir/install}

echo "Installing into $INSTALL_PATH"

BINARY="${MHM2_BINARY:=mhm2}"

rm -rf $INSTALL_PATH/bin/mhm2
rm -rf $INSTALL_PATH/bin/${BINARY}

if [ "$1" == "cleanall" ]; then
    rm -rf .build/*
    # if this isn't removed then the the rebuild will not work
    rm -rf $INSTALL_PATH/cmake
    exit 0
elif [ "$1" == "clean" ]; then
    rm -rf .build/bin .build/CMake* .build/lib* .build/makeVersionFile .build/src .build/test .build/cmake* .build/Makefile .build/make*
    # if this isn't removed then the the rebuild will not work
    rm -rf $INSTALL_PATH/cmake
    exit 0
else
    mkdir -p $rootdir/.build
    cd $rootdir/.build
    if [ "$1" == "Debug" ] || [ "$1" == "Release" ] || [ "$1" == "RelWithDebInfo" ]; then
        rm -rf *
        rm -rf $INSTALL_PATH/cmake
        cmake $rootdir -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=$1 -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH \
              -DMHM2_ENABLE_TESTING=0 $MHM2_CMAKE_EXTRAS $2
        #-DENABLE_CUDA=0
    fi
    make -j ${MHM2_BUILD_THREADS} all install
    #make VERBOSE=1 -j ${MHM2_BUILD_THREADS} all install
    # make -j ${MHM2_BUILD_THREADS} check
    if [ "$BINARY" != "mhm2" ]; then
        mv -f $INSTALL_PATH/bin/mhm2 $INSTALL_PATH/bin/${BINARY}
    fi
fi

echo "Build took $((SECONDS))s"
