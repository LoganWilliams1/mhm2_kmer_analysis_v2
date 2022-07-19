# Instructions on how to deploy a tagged release

## Choose version X.Y.Z
   * update ChangeLog.md
   * update VERSION to be X.Y.Z
   * git commit
   * git tag -a vX.Y.Z
   * git push
   * git push -a

## Make a new clone with submodules included
   * git clone git@bitbucket.org:berkeleylab/mhm2.git mhm2-vX.Y.Z
   * cd mhm2-vX.Y.Z
   * git submodule init
   * git submodule update
   * rm -rf $(find . -name '.git' )
   * cd ..
   * tar -czf mhm2-vX.Y.Z.tar.gz mhm2-vX.Y.Z/
   * upload tar.gz to Downloads section of bitbucket

## Test
   * tar -xzf mhm2-vX.Y.Z.tar.gz
   * cd mhm2-vX.Y.Z
   * mkdir build
   * cd build
   * cmake -DCMAKE_INSTALL_PREFIX=install ..
   * make -j install
   * ./install/bin/ci_asm_qual_test.sh
