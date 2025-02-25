# Instructions on how to deploy a tagged release

## Choose version X.Y.Z
   * git checkout -b VersionX.Y.Z
   * update ChangeLog.md
   * update VERSION to be X.Y.Z
   * git commit
   * git tag -a vX.Y.Z
   * git push
   * git push --tags
   * check CI and optionally merge to master

## Make a new clone with submodules included
   * git clone git@bitbucket.org:berkeleylab/mhm2.git mhm2-vX.Y.Z
   * git checkout VersionX.Y.Z
   * cd mhm2-vX.Y.Z
   * git submodule init
   * git submodule update
   * rm -rf $(find . -name '.git' )
   * cd ..
   * tar -czf mhm2-vX.Y.Z.tar.gz mhm2-vX.Y.Z/

## Test
   * tar -xzf mhm2-vX.Y.Z.tar.gz
   * cd mhm2-vX.Y.Z
   * mkdir build
   * cd build
   * cmake -DCMAKE_INSTALL_PREFIX=install ..
   * make -j install
   * ./install/bin/ci_asm_qual_test.sh
   
## Release
   * upload tar.gz to Downloads section of bitbucket
   * update the wiki page to announce and link to the release: https://bitbucket.org/berkeleylab/mhm2/wiki/Home
   * build and push docker image to robegan21/mhm2 with vX.Y.Z tag and latest

