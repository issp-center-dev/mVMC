#!/bin/sh
#
# #############  Note  ################################
# Before packing, you should clean the GIT directory as
# $ git clean -f -d -x
# #####################################################
#
# Clone all submodules
#
git submodule update -i -r
#
# Version ID
#
major=`awk '$2=="VERSION_MAJOR"{print $3}' src/mVMC/include/version.h`
minor=`awk '$2=="VERSION_MINOR"{print $3}' src/mVMC/include/version.h`
patch=`awk '$2=="VERSION_PATCH"{print $3}' src/mVMC/include/version.h`
vid=`echo ${major}.${minor}.${patch}`
#
mkdir mVMC-${vid}
#
cp -rf * mVMC-${vid}
#
# Build docments
#
cd mVMC-${vid}/doc/jp
make -f makefile_doc_jp
sed -i -e "s/mathjax/pngmath/g" conf.py
make latexpdfja
make html
cd ../en
sed -i -e "s/mathjax/pngmath/g" conf.py
make latexpdfja
make html
cd ../../
#
# Remove some files
#
find ./ -name ".git*" -delete
rm dist.sh
rm -rf mVMC-${vid}
#
# Pack
#
cd ../
tar czvf mVMC-${vid}.tar.gz mVMC-${vid}
rm -rf mVMC-${vid}
