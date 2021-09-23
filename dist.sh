#!/bin/sh

# #############  README  ##############
# This makes an archive file including static copy of submodules (StdFace)
# and PDF formatted documents docs/mVMC_(ja|en).pdf
# The output filename is mVMC-${vid}.tar.gz,
# where ${vid} is the version number such as 2.1.0 .
# Before using this, install the following python packages:
#   sphinx
#   sphinx_numfig
#   sphinxcontib_spelling
#   git-archive-all
# #####################################

set -e

# Retrieve Version ID
major=`awk '$2=="VERSION_MAJOR"{print $3}' src/mVMC/include/version.h`
minor=`awk '$2=="VERSION_MINOR"{print $3}' src/mVMC/include/version.h`
patch=`awk '$2=="VERSION_PATCH"{print $3}' src/mVMC/include/version.h`
vid=`echo ${major}.${minor}.${patch}`

# Build PDF docs
rm -rf build-docs
mkdir build-docs
cd build-docs
cmake -DDocument=ON ../
make doc-ja-pdf
make doc-en-pdf
cp doc/ja/source/pdf/mVMC.pdf ../doc/mVMC-${vid}_ja.pdf
cp doc/en/source/pdf/mVMC.pdf ../doc/mVMC-${vid}_en.pdf
cd ../
rm -rf build-docs

# Make archive
git-archive-all \
  --extra=doc/mVMC_ja.pdf \
  --extra=doc/mVMC_en.pdf \
  --prefix=mVMC-${vid} \
  mVMC-${vid}.tar.gz
