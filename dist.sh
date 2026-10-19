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
major=`awk '$2=="MVMC_VERSION_MAJOR"{print $3}' src/mVMC/include/version.h`
minor=`awk '$2=="MVMC_VERSION_MINOR"{print $3}' src/mVMC/include/version.h`
patch=`awk '$2=="MVMC_VERSION_PATCH"{print $3}' src/mVMC/include/version.h`
pre=`awk '$2=="MVMC_VERSION_PRERELEASE"{gsub(/"/,"",$3); print $3}' src/mVMC/include/version.h`
vid=`echo ${major}.${minor}.${patch}`
if [ -n "${pre}" ]; then
  vid=${vid}-${pre}
fi

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
  --extra=doc/mVMC-${vid}_ja.pdf \
  --extra=doc/mVMC-${vid}_en.pdf \
  --prefix=mVMC-${vid} \
  mVMC-${vid}.tar.gz

# Write the hash of the commit into cmake/git_archive.txt of the tarball,
# unless it is filled in already. "vmc.out -v" built from the tarball prints it.
hash=`git rev-parse HEAD`
tmpdir=`mktemp -d`
tar xzf mVMC-${vid}.tar.gz -C ${tmpdir}
if grep -q Format ${tmpdir}/mVMC-${vid}/cmake/git_archive.txt; then
  echo ${hash} > ${tmpdir}/mVMC-${vid}/cmake/git_archive.txt
  COPYFILE_DISABLE=1 tar czf mVMC-${vid}.tar.gz -C ${tmpdir} mVMC-${vid}
fi
rm -rf ${tmpdir}
if [ -n "`git status --porcelain --untracked-files=no`" ]; then
  echo 'WARNING: the source has changes which are not committed.'
  echo '         The hash written into the tarball is that of HEAD.'
fi
