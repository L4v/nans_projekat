#! /bin/bash
mkdir -p code build > /dev/null
pushd code > /dev/null
emacs -q -l ../misc/init.el & > /dev/null
popd > /dev/null
cd build
