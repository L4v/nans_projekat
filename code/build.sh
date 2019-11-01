#! /bin/bash
mkdir -p ../build
pushd ../build

CompilerFlags=(-Wall -Werror -Wl,-rpath,'$ORIGIN'
	       -Wno-unused-function -Wno-write-strings -Wno-unused-variable -g -Wno-null-dereference
	       -Wno-unused-but-set-variable)
LinkerFlags=(-lGL -lGLEW `sdl2-config --cflags --libs` -ldl)

g++ ${CompilerFlags[*]} ../code/sdl_nans.cpp -o nans ${LinkerFlags[*]}

popd
