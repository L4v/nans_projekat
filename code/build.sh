#! /bin/bash
mkdir -p ../build
pushd ../build

CompilerFlags=( -DSLOW_BUILD -DDRAW_SPHERES -DDRAW_CUBES -DDRAW_FLOOR=1 -DDRAW_COORDINATES=0
		-DDRAW_EPA=1 -DDRAW_MINKOWSKI=1 -DDRAW_WIRE=1 -DLOGGING=1
		-Wall -Werror -Wl,-rpath,'$ORIGIN'
		-Wno-unused-function -Wno-write-strings -Wno-unused-variable -g -Wno-null-dereference
		-Wno-unused-but-set-variable)
LinkerFlags=(-lGL -lGLEW `sdl2-config --cflags --libs` -ldl)

g++ ${CompilerFlags[*]} -shared -fpic ../code/nans.cpp -o nans.so ${LinkerFlags[*]}
g++ ${CompilerFlags[*]} ../code/sdl_nans.cpp -o nans ${LinkerFlags[*]}

popd
