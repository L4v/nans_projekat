#! /bin/bash
mkdir -p ../build
pushd ../build

CompilerFlags=( -DSLOW_BUILD -DDRAW_SPHERES -DDRAW_CUBES -DDRAW_FLOOR=1 -DDRAW_COORDINATES=0
		-DDRAW_EPA=0 -DDRAW_MINKOWSKI=0 -DDRAW_WIRE=0 -DLOGGING=1
		-DCC_COLLISIONS=1 -DFLOOR_COLLISIONS=1 -DCS_COLLISIONS=1 -DSS_COLLISIONS=1
		-Wall -Werror -Wl,-rpath,'$ORIGIN'
		-Wno-unused-function -Wno-write-strings -Wno-unused-variable -g -Wno-null-dereference
		-Wno-unused-but-set-variable)
LinkerFlags=(-lGL -lGLEW `sdl2-config --cflags --libs` -ldl -lassimp)

g++ ${CompilerFlags[*]} -shared -fpic ../code/nans.cpp -o nans.so
g++ ${CompilerFlags[*]} ../code/sdl_nans.cpp -o nans ${LinkerFlags[*]}

popd
