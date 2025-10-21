cmake_minimum_required(VERSION 3.12)
project(nans LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
if (NOT GLEW_LIBRARIES)
    message(STATUS "GLEW not found via find_package, falling back to pkg-config")
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GLEW REQUIRED glew)
    include_directories(${GLEW_INCLUDE_DIRS})
    link_directories(${GLEW_LIBRARY_DIRS})
endif()

find_package(SDL2 REQUIRED)
find_package(assimp REQUIRED)
find_package(Freetype REQUIRED)
find_package(glm REQUIRED)

set(CompilerFlags
    -DSLOW_BUILD -DDRAW_SPHERES -DDRAW_CUBES -DDRAW_FLOOR=1 -DDRAW_COORDINATES=0
    -DDRAW_EPA=0 -DDRAW_MINKOWSKI=0 -DDRAW_WIRE=0 -DLOGGING=1
    -DCC_COLLISIONS=1 -DFLOOR_COLLISIONS=1 -DCS_COLLISIONS=1 -DSS_COLLISIONS=1
    -Wall -Werror
    -Wvla-cxx-extension
    -Wno-unused-function -Wno-write-strings -Wno-unused-variable -g -Wno-null-dereference
    -Wno-unused-but-set-variable
)

if (APPLE)
    set(LinkerFlags -Wl,-rpath,@loader_path)
else()
    set(LinkerFlags -Wl,-rpath,'$ORIGIN')
endif()

add_library(nans SHARED ./code/nans.cpp)
target_compile_options(nans PRIVATE ${CompilerFlags})
target_link_options(nans PRIVATE ${LinkerFlags})
target_include_directories(nans PRIVATE
    ${OPENGL_INCLUDE_DIRS}
    ${SDL2_INCLUDE_DIRS}
    ${ASSIMP_INCLUDE_DIRS}
    ${FREETYPE_INCLUDE_DIRS}
    ${GLM_INCLUDE_DIRS}
    ${GLEW_INCLUDE_DIRS}
)
target_link_libraries(nans PRIVATE
    ${OPENGL_LIBRARIES}
    ${GLEW_LIBRARIES}
)

include(CMakePrintHelpers)
cmake_print_variables(GLEW_INCLUDE_DIRS GLEW_LIBRARIES)


add_executable(nans_exec ./code/sdl_nans.cpp )
target_compile_options(nans_exec PRIVATE ${CompilerFlags})
target_link_options(nans_exec PRIVATE ${LinkerFlags})
target_include_directories(nans_exec PRIVATE
    ${OPENGL_INCLUDE_DIRS}
    ${GLEW_INCLUDE_DIRS}
    ${SDL2_INCLUDE_DIRS}
    ${ASSIMP_INCLUDE_DIRS}
    ${FREETYPE_INCLUDE_DIRS}
    ${GLM_INCLUDE_DIRS}
)
target_link_libraries(nans_exec PRIVATE
    ${OPENGL_LIBRARIES}
    ${GLEW_LIBRARIES}
    ${SDL2_LIBRARIES}
    ${ASSIMP_LIBRARIES}
    ${FREETYPE_LIBRARIES}
    ${GLM_LIBRARIES}
    dl
)

add_dependencies(nans_exec nans)
