#ifndef NANS_H

#include <stdint.h>
#include <cstddef>
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <cstdio>

#define internal static
#define global_variable static
#define local_persist static

typedef int8_t int8;
typedef int16_t int16;
typedef int32_t int32;
typedef int64_t int64;
typedef int32 bool32;

typedef uint8_t uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;

typedef size_t memory_index;

typedef float real32;
typedef double real64;

// NOTE(Jovan): Just for testing purposes
// #include "matmath.h"

#define Kibibytes(Value) ((Value) * 1024LL)
#define Mebibytes(Value) (Kibibytes(Value) * 1024LL)
#define Gibibytes(Value) (Mebibytes(Value) * 1024LL)

#define ArrayCount(Array) (sizeof(Array) / sizeof((Array)[0]))


#define DEFAULT_WINDOW_WIDTH 1024
#define DEFAULT_WINDOW_HEIGHT 768
#define MAX_CUBE_COUNT 64


#if SLOW_BUILD
#define Assert(Expression)			\
  if(!(Expression)) {*(int*)0 = 0;}
#else
#define Assert(Expression)
#endif


struct sdl_mouse
{
  real32 Sensitivity;
  int32 XRel;
  int32 YRel;
  int32 X;
  int32 Y;
};

struct sdl_button_state
{
  int32 HalfTransitionCount;
  bool32 EndedDown;
};

struct sdl_keyboard
{
  union
  {
    sdl_button_state Buttons[5];
    struct{
      sdl_button_state MoveForward;
      sdl_button_state MoveBack;
      sdl_button_state MoveLeft;
      sdl_button_state MoveRight;
      sdl_button_state ShootAction;
    };
  };
};

struct sdl_input
{
  sdl_keyboard KeyboardController;
  sdl_mouse MouseController;
};

struct sdl_render
{
  uint32 Shaders[1];
  uint32 Textures[1];
  uint32 VAOs[1];
};

struct memory
{
  void* PermanentStorage;
  uint64 PermanentStorageSize;
  void* TransientStorage;
  uint64 TransientStorageSize;

  bool32 IsInitialized;
};

struct sdl_camera
{
  real32 FOV;
  real32 Pitch;
  real32 Yaw;
  
  glm::vec3 Position;
  glm::vec3 Target;
  glm::vec3 Direction;
  
  glm::vec3 Up;
  glm::vec3 Front;
  glm::vec3 Right;
};

struct sdl_state
{
  glm::vec3 Positions[10];
  sdl_camera Camera;
};

#define SIM_UPDATE_AND_RENDER(name) void name(memory* Memory, sdl_input* Input, sdl_render* Render, real32 dt)
typedef SIM_UPDATE_AND_RENDER(sim_update_and_render);
SIM_UPDATE_AND_RENDER(SimUpdateAndRenderStub)
{
}


#define NANS_H
#endif
