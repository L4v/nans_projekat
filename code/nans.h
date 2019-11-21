#ifndef NANS_H

#include <stdint.h>
#include <cstddef>
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/normal.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <cstdio>
#include <vector>
#include <time.h>

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
#define Pi32 3.14159265359f
#define ArrayCount(Array) (sizeof(Array) / sizeof((Array)[0]))

#define GRAVITY_ACCEL  9.81f
#define GLOBAL_FRICTION  1.5f

#define DEFAULT_WINDOW_WIDTH 1024
#define DEFAULT_WINDOW_HEIGHT 768
#define MAX_CUBE_COUNT 16
#define MAX_SPHERE_COUNT 16
#define MAX_GJK_ITERATIONS 64
#define MAX_EPA_ERROR 0.0001f
#define MAX_EPA_ITERATIONS 64

#if SLOW_BUILD
#define Assert(Expression)			\
  if(!(Expression)) {*(int*)0 = 0;}
#else
#define Assert(Expression)
#endif


  enum collision_type
    {
     /* NOTE(Jovan): Types of collisions
	CC - CUBE x CUBE
	CS - CUBE x SPHERE
	SS - SPHERE x SPHERE
     */
     CC = 0,
     CS,
     SS
    };

enum evolve_result
  {
   NoIntersection = 0,
   FoundIntersection,
   StillEvolving
  };
  
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
    sdl_button_state Buttons[12];
    struct{
      sdl_button_state MoveForward;
      sdl_button_state MoveBack;
      sdl_button_state MoveLeft;
      sdl_button_state MoveRight;
      sdl_button_state ShootAction;

      // NOTE(Jovan): Moving the cube for debugging
      sdl_button_state DebugUp;
      sdl_button_state DebugDown;
      sdl_button_state DebugLeft;
      sdl_button_state DebugRight;
      sdl_button_state DebugForward;
      sdl_button_state DebugBack;
      sdl_button_state DebugReset;
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
  uint32 Shaders[3];
  uint32 Textures[3];
  uint32 VAOs[4];
  uint32 VBOs[1];
  uint32 EBOs[1];
  uint32* Indices;
  uint32 Num;
};

struct memory
{
  void* PermanentStorage;
  uint64 PermanentStorageSize;
  void* TransientStorage;
  uint64 TransientStorageSize;

  bool32 IsInitialized;
};

struct memory_arena
{
  memory_index Size;
  uint8* Base;
  memory_index Used;
};

// NOTE(Jovan): Vertex struct, for EPA mostly
struct vertex
{
  glm::vec3 P;
  glm::vec3 SupA;
  
  inline bool operator==(const vertex& R)
    {
     return (this->P == R.P);
    }
};

struct simplex
{
  uint32 Count;
  vertex* Vertices;
};

struct edge
{
  uint32 Count;
  vertex* A;
  vertex* B;

};

struct triangle
{
  uint32 Count;
  vertex* A;
  vertex* B;
  vertex* C;
  glm::vec3* N;
};

struct sdl_camera
{
  real32 FOV;
  real32 Pitch;
  real32 Yaw;
  real32 Speed;
  
  glm::vec3 Position;
  glm::vec3 Target;
  glm::vec3 Direction;
  
  glm::vec3 Up;
  glm::vec3 Front;
  glm::vec3 Right;
};

struct cube
{
  glm::vec3 Position;
  glm::vec3 Velocity;
  glm::vec3 Forces;
  glm::vec3 Vertices[8];
  real32 Size;
  real32 Mass;
  real32 XAngle;
  real32 YAngle;
  real32 ZAngle;
};

struct sphere
{
  glm::vec3 Position;
  real32 Radius;
  real32 XAngle;
  real32 YAngle;
  real32 ZAngle;
};

struct sdl_state
{
  // NOTE(Jovan): Vertices for the simplex
  memory_arena SimplexArena;
  memory_arena EdgeArena;
  memory_arena TriangleArena;
  simplex* Simplex;
  edge *Edge;
  triangle *Triangle;
  
  cube Cubes[MAX_SPHERE_COUNT];
  sphere Spheres[MAX_SPHERE_COUNT];
  sdl_camera Camera;
  uint32 CubeCount;
  uint32 SphereCount;
  uint32 GJKIteration;
};

// NOTE(Jovan): Return value for ODEs and physics functions
// TODO(Jovan): Maybe use other method?
struct phys_return
{
  glm::vec3 X;
  glm::vec3 Y;
};

#define SIM_UPDATE_AND_RENDER(name) void name(memory* Memory, sdl_input* Input, sdl_render* Render, real32 dt)
typedef SIM_UPDATE_AND_RENDER(sim_update_and_render);
SIM_UPDATE_AND_RENDER(SimUpdateAndRenderStub)
{
}


#define NANS_H
#endif
