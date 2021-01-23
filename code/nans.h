#ifndef NANS_H

#include <stdint.h>
#include <cstddef>
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/normal.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/matrix_operation.hpp>
#include <cstdio>
#include <vector>
#include <time.h>
#include <math.h>
// TODO(Jovan): Fix this bug with freetype2
/* #define internal static */
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

#define Kibibytes(Value) ((Value)*1024LL)
#define Mebibytes(Value) (Kibibytes(Value) * 1024LL)
#define Gibibytes(Value) (Mebibytes(Value) * 1024LL)
#define Pi32 3.14159265359f
#define ArrayCount(Array) (sizeof(Array) / sizeof((Array)[0]))

#define GRAVITY_ACCEL 9.81f
#define GLOBAL_FRICTION 1.5f

#define SHOOT_FORCE 4001.0f

#define DEFAULT_WINDOW_WIDTH 1920
#define DEFAULT_WINDOW_HEIGHT 1080
#define MAX_CUBE_COUNT 16
#define MAX_SPHERE_COUNT 16
#define MAX_GJK_ITERATIONS 64
#define MAX_EPA_ERROR 0.001f
#define MAX_EPA_ITERATIONS 64

#define POINT_LIGHT_COUNT 4


#if SLOW_BUILD
#define Assert(Expression) \
    if (!(Expression))     \
    {                      \
        *(int *)0 = 0;     \
    }
#else
#define Assert(Expression)
#endif

#define PrintVector(Vector) printf("%f %f %f\n", Vector.x, Vector.y, Vector.z)

enum collision_type
{
    /* NOTE(Jovan): Types of collisions
	CC - CUBE x CUBE
	CS - CUBE x SPHERE
	SS - SPHERE x SPHERE
	CF - CUBE x FLOOR
	SF - SPHERE x FLOOR
	
     */
    CC = 0,
    CS,
    CF,
    //
    SS,
    SF
};

enum evolve_result
{
    NoIntersection = 0,
    FoundIntersection,
    StillEvolving
};

struct light
{
    glm::vec3 Position;
    real32 Size;
    
    real32 Kc;
    real32 Kl;
    real32 Kq;

    glm::vec3 Ambient;
    glm::vec3 Diffuse;
    glm::vec3 Specular;
};

struct character
{
    uint32 TextureId;
    glm::ivec2 Size;
    glm::ivec2 Bearing;
    long int Advance;
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
        sdl_button_state Buttons[13];
        struct
        {
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
            sdl_button_state DebugContinue;
        };
    };
};

struct sdl_input
{
    sdl_keyboard KeyboardController;
    sdl_mouse MouseController;
};

enum sdl_texture
{
    CONTAINER_DIFFUSE = 0,
    CONTAINER_SPECULAR,
    EARTH,
    CHECKERBOARD,
    AMONG_US_DIFFUSE,
    AMONG_US_SPECULAR,
    METAL_ALBEDO,
    METAL_SPECULAR,
    TEXTURE_COUNT
};

enum sdl_vao
{
    CUBEVAO = 0,
    SPHEREVAO,
    FLOORVAO,
    TEXTVAO,
    LIGHTVAO,
    MODELVAO,
    VAO_COUNT
};

enum sdl_vbo
{
    CUBEVBO = 0,
    SPHEREVBO,
    FLOORVBO,
    TEXTVBO,
    LIGHTVBO,
    MODELVBO,
    VBO_COUNT
};

enum sdl_shader
{
    SIMPLE_COLOR_SH,
    TEXT_SH,
    LIGHTING_SH,
    SHADER_COUNT
};

struct sdl_render
{
    uint32 Shaders[SHADER_COUNT];
    uint32 Textures[TEXTURE_COUNT];
    uint32 VAOs[VAO_COUNT];
    uint32 VBOs[VBO_COUNT];
    uint32 *Indices;
    uint32 *ModelIndices;
    uint32 Num;
    uint32 ModelNum;

    uint32 LightVAO;

    glm::mat4 View;
    glm::mat4 Projection;
};

struct memory
{
    void *PermanentStorage;
    uint64 PermanentStorageSize;
    void *TransientStorage;
    uint64 TransientStorageSize;

    bool32 IsInitialized;
};

struct memory_arena
{
    memory_index Size;
    uint8 *Base;
    memory_index Used;
};

// NOTE(Jovan): Vertex struct, for EPA mostly
struct vertex
{
    glm::vec3 P;
    glm::vec3 SupA;
    glm::vec3 SupB;

    inline bool operator==(const vertex &R)
    {
        return (this->P == R.P);
    }
};

struct simplex
{
    uint32 Count;
    vertex *Vertices;
};

struct edge
{
    uint32 Count;
    vertex A;
    vertex B;
};

struct triangle
{
    uint32 Count;
    vertex A;
    vertex B;
    vertex C;
    glm::vec3 N;
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

struct plane
{
    glm::mat4 Model;
    glm::vec3 Position;
    glm::vec3 Vertices[4];
    real32 Size;
};

struct cube
{
    // TODO(Jovan): Maybe better to keep inverse model matrix?
    glm::mat4 Model;
    glm::vec3 Vertices[8];

    glm::vec3 Position;
    glm::vec3 V;
    glm::vec3 Forces;

    glm::vec3 Angles;
    glm::vec3 W;
    glm::vec3 Torque;
    real32 Size;
    real32 Mass;
    real32 MOI;
};

struct sphere
{
    // TODO(Jovan): Is the model required?
    glm::mat4 Model;
    glm::vec3 Position;
    glm::vec3 V;
    glm::vec3 Forces;

    glm::vec3 Angles;
    glm::vec3 W;
    glm::vec3 Torque;

    real32 Radius;
    real32 Mass;
    real32 MOI;
};

// NOTE(Jovan): Serves as info for a pair of 2 objects
struct contact_pair
{
    // NOTE(Jovan): For identifying the type of collision
    // so we know which types of objects to look at
    collision_type Type;

    // NOTE(Jovan): Delta Lambda Normal
    real32 DLNormal;
    real32 DLNormalSum;

    // NOTE(Jovan): Delta Lambda Tangent
    real32 DLTangent1;
    real32 DLTangent1Sum;

    // NOTE(Jovan): Delta Lambda Tangent
    real32 DLTangent2;
    real32 DLTangent2Sum;

    // NOTE(Jovan): The intersection points
    glm::vec3 PointA;
    glm::vec3 PointB;

    // NOTE(Jovan): Collision normal and tangents
    glm::vec3 N;
    glm::vec3 T1;
    glm::vec3 T2;

    // NOTE(Jovan): Indices of objects
    int32 IndexA;
    int32 IndexB;

    /* glm::vec3 PosA; */
    /* glm::vec3 PosB; */
};

struct sdl_state
{
    // TODO(jovan): Memory arenas as well, maybe into a "world" struct?
    cube Cubes[MAX_SPHERE_COUNT];
    sphere Spheres[MAX_SPHERE_COUNT];
    // NOTE(Jovan): The floor is just a giant cube
    cube Floor;
    sdl_camera Camera;
    light Lights[POINT_LIGHT_COUNT];
    uint32 CubeCount;
    uint32 SphereCount;
    std::vector<contact_pair> Pairs;
};

// NOTE(Jovan): Return value for ODEs and physics functions
// IMPORTANT(Jovan): No need for this as we're not using X anywhere
/* struct phys_return */
/* { */
/*   glm::vec3 X; */
/*   glm::vec3 Y; */
/* }; */

#define SIM_UPDATE_AND_RENDER(name) void name(memory *Memory, sdl_input *Input, sdl_render *Render, real32 dt)
typedef SIM_UPDATE_AND_RENDER(sim_update_and_render);
SIM_UPDATE_AND_RENDER(SimUpdateAndRenderStub)
{
}

#define NANS_H
#endif
