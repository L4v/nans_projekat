#include <iostream>
#include <GL/glew.h>
#include <SDL2/SDL.h>
#include <stdint.h>
#include <sys/mman.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

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

#define Kibibytes(Value) ((Value) * 1024LL)
#define Mebibytes(Value) (Kibibytes(Value) * 1024LL)
#define Gibibytes(Value) (Mebibytes(Value) * 1024LL)

#define ArrayCount(Array) (sizeof(Array) / sizeof((Array)[0]))

#if SLOW_BUILD
#define Assert(Expression)			\
  if(!(Expression)) {*(int*)0 = 0;}
#else
#define Assert(Expression)
#endif

#define DEFAULT_WINDOW_WIDTH 1024
#define DEFAULT_WINDOW_HEIGHT 768
#define MAX_CUBE_COUNT 64
#define internal static
#define global_variable static
#define local_persist static

global_variable uint64 GlobalPerfCountFrequency;

struct sdl_state
{
  glm::vec3 Positions[64];
  real32 FOV;
  real32 Pitch;
  real32 Yaw;
};

struct memory
{
  void* PermanentStorage;
  uint64 PermanentStorageSize;
  void* TransientStorage;
  uint64 TransientStorageSize;

  bool32 IsInitialized;
};

internal const char*
LoadShader(const char* Path)
{
  char* shaderText = 0;
  int64 length;

  FILE* file = fopen(Path, "rb");
  
  if(file)
    {
      fseek(file, 0, SEEK_END);
      length = ftell(file);
      fseek(file, 0, SEEK_SET);
      shaderText = (char*)malloc(length);
      if(shaderText)
	{
	  fread(shaderText, 1, length, file);
	}
      fclose(file);
    }
  return shaderText;
}

enum shader_type{
     Vertex = 0,
     Fragment
};

  
internal uint32
LoadTexture(char* Path)
{
  int32 Width, Height, NChannels;
  uint8* Data = stbi_load("../res/texture/container.jpg", &Width, &Height, &NChannels, 0);
  GLenum Format;
  uint32 TextureID;
  glGenTextures(1, &TextureID);
  
  if(Data)
    {
      
      if(NChannels == 1)
	{
	  Format = GL_RED;
	}
      else if(NChannels == 3)
	{
	  Format = GL_RGB;
	}
      else if(NChannels == 4)
	{
	  Format = GL_RGBA;
	}
      
      glBindTexture(GL_TEXTURE_2D, TextureID);
      glTexImage2D(GL_TEXTURE_2D, 0, Format, Width, Height, 0, Format, GL_UNSIGNED_BYTE, Data);
      glGenerateMipmap(GL_TEXTURE_2D);
      
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    }
  else
    {
      printf("ERROR::TEXTURE::Failed to load texure!\n");
    }
  stbi_image_free(Data);

  return TextureID;
}

internal void
CheckShaderCompilation(uint32 Shader, shader_type Type)
{
  // NOTE(Jovan): Check if shader compilation failed
  int32 success;
  GLchar infoLog[512];
  glGetShaderiv(Shader, GL_COMPILE_STATUS, &success);
  if(!success)
    {
      glGetShaderInfoLog(Shader, 512, NULL, infoLog);
      std::cout <<
	(Type == Vertex ? "ERROR::SHADER::VERTEX:COMPILATION_FAILED"
	 : "ERROR::SHADER::FRAGMENT:COMPILATION_FAILED")
		<< std::endl << infoLog << std::endl;
    }  
}

internal void
CheckShaderLink(uint32 Program)
{
  int32 success;
  GLchar infoLog[512];
  // NOTE(l4v): Check for program linking failure
  glGetProgramiv(Program, GL_LINK_STATUS, &success);
  if(!success)
    {
      glGetProgramInfoLog(Program, 512, NULL, infoLog);
      std::cout << "ERROR::SHADER_PROGRAM::LINKING_FAILED" << std::endl << infoLog << std::endl;
    }
}

internal void
SetUniformM4(uint32 ID, char* Uniform, const glm::mat4 &Mat4)
{
  glUniformMatrix4fv(glGetUniformLocation(ID, Uniform), 1, GL_FALSE, &Mat4[0][0]);
}

void
SDLWindowResize(int32 Width, int32 Height)
{
  glViewport(0, 0, Width, Height);
}

bool32
SDLHandleEvent(SDL_Event* Event)
{
  bool32 ShouldQuit = 0;
  switch(Event->type)
    {
      // NOTE(Jovan): Quit event
    case SDL_QUIT:
      {
	ShouldQuit = 1;
      }break;

      // NOTE(Jovan): Window events
    case SDL_WINDOWEVENT:
      {
	switch(Event->window.event)
	  {
	  case SDL_WINDOWEVENT_RESIZED:
	    {
	      SDLWindowResize(Event->window.data1, Event->window.data2);
	    }break;
	  }
      }break;

      // NOTE(Jovan): Keyboard events
    case SDL_KEYDOWN:
    case SDL_KEYUP:
      {
	SDL_Keycode KeyCode = Event->key.keysym.sym;
	bool32 IsDown = (Event->key.state == SDL_PRESSED);

	if(!(Event->key.repeat))
	  {
	    if(KeyCode == SDLK_ESCAPE)
	      {
		ShouldQuit = 1;
	      }
	  }
      }break;
      
    }

  return ShouldQuit;
}

inline internal int64
SDLGetWallClock()
{
  int64 Result = SDL_GetPerformanceCounter();
  return Result;
}

inline internal real32
SDLGetSecondsElapsed(int64 Start, int64 End)
{
  real32 Result = (real32)(End - Start) / (real32)GlobalPerfCountFrequency;
  return Result;
}

int main()
{

  // NOTE(Jovan): Memory allocation
  memory SimMemory = {};
  SimMemory.PermanentStorageSize = Mebibytes(64);
  SimMemory.TransientStorageSize = Mebibytes(64);
  uint64 TotalStorageSize = SimMemory.PermanentStorageSize +
    SimMemory.TransientStorageSize;
  SimMemory.PermanentStorage = mmap(0,
				    TotalStorageSize,
				    PROT_READ | PROT_WRITE,
				    MAP_ANON | MAP_PRIVATE,
				    -1,
				    0);
  SimMemory.TransientStorage = ((uint8*)SimMemory.PermanentStorage +
				SimMemory.PermanentStorageSize);
  Assert(SimMemory.PermanentStorage != (void*)-1);
  Assert(SimMemory.TransientStorage != (void*)-1);
  // NOTE(Jovan): SDL stuff
  // --------------------
  SDL_Window* Window = 0;
  SDL_GLContext GLContext = 0;
  GlobalPerfCountFrequency = SDL_GetPerformanceFrequency();
  real32 SimUpdatehz = 60.0f;
  real32 TargetSecondsPerFrame = 1.0f / SimUpdatehz;


  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 0);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

  if(SDL_Init(SDL_INIT_VIDEO) > 0)
    {
      printf("ERROR::SDL::Failed to init!\n");
      return 1;
    }

  uint32 Flags = SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE;
  
  Window = SDL_CreateWindow("NANS Projekat",
			    SDL_WINDOWPOS_CENTERED,
			    SDL_WINDOWPOS_CENTERED,
			    DEFAULT_WINDOW_WIDTH,
			    DEFAULT_WINDOW_HEIGHT,
			    Flags);
  GLContext = SDL_GL_CreateContext(Window);

  if(!GLContext)
    {
      printf("ERROR::SDL::Failed to create GL context!\n");
      return 1;
    }

  // NOTE(Jovan): Get window size and set viewport
  int32 Width, Height;
  SDL_GetWindowSize(Window, &Width, &Height);
  glewInit();
  glViewport(0, 0, Width, Height);  
  
  // NOTE(Jovan): End of SDL stuff
  // -----------------------------

  // NOTE(Jovan): GL modeling and buffering
  // --------------------------------------
  
  const char* CubeFragmentShaderSource = LoadShader("../shaders/cube.fs");
  const char* CubeVertexShaderSource = LoadShader("../shaders/cube.vs");

  // real32 RectVertices[] = {
  // 			   // X  |  Y   |  Z  | Tex coords 
  // 			   -0.5f, -0.5f, 0.0f, 0.0f, 0.0f, // LL
  // 			   -0.5f,  0.5f, 0.0f, 0.0f, 1.0f, // UL
  // 			    0.5f, -0.5f, 0.0f, 1.0f, 0.0f, // LR
  // 			    0.5f,  0.5f, 0.0f, 1.0f, 1.0f  // UR
  // };

  real32 RectVertices[] = {
			   // X  |  Y   |  Z  | Tex coords
			   -0.5f, -0.5f, -0.5f,  0.0f, 0.0f,
			   0.5f, -0.5f, -0.5f,  1.0f, 0.0f,
			   0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
			   0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
			   -0.5f,  0.5f, -0.5f,  0.0f, 1.0f,
			   -0.5f, -0.5f, -0.5f,  0.0f, 0.0f,

			   -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
			   0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
			   0.5f,  0.5f,  0.5f,  1.0f, 1.0f,
			   0.5f,  0.5f,  0.5f,  1.0f, 1.0f,
			   -0.5f,  0.5f,  0.5f,  0.0f, 1.0f,
			   -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,

			   -0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
			   -0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
			   -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
			   -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
			   -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
			   -0.5f,  0.5f,  0.5f,  1.0f, 0.0f,

			   0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
			   0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
			   0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
			   0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
			   0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
			   0.5f,  0.5f,  0.5f,  1.0f, 0.0f,

			   -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
			   0.5f, -0.5f, -0.5f,  1.0f, 1.0f,
			   0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
			   0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
			   -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
			   -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,

			   -0.5f,  0.5f, -0.5f,  0.0f, 1.0f,
			   0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
			   0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
			   0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
			   -0.5f,  0.5f,  0.5f,  0.0f, 0.0f,
			   -0.5f,  0.5f, -0.5f,  0.0f, 1.0f
  };

  uint32 RectIndices[] = {
			  0, 1, 3,
			  0, 2, 3
  };
  
  // NOTE(Jovan): Creating shaders and shader programs
  uint32 CubeVertexShader, CubeFragmentShader;
  CubeVertexShader = glCreateShader(GL_VERTEX_SHADER);
  CubeFragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

  glShaderSource(CubeVertexShader, 1, &CubeVertexShaderSource, 0);
  glShaderSource(CubeFragmentShader, 1, &CubeFragmentShaderSource, 0);

  glCompileShader(CubeVertexShader);
  CheckShaderCompilation(CubeVertexShader, Vertex);
  glCompileShader(CubeFragmentShader);
  CheckShaderCompilation(CubeVertexShader, Fragment);

  uint32 CubeShaderProgram;
  CubeShaderProgram = glCreateProgram();
  glAttachShader(CubeShaderProgram, CubeVertexShader);
  glAttachShader(CubeShaderProgram, CubeFragmentShader);
  glLinkProgram(CubeShaderProgram);
  CheckShaderLink(CubeShaderProgram);

  glDeleteShader(CubeVertexShader);
  glDeleteShader(CubeFragmentShader);
  
  // NOTE(Jovan): VAO, EBO, VBO
  uint32 CubeVAO, CubeVBO, CubeEBO;
  
  glGenVertexArrays(1, &CubeVAO);
  glBindVertexArray(CubeVAO);
  glGenBuffers(1, &CubeVBO);
  glGenBuffers(1, &CubeEBO);

  glBindBuffer(GL_ARRAY_BUFFER, CubeVBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(RectVertices), RectVertices,
	       GL_DYNAMIC_DRAW);
  
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, CubeEBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(RectIndices), RectIndices,
	       GL_DYNAMIC_DRAW);
  
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(real32),
			(void*) 0);
  glEnableVertexAttribArray(0);
  
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(real32),
			(void*) (3 * sizeof(real32)));
  glEnableVertexAttribArray(1);
  glBindVertexArray(0);


  // NOTE(Jovan): Textures
  uint32 CubeTexture;

  CubeTexture = LoadTexture("../res/texture/container.jpg");

  
  // NOTE(Jovan): End of GL modeling and buffering
  // ---------------------------------------------

  // NOTE(Jovan): Coordinate systems
  // -------------------------------
  glm::mat4 Projection = glm::perspective(glm::radians(45.0f),
					  (real32)Width / (real32)Height,
					  0.1f,
					  100.0f);
  glm::mat4 Model = glm::mat4(1.0f);
  glm::mat4 View = glm::mat4(1.0f);
  View = glm::translate(View, glm::vec3(0.0f, 0.0f, -3.0f));
  glm::vec3 CameraPos = glm::vec3(0.0f, 0.0f, 3.0f);
  glm::vec3 CameraTarget = glm::vec3(0.0f, 0.0f, 0.0f);
  glm::vec3 CameraDirection = glm::normalize(CameraPos - CameraTarget);
  glm::vec3 Up = glm::vec3(0.0f, 1.0f, 0.0f);
  glm::vec3 CameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
  glm::vec3 CameraRight = glm::normalize(glm::cross(Up, CameraDirection));
  glm::vec3 CameraUp = glm::cross(CameraDirection, CameraRight);
  View = glm::lookAt(CameraPos, CameraPos + CameraFront, CameraUp);

  // NOTE(Jovan): End of coordinate systems
  // --------------------------------------

  glEnable(GL_DEPTH_TEST);
  
  // NOTE(Jovan): Main loop
  bool32 Running = 1;
  uint64 LastCounter = SDLGetWallClock();
  real32 dt = 0.0f;
  real32 MouseSensitivity = 0.5f;
  // NOTE(Jovan): Grab mouse
  SDL_SetRelativeMouseMode(SDL_TRUE);
  while(Running)
    {

      // TODO(Jovan): Move to sim update / render
      sdl_state* SimState = (sdl_state*)SimMemory.PermanentStorage;
      if(!SimMemory.IsInitialized)
	{
	  SimState->FOV = 45.0f; 
	  SimState->Pitch = 0.0f;
	  SimState->Yaw = -90.0f;
	  for(int CubeIndex = 0;
	      CubeIndex < 3;
	      ++CubeIndex)
	    {
	      SimState->Positions[CubeIndex] = glm::vec3(0.2*CubeIndex,
							 0.3*CubeIndex,
							 0.4*CubeIndex);
	    }
	  SimMemory.IsInitialized = 1;
	}
      
      SDL_Event Event;
      while(SDL_PollEvent(&Event))
	{
	  // TODO(Jovan): Move this somewhere else maybe?
	  if(Event.type == SDL_MOUSEMOTION)
	    {
	      SimState->Yaw += Event.motion.xrel * MouseSensitivity;
	      SimState->Pitch += -Event.motion.yrel * MouseSensitivity;
	      if(SimState->Pitch > 89.0f)
		{
		  SimState->Pitch = 89.0f;
		}
	      if(SimState->Pitch < -89.0f)
		{
		  SimState->Pitch = -89.0f;
		}
	    }
	  // NOTE(Jovan): Check for exit
	  if(SDLHandleEvent(&Event))
	    {
	      Running = 0;
	    }
	}

      glm::vec3 Front;
      Front.x = cos(glm::radians(SimState->Yaw)) * cos(glm::radians(SimState->Pitch));
      Front.y = sin(glm::radians(SimState->Pitch));
      Front.z = sin(glm::radians(SimState->Yaw)) * cos(glm::radians(SimState->Pitch));
      CameraFront = glm::normalize(Front);
      View = glm::lookAt(CameraPos, CameraPos + CameraFront, CameraUp);
      
      // NOTE(Jovan): Work timing
      int64 WorkCounter = SDLGetWallClock();
      real32 WorkSecondsElapsed = SDLGetSecondsElapsed(LastCounter, WorkCounter);
      real32 SecondsElapsedForFrame = WorkSecondsElapsed;
      if(SecondsElapsedForFrame < TargetSecondsPerFrame)
	{
	  while(SecondsElapsedForFrame < TargetSecondsPerFrame)
	    {
	      // TODO(Jovan): Check sleep granulaity on linux
	      uint32 MSToSleepFor = (uint32)(1000.0f * (TargetSecondsPerFrame -
							SecondsElapsedForFrame));
	      SDL_Delay(MSToSleepFor);
	      SecondsElapsedForFrame = SDLGetSecondsElapsed(LastCounter, SDLGetWallClock());
	    }
	}
      else
	{
	  // TODO(Jovan): Missed frame rate???
	}

      glClearColor(0.8f, 0.0f, 0.8f, 1.0f);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      // NOTE(Jovan): Cube drawing
      // -------------------------
      glUseProgram(CubeShaderProgram);
      glBindTexture(GL_TEXTURE_2D, CubeTexture);

      SetUniformM4(CubeShaderProgram, "View", View);
      SetUniformM4(CubeShaderProgram, "Projection", Projection);

      for(int CubeIndex = 0;
	  CubeIndex < 10;
	  ++CubeIndex)
	{
	  int32 CubeSize = 1.0f;
	  real32 Angle = 50.0f;
	  Model = glm::mat4(1.0f);
	  Model = glm::scale(Model, glm::vec3(CubeSize, CubeSize, CubeSize));
	  Model = glm::translate(Model, glm::vec3(0.0, 0.0, 0.0));//SimState->Positions[CubeIndex]);
	  Model = glm::rotate(Model, glm::radians(Angle), glm::vec3(1.0f, 1.0f, 1.0f));
	  SetUniformM4(CubeShaderProgram, "Model", Model);
	  glBindVertexArray(CubeVAO);
	  glDrawArrays(GL_TRIANGLES, 0, 36);
	}

      // NOTE(Jovan): End cube drawing
      // -----------------------------
      
      glBindVertexArray(0);

      // NOTE(Jovan): Swap buffers
      SDL_GL_SwapWindow(Window);

      // NOTE(Jovan): Timing
      int64 EndCounter = SDLGetWallClock();
      real32 SPerFrame = SDLGetSecondsElapsed(LastCounter, EndCounter);
      real32 FPS = 1.0f / SPerFrame;
      real32 MSPerFrame = 1000.0f * SDLGetSecondsElapsed(LastCounter, EndCounter);
      dt = SPerFrame;
      LastCounter = EndCounter;

      printf("MSPerFrame = %f, FPS = %f\n", MSPerFrame, FPS);
      
    }
  
  return 0;
}
