#include "nans.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <SDL2/SDL.h>
#include <sys/mman.h>
#include <dlfcn.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>

#include "sdl_nans.h"

global_variable uint64 GlobalPerfCountFrequency;

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
      if(Type == Vertex)
	{
	  printf("ERROR::SHADER::VERTEX:COMPILATION_FAILED\n");
	}
      else
	{
	  printf("ERROR::SHADER::FRAGMENT:COMPILATION_FAILED\n");
	}
      printf("%s\n", infoLog);
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
      printf("ERROR::SHADER_PROGRAM::LINKING_FAILED\n%s\n", infoLog);
    }
}

internal void
SDLProcessSimKeyboardButton(sdl_button_state* NewState, bool32 IsDown)
{
  Assert(NewState->EndedDown != IsDown);
  NewState->EndedDown = IsDown;
  ++NewState->HalfTransitionCount;
}

void
SDLWindowResize(int32 Width, int32 Height)
{
  glViewport(0, 0, Width, Height);
}

internal bool32
SDLHandleEvent(SDL_Event* Event, sdl_input* Input) 
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
    case SDL_MOUSEMOTION:
      {
	Input->MouseController.X = Event->motion.x;
	Input->MouseController.Y = Event->motion.y;
	Input->MouseController.XRel = Event->motion.xrel;
	Input->MouseController.YRel = Event->motion.yrel;
      }break;

      // NOTE(Jovan): Keyboard events
    case SDL_KEYDOWN:
    case SDL_KEYUP:
      {
	SDL_Keycode KeyCode = Event->key.keysym.sym;
	bool32 IsDown = (Event->key.state == SDL_PRESSED);

	if(Event->key.state == SDL_RELEASED)
	  {

	  }
	else if(Event->key.repeat)
	  {
	    
	  }

	if(!(Event->key.repeat))
	  {
	    if(KeyCode == SDLK_w)
	      {
		SDLProcessSimKeyboardButton(&Input->KeyboardController.MoveForward, IsDown);
	      }
	    else if(KeyCode == SDLK_a)
	      {
		SDLProcessSimKeyboardButton(&Input->KeyboardController.MoveLeft, IsDown);
	      }
	    else if(KeyCode == SDLK_d)
	      {
		SDLProcessSimKeyboardButton(&Input->KeyboardController.MoveRight, IsDown);
	      }
	    else if(KeyCode == SDLK_s)
	      {
		SDLProcessSimKeyboardButton(&Input->KeyboardController.MoveBack, IsDown);
	      }
	    else if(KeyCode == SDLK_SPACE)
	      {
		SDLProcessSimKeyboardButton(&Input->KeyboardController.ShootAction, IsDown);
	      }
	    else if(KeyCode == SDLK_ESCAPE)
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

internal inline time_t
SDLGetLastWriteTime(char* Filename)
{
  struct stat FileInfo = {};
  time_t Result = 0;
  if(stat(Filename, &FileInfo) != -1)
    {
      Result = FileInfo.st_mtim.tv_sec;
    }

  return Result;
}

internal sdl_sim_code
SDLLoadSimCode(char* DynLibName)
{
  sdl_sim_code Result = {};
  Result.SimCodeDynLib = dlopen(DynLibName, RTLD_NOW | RTLD_GLOBAL);
  Result.DynLibLastWriteTime = SDLGetLastWriteTime(DynLibName);

  if(Result.SimCodeDynLib)
    {
      Result.UpdateAndRender = (sim_update_and_render*)
	dlsym(Result.SimCodeDynLib, "SimUpdateAndRender");

      Result.IsValid = (Result.UpdateAndRender != 0);
    }
  
  if(!Result.IsValid)
    {
      printf("WARNING::SimCode::Not loaded properly, using stub.\n");
      Result.UpdateAndRender = SimUpdateAndRenderStub;
    }
  
  return Result;
}

internal void
SDLUnloadSimCode(sdl_sim_code* SimCode)
{
  if(dlclose(SimCode->SimCodeDynLib))
    {
      printf("ERROR::SimCode::Not unloaded properly.\n");
      char* Error = dlerror();
      printf("%s\n", Error);
    }
  SimCode->SimCodeDynLib = 0;
  SimCode->IsValid = 0;
  SimCode->UpdateAndRender = SimUpdateAndRenderStub;
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
  const char* SphereVSSource = LoadShader("../shaders/sphere.vs");
  const char* SphereFSSource = LoadShader("../shaders/sphere.fs");

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

  int32 Stacks = 20;
  int32 Slices = 20;
  
  real32 Vertices[(Stacks + 1) * (Slices + 1) * 3];
  uint32 Indices[(Slices * Stacks + Slices) * 6];

  int32 Index = 0;
  for(int i = 0; i <= Stacks; i++)
    {
      real32 V = (real32)i / (real32)Stacks;
      real32 Phi = V * Pi32;
      
      for(int j = 0; j <= Slices; j++)
	{
	  real32 U = (real32)j / (real32)Slices;
	  real32 Theta = U * (Pi32 * 2);

	  real32 x = cos(Theta) * sin(Phi);
	  real32 y = cos(Phi);
	  real32 z = sin(Theta) * sin(Phi);

	  Vertices[Index++] = x;
	  Vertices[Index++] = y;
	  Vertices[Index++] = z;
	}
    }
  Index = 0;
  for(int32 i = 0; i < Slices * Stacks + Slices; i++)
    {
      Indices[Index++] = i;
      Indices[Index++] = i + Slices + 1;
      Indices[Index++] = i + Slices;

      Indices[Index++] = i + Slices + 1;
      Indices[Index++] = i;
      Indices[Index++] = i + 1;
    }
  
  // NOTE(Jovan): Creating shaders and shader programs
  uint32 CubeVertexShader, CubeFragmentShader,
    SphereVS, SphereFS;
  
  CubeVertexShader = glCreateShader(GL_VERTEX_SHADER);
  CubeFragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
  SphereVS  = glCreateShader(GL_VERTEX_SHADER);
  SphereFS = glCreateShader(GL_FRAGMENT_SHADER);

  glShaderSource(CubeVertexShader, 1, &CubeVertexShaderSource, 0);
  glShaderSource(CubeFragmentShader, 1, &CubeFragmentShaderSource, 0);
  glShaderSource(SphereVS, 1, &SphereVSSource, 0);
  glShaderSource(SphereFS, 1, &SphereFSSource, 0);
  glCompileShader(SphereVS);
  CheckShaderCompilation(SphereVS, Vertex);
  glCompileShader(SphereFS);
  CheckShaderCompilation(SphereFS, Fragment);

  glCompileShader(CubeVertexShader);
  CheckShaderCompilation(CubeVertexShader, Vertex);
  glCompileShader(CubeFragmentShader);
  CheckShaderCompilation(CubeVertexShader, Fragment);

  uint32 CubeShaderProgram, SphereShaderProgram;
  CubeShaderProgram = glCreateProgram();
  SphereShaderProgram = glCreateProgram();
  
  glAttachShader(CubeShaderProgram, CubeVertexShader);
  glAttachShader(CubeShaderProgram, CubeFragmentShader);
  glAttachShader(SphereShaderProgram, SphereVS);
  glAttachShader(SphereShaderProgram, SphereFS);
  
  glLinkProgram(CubeShaderProgram);
  CheckShaderLink(CubeShaderProgram);
  
  glLinkProgram(SphereShaderProgram);
  CheckShaderLink(SphereShaderProgram);

  glDeleteShader(CubeVertexShader);
  glDeleteShader(CubeFragmentShader);
  glDeleteShader(SphereVS);
  glDeleteShader(SphereFS);
  
  // NOTE(Jovan): VAO, EBO, VBO
  uint32 CubeVAO, CubeVBO, CubeEBO,
    SphereVAO, SphereVBO, SphereEBO;

  // NOTE(Jovan): Cube data
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

  // NOTE(Jovan): Sphere data
  glGenVertexArrays(1, &SphereVAO);
  glBindVertexArray(SphereVAO);
  
  glGenBuffers(1, &SphereVBO);
  glBindBuffer(GL_ARRAY_BUFFER, SphereVBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(Vertices), Vertices, GL_DYNAMIC_DRAW);

  glGenBuffers(1, &SphereEBO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, SphereEBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(Indices), Indices, GL_DYNAMIC_DRAW);
  
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3  * sizeof(real32), (void*)0);
  glEnableVertexAttribArray(0);
  
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
  // NOTE(Jovan): Textures
  uint32 CubeTexture;

  CubeTexture = LoadTexture("../res/texture/container.jpg");
  
  // NOTE(Jovan): End of GL modeling and buffering
  // ---------------------------------------------

  glEnable(GL_DEPTH_TEST);

  // NOTE(Jovan): Main loop
  bool32 Running = 1;
  uint64 LastCounter = SDLGetWallClock();
  real32 dt = 0.0f;
  sdl_input SimInput[2];
  sdl_input* NewInput = &SimInput[0];
  sdl_input* OldInput = &SimInput[1];
  *NewInput = {};
  *OldInput = {};
  
  // NOTE(Jovan): Grab mouse
  SDL_SetRelativeMouseMode(SDL_TRUE);

  sdl_sim_code Sim = SDLLoadSimCode("nans.so");

  // NOTE(JOVAN): SPHERE
  // -------------------


  // NOTE(JOVAN): SPHERE END
  // ----------------------

  
  // NOTE(Jovan): Renderer
  sdl_render Render = {};
  Render.Shaders[0] = CubeShaderProgram;
  Render.Textures[0] = CubeTexture;
  Render.VAOs[0] = CubeVAO;
  
  Render.Shaders[1] = SphereShaderProgram;
  Render.VAOs[1] = SphereVAO;
  Render.VBOs[0] = SphereVBO;
  Render.Indices = Indices;
  Render.Num = ArrayCount(Indices);
  
  while(Running)
    {

      time_t NewDynLibWriteTime = SDLGetLastWriteTime("nans.so");
      if((NewDynLibWriteTime != Sim.DynLibLastWriteTime))
	{
	  printf("nans.so difference: %ld\n", NewDynLibWriteTime - Sim.DynLibLastWriteTime);
	  printf("Code changed!\n");
	  SDLUnloadSimCode(&Sim);
	  SDL_Delay(100);
	  Sim = SDLLoadSimCode("nans.so");
	}

      sdl_keyboard* OldKeyboardController = &OldInput->KeyboardController;
      sdl_keyboard* NewKeyboardController = &NewInput->KeyboardController;
      *NewKeyboardController = {};

      sdl_mouse* OldMouseController = &OldInput->MouseController;
      sdl_mouse* NewMouseController = &NewInput->MouseController;
      *NewMouseController = {};
      NewMouseController->Sensitivity = 0.5f;

      for(uint32 ButtonIndex = 0;
	  ButtonIndex < ArrayCount(NewKeyboardController->Buttons);
	  ++ButtonIndex)
	{
	  NewKeyboardController->Buttons[ButtonIndex].EndedDown =
	    OldKeyboardController->Buttons[ButtonIndex].EndedDown;
	}

      NewInput->MouseController.X = OldInput->MouseController.X;
      NewInput->MouseController.Y = OldInput->MouseController.Y;
      NewInput->MouseController.XRel = OldInput->MouseController.XRel;
      NewInput->MouseController.YRel = OldInput->MouseController.YRel;
      
      SDL_Event Event;
      while(SDL_PollEvent(&Event))
	{
	  // NOTE(Jovan): Check for exit
	  if(SDLHandleEvent(&Event, NewInput))
	    {
	      Running = 0;
	    }
	}
      
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
      
      glClearColor(0.8f, 0.8f, 0.8f, 1.0f);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      Sim.UpdateAndRender(&SimMemory, NewInput, &Render, dt);
      glBindVertexArray(0);

      // NOTE(Jovan): SPHERE
      // -------------------

      // NOTE(Jovan): END SPHERE
      // -----------------------

      
      // NOTE(Jovan): Swap buffers
      SDL_GL_SwapWindow(Window);

      // NOTE(Jovan): Timing
      int64 EndCounter = SDLGetWallClock();
      real32 SPerFrame = SDLGetSecondsElapsed(LastCounter, EndCounter);
      real32 FPS = 1.0f / SPerFrame;
      real32 MSPerFrame = 1000.0f * SDLGetSecondsElapsed(LastCounter, EndCounter);
      dt = SPerFrame;
      LastCounter = EndCounter;
      char Buffer[256];
      sprintf(Buffer, "MSPerFrame = %f, FPS = %f", MSPerFrame, FPS);
      SDL_SetWindowTitle(Window, Buffer );

      sdl_input* Temp = NewInput;
      NewInput = OldInput;
      OldInput = Temp;
    }
  
  return 0;
}
