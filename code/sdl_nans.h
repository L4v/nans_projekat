enum shader_type{
     Vertex = 0,
     Fragment
};

struct memory
{
  void* PermanentStorage;
  uint64 PermanentStorageSize;
  void* TransientStorage;
  uint64 TransientStorageSize;

  bool32 IsInitialized;
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

struct sdl_camera
{
  real32 FOV;
  real32 Pitch;
  real32 Yaw;
  real32 X;
  real32 Y;
  real32 Z;
};

struct sdl_state
{
  glm::vec3 Positions[10];
  sdl_camera Camera;
};
