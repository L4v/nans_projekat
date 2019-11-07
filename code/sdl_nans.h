enum shader_type{
     Vertex = 0,
     Fragment
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
  
  glm::vec3 Position;
  glm::vec3 Target;
  glm::vec3 Direction;
  
  glm::vec3 Up;
  glm::vec3 Front;
  glm::vec3 Right;
};

struct memory
{
  void* PermanentStorage;
  uint64 PermanentStorageSize;
  void* TransientStorage;
  uint64 TransientStorageSize;

  bool32 IsInitialized;
};

struct sdl_render
{
  uint32 Shaders[1];
  uint32 Textures[1];
  uint32 VAOs[1];
};

struct sdl_state
{
  glm::vec3 Positions[10];
  sdl_camera Camera;
};

