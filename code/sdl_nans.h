enum shader_type{
		 Vertex = 0,
		 Fragment
};

struct sdl_sim_code
{
  void* SimCodeDynLib;
  time_t DynLibLastWriteTime;
  sim_update_and_render* UpdateAndRender;

  bool32 IsValid;
};
