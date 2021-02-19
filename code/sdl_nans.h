enum shader_type{
		 VERTEX = 0,
		 FRAGMENT
};

struct sdl_sim_code
{
  void* SimCodeDynLib;
  time_t DynLibLastWriteTime;
  sim_update_and_render* UpdateAndRender;

  bool32 IsValid;
};
