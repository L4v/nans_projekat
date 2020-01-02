internal void
SetUniformM4(uint32 ID, char* Uniform, const glm::mat4 &Mat4)
{
  glUniformMatrix4fv(glGetUniformLocation(ID, Uniform), 1, GL_FALSE, &Mat4[0][0]);
}

internal void
SetUniformV3(uint32 ID, char* Uniform, const glm::vec3 &Vec3)
{
  glUniform3fv(glGetUniformLocation(ID, Uniform), 1, &Vec3[0]);
}

internal void
InitializeArena(memory_arena* Arena, memory_index Size, uint8* Base)
{
  Arena->Size = Size;
  Arena->Base = Base;
  Arena->Used = 0;
}

#define PushStruct(Arena, Type) (Type*)PushSize_(Arena, sizeof(Type))
#define PushArray(Arena, Count, Type) (Type*)PushSize_(Arena, (Count) * sizeof(Type))
internal void*
PushSize_(memory_arena* Arena, memory_index Size)
{
  Assert((Arena->Used + Size) <= Arena->Size);
  void* Result = Arena->Base + Arena->Used;
  Arena->Used += Size;

  return Result;
}
