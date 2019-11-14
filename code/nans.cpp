#include "nans.h"
internal void
SetUniformM4(uint32 ID, char* Uniform, const glm::mat4 &Mat4)
{
  glUniformMatrix4fv(glGetUniformLocation(ID, Uniform), 1, GL_FALSE, &Mat4[0][0]);
}

internal void
UpdateCamera(sdl_state* State, sdl_input* Input)
{

  State->Camera.Yaw += Input->MouseController.XRel * Input->MouseController.Sensitivity;
  State->Camera.Pitch += -Input->MouseController.YRel * Input->MouseController.Sensitivity;
  Input->MouseController.XRel = 0.0f;
  Input->MouseController.YRel = 0.0f;
  if(State->Camera.Pitch > 89.0f)
    {
      State->Camera.Pitch = 89.0f;
    }
  if(State->Camera.Pitch < -89.0f)
    {
      State->Camera.Pitch = -89.0f;
    }
}

internal phys_return
Euler(phys_return (*F)(glm::vec3, glm::vec3, real32), real32 dt,
      phys_return Y0, glm::vec3 SummedForce, real32 Mass)
{
  phys_return Result = {};
  phys_return Tmp = F(Y0.Y, SummedForce, Mass);
  Result.X = Y0.X + dt * Tmp.X;
  Result.Y = Y0.Y + dt* Tmp.Y;
  return Result;
}

internal phys_return
MovementFunction(glm::vec3 Velocity, glm::vec3 SummedForces, real32 Mass)
{
  phys_return Result = {};
  Result.X = Velocity;
  Result.Y = (real32)(1.0 / Mass) * (SummedForces - GLOBAL_FRICTION * Velocity);
  return Result;
}

internal void
RotationFunction()
{
  // TODO(Jovan): Implement
  // ubrzanje = 1/m (...);
}

internal void
CubeAddForce(sdl_state* State, int32 CubeIndex, glm::vec3 Force)
{
  State->Cubes[CubeIndex].Forces += Force;
};

internal void
CubeClearForces(sdl_state* State, int32 CubeIndex)
{
  State->Cubes[CubeIndex].Forces = glm::vec3(0.0);
};

internal void
ResolveCollision()
{
  // TODO(Jovan): Implement
}
internal void
RemoveVertex(sdl_state* State, uint32 VertexIndex)
{
  Assert(VertexIndex < State->VertexCount);
  for(uint32 i = VertexIndex;
      i < State->VertexCount - 1;
      ++i)
    {
      State->Vertices[i] = State->Vertices[i + 1];
    }
  State->VertexCount--;
}
internal inline void
PushVertex(sdl_state* State, glm::vec3 Vertex)
{
  Assert(State->VertexCount < 4);
  State->Vertices[State->VertexCount++] = Vertex;
}

internal inline void
ResetVertices(sdl_state* State)
{
  State->VertexCount = 0;
}

internal void
UpdateVertices(sdl_state* State, int32 CubeIndex)
{
  real32 Size = State->Cubes[CubeIndex].Size;
  State->Cubes[CubeIndex].Vertices[0] = State->Cubes[CubeIndex].Position +
    Size / 2.0f * glm::vec3(1.0, 1.0, 1.0);
  State->Cubes[CubeIndex].Vertices[1] = State->Cubes[CubeIndex].Position +
    Size / 2.0f * glm::vec3(1.0, 1.0, -1.0);
  State->Cubes[CubeIndex].Vertices[2] = State->Cubes[CubeIndex].Position +
    Size / 2.0f * glm::vec3(-1.0, 1.0, 1.0);
  State->Cubes[CubeIndex].Vertices[3] = State->Cubes[CubeIndex].Position +
    Size / 2.0f * glm::vec3(-1.0, 1.0, -1.0);
  State->Cubes[CubeIndex].Vertices[4] = State->Cubes[CubeIndex].Position +
    Size / 2.0f * glm::vec3(1.0, -1.0, 1.0);
  State->Cubes[CubeIndex].Vertices[5] = State->Cubes[CubeIndex].Position +
    Size / 2.0f * glm::vec3(1.0, -1.0, -1.0);
  State->Cubes[CubeIndex].Vertices[6] = State->Cubes[CubeIndex].Position +
    Size / 2.0f * glm::vec3(-1.0, -1.0, 1.0);
  State->Cubes[CubeIndex].Vertices[7] = State->Cubes[CubeIndex].Position +
    Size / 2.0f * glm::vec3(-1.0, -1.0, -1.0);
}

internal glm::vec3
GetCubeSupport(sdl_state* State, int32 CubeIndex, glm::vec3 Direction)
{
  real32 MaxDistance = -FLT_MAX;
  glm::vec3 Result = glm::vec3(0.0f);

  // NOTE(Jovan) 8 vertices in a cube
  for(uint32 VertexIndex = 0;
      VertexIndex < 8;
      ++VertexIndex)
    {
      glm::vec3 CurrVertex = State->Cubes[CubeIndex].Vertices[VertexIndex];
      real32 Distance = glm::dot(CurrVertex, Direction);
      if(Distance > MaxDistance)
	{
	  MaxDistance = Distance;
	  Result = CurrVertex;
	}
    }

  return Result;
}

internal glm::vec3
CalculateSupport(sdl_state* State, int32 AIndex, int32 BIndex, glm::vec3 Direction, collision_type Type)
{
  glm::vec3 Result = glm::vec3(0.0f);
  switch(Type)
    {
    case CC:
      {
	glm::vec3 SupportA = GetCubeSupport(State, AIndex, Direction);
	glm::vec3 SupportB = GetCubeSupport(State, BIndex, -1.0f * Direction);
	Result = SupportA - SupportB;
      }break;
    case CS:
      {
	// TODO(Jovan): Implement
      }break;
    case SS:
      {
	// TODO(Jovan): Implement
      }break;
    }

  return Result;
}

internal inline glm::vec3
DoubleCross(glm::vec3 A, glm::vec3 B)
{
  glm::vec3 Result = glm::cross(glm::cross(A, B), A);
  return Result;
}

internal bool32
AddSupport(sdl_state* State, int32 AIndex, int32 BIndex, glm::vec3 Direction)
{
  bool32 Result = 0;
  glm::vec3 NewVertex = CalculateSupport(State, AIndex, BIndex, Direction, CC);
  PushVertex(State, NewVertex);
  if(glm::dot(Direction, NewVertex) >= 0)
    {
      Result = 1;
    }
  else
    {
      Result = 0;
    }
  return Result;
}

internal evolve_result
EvolveSimplex(sdl_state* State, int32 AIndex, int32 BIndex)
{
  // TODO(Jovan): For now, under the assumption that it's 2 cubes
  evolve_result Result = NoIntersection;
  glm::vec3 Direction = glm::vec3(0.0f);
  cube* ShapeA = &State->Cubes[AIndex];
  cube* ShapeB = &State->Cubes[BIndex];
  uint32 NoVertices = State->VertexCount;
  switch(NoVertices)
    {
    case 0:
      {
	// NOTE(Jovan): Best presumption of dir is the relative direction
	// of the 2 shapes
	Direction = ShapeB->Position - ShapeA->Position;
      }break;
    case 1:
      {
	// NOTE(Jovan): Flip the direction to point to the origin
	Direction *= -1;
      }break;
    case 2:
      {
	// NOTE(Jovan): Form line from first 2 vertices
	glm::vec3 AB = State->Vertices[1] - State->Vertices[0];
	// NOTE(Jovan): Form line from origin to A
	glm::vec3 A0 = -1.0f * State->Vertices[0];

	// NOTE(Jovan): Direction perpendicular to AB in the direction
	// of the origin
	Direction = DoubleCross(AB, A0);
      }break;
    case 3:
      {
	glm::vec3 AC = State->Vertices[2] - State->Vertices[0];
	glm::vec3 AB = State->Vertices[1] - State->Vertices[0];
	Direction = glm::cross(AC, AB);

	// NOTE(Jovan): Ensure that Direction points to the origin
	glm::vec3 A0 = -1.0f * State->Vertices[0];
	if(glm::dot(Direction, A0) < 0)
	  {
	    Direction *= -1.0f;
	  }
      }break;
    case 4:
      {
	// NOTE(Jovan): 3 edges of interest
	glm::vec3 DA = State->Vertices[3] - State->Vertices[0];
	glm::vec3 DB = State->Vertices[3] - State->Vertices[1];
	glm::vec3 DC = State->Vertices[3] - State->Vertices[2];

	// NOTE(Jovan): Dir to the origin
	glm::vec3 D0 = -1.0f * State->Vertices[3];

	// NOTE(Jovan): Check triangles ABD, BCD, CAD
	glm::vec3 ABDNorm = glm::cross(DA, DB);
	glm::vec3 BCDNorm = glm::cross(DB, DC);
	glm::vec3 CADNorm = glm::cross(DC, DA);
	if(glm::dot(ABDNorm, D0) > 0)
	  {
	    // NOTE(Jovan): Origin outside of ABD -> remove C
	    RemoveVertex(State, 3);
	    Direction = ABDNorm;
	  }
	else if(glm::dot(BCDNorm, D0) > 0)
	  {
	    // NOTE(Jovan): Origin outside of BCD -> remove A
	    RemoveVertex(State, 0);
	    Direction = BCDNorm;
	  }
	else if(glm::dot(CADNorm, D0) > 0)
	  {
	    // NOTE(Jovan): Origin outside of CAD -> remove B
	    RemoveVertex(State, 1);
	    Direction = CADNorm;
	  }
	else
	  {
	    // NOTE(Jovan): Origin is inside of the tetrahedron
	    Result = FoundIntersection;
	    return Result;
	  }
	
      }break;
    default:
      {
	printf("ERROR::There can't be more than 4 vertices\n");
      }
    }
  if(AddSupport(State, AIndex, BIndex, Direction))
    {
      Result = StillEvolving;
    }
  else
    {
      Result = NoIntersection;
    }
  return Result;
}

internal bool32
DetectCollisions(sdl_state* State, int32 AIndex, int32 BIndex, collision_type Type)
{
  bool32 Result = 0;
  evolve_result EvolutionResult = StillEvolving;

  // TODO(Jovan): Possibly infinite while loop
  while(EvolutionResult == StillEvolving && State->GJKIteration++ <= MAX_GJK_ITERATIONS)
    {
      EvolutionResult = EvolveSimplex(State, AIndex, BIndex);
    }
  if(EvolutionResult == FoundIntersection)
    {
      Result = 1;
    }
  if(EvolutionResult == NoIntersection)
    {
      Result = 0;
    }
  State->GJKIteration = 0;
  return Result;
}

extern "C" SIM_UPDATE_AND_RENDER(SimUpdateAndRender)
{
  sdl_state* SimState = (sdl_state*)Memory->PermanentStorage;

  // NOTE(Jovan): Init
  // -----------------
  if(!Memory->IsInitialized)
    {

      // NOTE(Jovan): Init random seed
      srand(time(0));

      // NOTE(Jovan): GJK init stuff
      SimState->GJKIteration = 0;
      for(uint32 VertexIndex = 0;
	  VertexIndex < ArrayCount(SimState->Vertices);
	  ++VertexIndex)
	{
	  SimState->Vertices[VertexIndex] = glm::vec3(0.0f);
	}
      SimState->VertexCount = 0;

      
      // NOTE(Jovan): Camera init
      SimState->Camera.FOV = 45.0f; 
      SimState->Camera.Pitch = 0.0f;
      SimState->Camera.Yaw = -90.0f;
      SimState->Camera.Position = glm::vec3(0.0f, 0.0f, 3.0f);
      glm::vec3 Up = glm::vec3(0.0f, 1.0f, 0.0f);
      SimState->Camera.Target = glm::vec3(0.0f, 0.0f, 0.0f);
      SimState->Camera.Direction = glm::normalize(SimState->Camera.Position - SimState->Camera.Target);
      SimState->Camera.Front = glm::vec3(0.0f, 0.0f, -1.0f);
      SimState->Camera.Right = glm::normalize(glm::cross(Up, SimState->Camera.Direction));
      SimState->Camera.Up = glm::cross(SimState->Camera.Direction, SimState->Camera.Right);
      
      // NOTE(Jovan): Cube init
      SimState->CubeCount = 2;
      SimState->Cubes[0].Position = glm::vec3(0.0,
					      2.0,
					      0.0);
      SimState->Cubes[0].Velocity = glm::vec3(0.0);
      SimState->Cubes[0].Forces = glm::vec3(0.0);
      SimState->Cubes[0].XAngle = 0.0f;
      SimState->Cubes[0].YAngle = 0.0f;
      SimState->Cubes[0].ZAngle = 0.0f;
      SimState->Cubes[0].Size = 1.0f;
      SimState->Cubes[0].Mass = 1.0f;
      UpdateVertices(SimState, 0);

      
      SimState->Cubes[1].Position = glm::vec3(2.0,
					      2.0,
					      2.0);
      SimState->Cubes[1].Velocity = glm::vec3(0.0);
      SimState->Cubes[1].Forces = glm::vec3(0.0);
      SimState->Cubes[1].XAngle = 0.0f;
      SimState->Cubes[1].YAngle = 0.0f;
      SimState->Cubes[1].ZAngle = 0.0f;
      SimState->Cubes[1].Size = 1.0f;
      SimState->Cubes[1].Mass = 1.0f;
      UpdateVertices(SimState, 1);

      // NOTE(Jovan): Sphere init
      SimState->SphereCount = 1;
      SimState->Spheres[0].Position = glm::vec3(2.0,
						2.0,
						0.0);
      SimState->Spheres[0].XAngle = 0.0f;
      SimState->Spheres[0].YAngle = 0.0f;
      SimState->Spheres[0].ZAngle = 0.0f;
      SimState->Spheres[0].Radius = 0.5f;
      
      Memory->IsInitialized = 1;
    }
  // NOTE(Jovan): End init
  // ---------------------

  // NOTE(Jovan): Coordinate systems
  // -------------------------------
  glm::mat4 Projection = glm::perspective(glm::radians(45.0f),
					  (real32)DEFAULT_WINDOW_WIDTH / (real32)DEFAULT_WINDOW_HEIGHT,
					  0.1f,
					  100.0f);
  glm::mat4 Model = glm::mat4(1.0f);
  glm::mat4 View = glm::mat4(1.0f);
  View = glm::translate(View, glm::vec3(0.0f, 0.0f, -3.0f));
  View = glm::lookAt(SimState->Camera.Position, SimState->Camera.Position + SimState->Camera.Front,
		     SimState->Camera.Up);

  // NOTE(Jovan): End coordinate systems
  // ------------------------------------
  
  UpdateCamera(SimState, Input);

  // NOTE(Jovan): Input
  // ------------------
  real32 CameraSpeed = 0.05f;
  if(Input->KeyboardController.MoveForward.EndedDown)
    {
      SimState->Camera.Position += CameraSpeed * SimState->Camera.Front;
    }
  if(Input->KeyboardController.MoveLeft.EndedDown)
    {
      SimState->Camera.Position -= glm::normalize(glm::cross(SimState->Camera.Front, SimState->Camera.Up)) * CameraSpeed;
    }
  if(Input->KeyboardController.MoveBack.EndedDown)
    {
      SimState->Camera.Position -= CameraSpeed * SimState->Camera.Front;
    }
  if(Input->KeyboardController.MoveRight.EndedDown)
    {
      SimState->Camera.Position += glm::normalize(glm::cross(SimState->Camera.Front, SimState->Camera.Up)) * CameraSpeed;
    }
  if(Input->KeyboardController.ShootAction.EndedDown)
    {
      CubeAddForce(SimState, 0, 20.0f * SimState->Camera.Front);
    }
  if(Input->KeyboardController.DebugLeft.EndedDown)
    {
      SimState->Cubes[0].Position += dt * glm::vec3(1.0, 0.0, 0.0);
    }
  if(Input->KeyboardController.DebugRight.EndedDown)
    {
      SimState->Cubes[0].Position += dt * glm::vec3(-1.0, 0.0, 0.0);
    }
  if(Input->KeyboardController.DebugUp.EndedDown)
    {
      SimState->Cubes[0].Position += dt * glm::vec3(0.0, 1.0, 0.0);
    }
  if(Input->KeyboardController.DebugDown.EndedDown)
    {
      SimState->Cubes[0].Position += dt * glm::vec3(0.0, -1.0, 0.0);
    }
  if(Input->KeyboardController.DebugForward.EndedDown)
    {
      SimState->Cubes[0].Position += dt * glm::vec3(0.0, 0.0, 1.0);
    }
  if(Input->KeyboardController.DebugBack.EndedDown)
    {
      SimState->Cubes[0].Position += dt * glm::vec3(0.0, 0.0, -1.0);
    }

  // TODO(Add normal constraint)
  if(SimState->Camera.Position.y < 0.5)
    {
      SimState->Camera.Position.y = 0.5;
    }
  
  // NOTE(Jovan) End input
  // ---------------------

  // NOTE(Jovan): Physics stuff
  // --------------------------
  // TODO(Jovan): Implement
  bool32 CollisionHappened = DetectCollisions(SimState, 0, 1, CC);
  // printf("A: %f, %f, %f | B: %f, %f, %f | C:%d\n",
  // 	 SimState->Cubes[0].Position.x,
  // 	 SimState->Cubes[0].Position.y,
  // 	 SimState->Cubes[0].Position.z,
  // 	 SimState->Cubes[1].Position.x,
  // 	 SimState->Cubes[1].Position.y,
  // 	 SimState->Cubes[1].Position.z,
  // 	 CollisionHappened);
  // TODO(Jovan): POSSIBLE MEMORY DEATH, but collision detection works
  // TODO(Jovan): Write own, memory safe, push function for array
  ResetVertices(SimState);
  ResolveCollision();

  for(uint32 CubeIndex = 0;
      CubeIndex < SimState->CubeCount;
      ++CubeIndex)
    {
      phys_return Y0 = {};
      phys_return Y = {};
      Y0.X = SimState->Cubes[CubeIndex].Position;
      Y0.Y = SimState->Cubes[CubeIndex].Velocity;
      Y = Euler(MovementFunction, dt, Y0, SimState->Cubes[CubeIndex].Forces, SimState->Cubes[CubeIndex].Mass);
      SimState->Cubes[CubeIndex].Position = Y.X;
      SimState->Cubes[CubeIndex].Velocity = Y.Y;
      UpdateVertices(SimState, CubeIndex);
      CubeClearForces(SimState, CubeIndex);
    }
  
  // NOTE(Jovan): End physics stuff
  // ------------------------------
  
  // TODO(Jovan): Maybe move to sdl_camera???
  glm::vec3 Front;
  Front.x = cos(glm::radians(SimState->Camera.Yaw)) * cos(glm::radians(SimState->Camera.Pitch));
  Front.y = sin(glm::radians(SimState->Camera.Pitch));
  Front.z = sin(glm::radians(SimState->Camera.Yaw)) * cos(glm::radians(SimState->Camera.Pitch));
  SimState->Camera.Front = glm::normalize(Front);
  View = glm::lookAt(SimState->Camera.Position, SimState->Camera.Position + SimState->Camera.Front,
		     SimState->Camera.Up);

  // NOTE(Jovan): Coordinate drawing
  glUseProgram(Render->Shaders[2]);
  SetUniformM4(Render->Shaders[2], "View", View);
  SetUniformM4(Render->Shaders[2], "Projection", Projection);
  real32 LineLength = 100.0f;
  for(uint32 CubeIndex = 0;
      CubeIndex < SimState->CubeCount;
      ++CubeIndex)
    {
      for(uint32 i = 0;
	  i < 8;
	  ++i)
	{
	  Model = glm::mat4(1.0);
	  Model = glm::translate(Model, SimState->Cubes[CubeIndex].Vertices[i]);
	  Model = glm::scale(Model, glm::vec3(LineLength));
	  SetUniformM4(Render->Shaders[2], "Model", Model);
	  glBindVertexArray(Render->VAOs[3]);
	  glDrawArrays(GL_LINES, 0, 6);
	}
    }
  glBindVertexArray(0);
  
  // NOTE(Jovan): Floor drawing
  // --------------------------

  glUseProgram(Render->Shaders[0]);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, Render->Textures[2]);
  SetUniformM4(Render->Shaders[0], "View", View);
  SetUniformM4(Render->Shaders[0], "Projection", Projection);
  real32 FloorSize = 1;
  Model = glm::mat4(1.0);
  Model = glm::translate(Model, glm::vec3(0.0, 0.0, 0.0));
  Model = glm::scale(Model, glm::vec3(FloorSize, 0, FloorSize));
  SetUniformM4(Render->Shaders[0], "Model", Model);
  glBindVertexArray(Render->VAOs[2]);
  glDrawArrays(GL_TRIANGLES, 0, 6);
  glBindVertexArray(0);
			 
  // NOTE(Jovan): End floor drawing
  // ------------------------------
  SetUniformM4(Render->Shaders[0], "View", View);
  SetUniformM4(Render->Shaders[0], "Projection", Projection);
  
  // NOTE(Jovan): Cube drawing
  // -------------------------
  glUseProgram(Render->Shaders[0]);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, Render->Textures[0]);

  SetUniformM4(Render->Shaders[0], "View", View);
  SetUniformM4(Render->Shaders[0], "Projection", Projection);
#if 1
  for(uint32 CubeIndex = 0;
      CubeIndex < SimState->CubeCount;
      ++CubeIndex)
    {
      Model = glm::mat4(1.0f);
      Model = glm::translate(Model, SimState->Cubes[CubeIndex].Position);
      Model = glm::rotate(Model, glm::radians(SimState->Cubes[CubeIndex].XAngle), glm::vec3(1.0f, 0.0f, 0.0f));
      Model = glm::rotate(Model, glm::radians(SimState->Cubes[CubeIndex].YAngle), glm::vec3(0.0f, 1.0f, 0.0f));
      Model = glm::rotate(Model, glm::radians(SimState->Cubes[CubeIndex].ZAngle), glm::vec3(0.0f, 0.0f, 1.0f));
      Model = glm::scale(Model, glm::vec3(SimState->Cubes[CubeIndex].Size,
					  SimState->Cubes[CubeIndex].Size,
					  SimState->Cubes[CubeIndex].Size));
      SetUniformM4(Render->Shaders[0], "Model", Model);
      glBindVertexArray(Render->VAOs[0]);
      glDrawArrays(GL_TRIANGLES, 0, 36);
    }
#endif
  glBindTexture(GL_TEXTURE_2D, 0);
  glBindVertexArray(0);
  // NOTE(Jovan): End cube drawing
  // -----------------------------

  // NOTE(Jovan): Sphere drawing
  // ---------------------------
#if 1
  glUseProgram(Render->Shaders[1]);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, Render->Textures[1]);
  SetUniformM4(Render->Shaders[1], "View", View);
  SetUniformM4(Render->Shaders[1], "Projection", Projection);
  for(uint32 SphereIndex = 0;
      SphereIndex < SimState->SphereCount;
      ++SphereIndex)
    {
      Model = glm::mat4(1.0);
      Model = glm::translate(Model, SimState->Spheres[SphereIndex].Position);
      Model = glm::rotate(Model, glm::radians(SimState->Spheres[SphereIndex].XAngle),
			  glm::vec3(1.0f, 0.0f, 0.0f));
      Model = glm::rotate(Model, glm::radians(SimState->Spheres[SphereIndex].YAngle),
			  glm::vec3(0.0f, 1.0f, 0.0f));
      Model = glm::rotate(Model, glm::radians(SimState->Spheres[SphereIndex].ZAngle),
			  glm::vec3(0.0f, 0.0f, 1.0f));
      Model = glm::scale(Model, glm::vec3(SimState->Spheres[SphereIndex].Radius,
					  SimState->Spheres[SphereIndex].Radius,
					  SimState->Spheres[SphereIndex].Radius));
      SetUniformM4(Render->Shaders[1], "Model", Model);
      glBindVertexArray(Render->VAOs[1]);
      glDrawElements(GL_TRIANGLES, Render->Num,  GL_UNSIGNED_INT, Render->Indices);
    }
  glBindVertexArray(0);
#endif
      
  // NOTE(Jovan): End sphere drawing
  // -------------------------------

  // NOTE(Jovan): Logging
  // --------------------

  // TODO(Jovan): Make better logging

  printf("Collision: %d\n", CollisionHappened);
  
  // NOTE(Jovan): End logging
  // ------------------------
        
}
