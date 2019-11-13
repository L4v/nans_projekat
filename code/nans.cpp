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
DetectCollisions(sdl_state* State, real32 dt)
{
  // TODO(Jovan): Implement
}

internal void
ResolveCollision()
{
  // TODO(Jovan): Implement
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

extern "C" SIM_UPDATE_AND_RENDER(SimUpdateAndRender)
{
  sdl_state* SimState = (sdl_state*)Memory->PermanentStorage;

  // NOTE(Jovan): Init
  // -----------------
  if(!Memory->IsInitialized)
    {

      // NOTE(Jovan): Init random seed
      srand(time(0));
      
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
      SimState->CubeCount = 1;
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
  DetectCollisions(SimState, dt);
  ResolveCollision();

  phys_return Y0 = {};
  phys_return Y = {};
  Y0.X = SimState->Cubes[0].Position;
  Y0.Y = SimState->Cubes[0].Velocity;
  Y = Euler(MovementFunction, dt, Y0, SimState->Cubes[0].Forces, SimState->Cubes[0].Mass);
  SimState->Cubes[0].Position = Y.X;
  SimState->Cubes[0].Velocity = Y.Y;
  UpdateVertices(SimState, 0);
  CubeClearForces(SimState, 0);

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
  for(uint32 i = 0;
      i < 8;
      ++i)
    {
      Model = glm::mat4(1.0);
      Model = glm::translate(Model, SimState->Cubes[0].Vertices[i]);
      Model = glm::scale(Model, glm::vec3(LineLength));
      SetUniformM4(Render->Shaders[2], "Model", Model);
      glBindVertexArray(Render->VAOs[3]);
      glDrawArrays(GL_LINES, 0, 6);
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
        
}
