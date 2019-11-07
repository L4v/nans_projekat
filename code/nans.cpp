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

extern "C" SIM_UPDATE_AND_RENDER(SimUpdateAndRender)
{
  sdl_state* SimState = (sdl_state*)Memory->PermanentStorage;
  if(!Memory->IsInitialized)
    {
	  
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
	  
      // TODO(Jovan): Remove, only for testing
      for(int CubeIndex = 0;
	  CubeIndex < 10;
	  ++CubeIndex)
	{
	  SimState->Positions[CubeIndex] = glm::vec3(0.5*CubeIndex,
						     1.0*CubeIndex,
						     1.5*CubeIndex);
	}
      
      Memory->IsInitialized = 1;
    }

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
  
  // NOTE(Jovan): Camera update
  UpdateCamera(SimState, Input);
  
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

  // TODO(Jovan): Maybe move to sdl_camera???
  glm::vec3 Front;
  Front.x = cos(glm::radians(SimState->Camera.Yaw)) * cos(glm::radians(SimState->Camera.Pitch));
  Front.y = sin(glm::radians(SimState->Camera.Pitch));
  Front.z = sin(glm::radians(SimState->Camera.Yaw)) * cos(glm::radians(SimState->Camera.Pitch));
  SimState->Camera.Front = glm::normalize(Front);
  View = glm::lookAt(SimState->Camera.Position, SimState->Camera.Position + SimState->Camera.Front,
		     SimState->Camera.Up);
  
  // NOTE(Jovan): Cube drawing
  // -------------------------
  glUseProgram(Render->Shaders[0]);
  glBindTexture(GL_TEXTURE_2D, Render->Textures[0]);

  SetUniformM4(Render->Shaders[0], "View", View);
  SetUniformM4(Render->Shaders[0], "Projection", Projection);

  for(int CubeIndex = 0;
      CubeIndex < 10;
      ++CubeIndex)
    {
      int32 CubeSize = 1.0f;
      local_persist real32 Angle = 50.0f;
      Angle += dt * 2.0f;
      Model = glm::mat4(1.0f);
      Model = glm::scale(Model, glm::vec3(CubeSize, CubeSize, CubeSize));
      Model = glm::translate(Model, SimState->Positions[CubeIndex]);
      Model = glm::rotate(Model, glm::radians(Angle), glm::vec3(1.0f, 1.0f, 1.0f));
      SetUniformM4(Render->Shaders[0], "Model", Model);
      glBindVertexArray(Render->VAOs[0]);
      glDrawArrays(GL_TRIANGLES, 0, 36);
    }

  // NOTE(Jovan): End cube drawing
  // -----------------------------
        
}
