#include "nans.h"
#define Pi32 3.14159265359f
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
  // NOTE(Jovan): Sphere testing
  // ---------------------------

  const int32 na = 36;
  const int32 nb = 18;
  const int32 na3 = na * 3;
  const int32 nn = nb*na3;
  real32 SpherePos[nn];
  real32 SphereNor[nn];
  uint32 SphereIX[na*(nb-1)*6];
  uint32 SphereVBO[4];
  uint32 SphereVAO[4];
  
  // NOTE(Jovan): End sphere testing
  // -------------------------------
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

      // NOTE(Jovan): SPHERE
      // -------------------
#if 0
      real32 x, y, z, a, b, da, db, r = 3.5f;
      int32 ia, ib, ix, iy;
      da = 2.0 * Pi32 / real32(na);
      db = Pi32 / real32(nb - 1);

      for(ix = 0, b = -0.5f * Pi32, ib = 0;
	  ib < nb;
	  ++ib, b+=db)
	{
	  for(a = 0.0, ia = 0;
	      ia < na;
	      ++ia, a += da, ix += 3)
	    {
	      x = cos(b) * cos(a);
	      y = cos(b) * sin(a);
	      z = sin(b);

	      SpherePos[ix + 0] = x * r;
	      SpherePos[ix + 1] = y * r;
	      SpherePos[ix + 2] = z * r;
	      SphereNor[ix + 0] = x;
	      SphereNor[ix + 1] = y;
	      SphereNor[ix + 2] = z;
	    }
	}
      for(ix = 0, iy = 0, ib = 1;
	  ib < nb;
	  ++ib)
	{
	  for(ia = 1;
	      ia < na;
	      ++ia, ++iy)
	    {
	      SphereIX[ix] = iy; ix++;
	      SphereIX[ix] = iy + 1; ix++;
	      SphereIX[ix] = iy + na; ix++;
		  
	      SphereIX[ix] = iy + na; ix++;
	      SphereIX[ix] = iy + 1; ix++;
	      SphereIX[ix] = iy + na + 1; ix++;
		  
	    }
	  SphereIX[ix] = iy; ix++;
	  SphereIX[ix] = iy + 1 - na; ix++;
	  SphereIX[ix] = iy + na; ix++;
		  
	  SphereIX[ix] = iy + na; ix++;
	  SphereIX[ix] = iy - na + 1; ix++;
	  SphereIX[ix] = iy + 1; ix++;
	  iy++;
	}

      uint32 i;
      glGenVertexArrays(4, SphereVAO);
      glGenBuffers(4, SphereVBO);
      glBindVertexArray(SphereVAO[0]);
      i = 0;
      glBindBuffer(GL_ARRAY_BUFFER,SphereVBO[i]);
      glBufferData(GL_ARRAY_BUFFER,sizeof(SpherePos),SpherePos,GL_STATIC_DRAW);
      glEnableVertexAttribArray(i);
      glVertexAttribPointer(i,3,GL_FLOAT,GL_FALSE,0,0);
      i=1; // indices
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,SphereVBO[i]);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(SphereIX),SphereIX,GL_STATIC_DRAW);
      glEnableVertexAttribArray(i);
      glVertexAttribPointer(i,4,GL_UNSIGNED_INT,GL_FALSE,0,0);
      i=2; // normal
      glBindBuffer(GL_ARRAY_BUFFER,SphereVBO[i]);
      glBufferData(GL_ARRAY_BUFFER,sizeof(SphereNor),SphereNor,GL_STATIC_DRAW);
      glEnableVertexAttribArray(i);
      glVertexAttribPointer(i,3,GL_FLOAT,GL_FALSE,0,0);

      
      glBindVertexArray(0);
      glBindBuffer(GL_ARRAY_BUFFER,0);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
      glDisableVertexAttribArray(0);
      glDisableVertexAttribArray(1);
      glDisableVertexAttribArray(2);
      glDisableVertexAttribArray(3);
#endif
      // NOTE(Jovan): END SPHERE
      // -----------------------
      
      
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

  // NOTE(Jovan): SPHERE
  // -------------------

#if 0
  float aspect=float(DEFAULT_WINDOW_WIDTH)/float(DEFAULT_WINDOW_HEIGHT);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0/aspect,aspect,0.1,100.0);
  glMatrixMode(GL_TEXTURE);
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(0.0,0.0,-10.0);

  glEnable(GL_DEPTH_TEST);
  glDisable(GL_TEXTURE_2D);
  
  glEnable(GL_CULL_FACE);
  glFrontFace(GL_CCW);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  glBindVertexArray(SphereVAO[0]);
  glDrawElements(GL_TRIANGLES,sizeof(SphereIX)/sizeof(GLuint),GL_UNSIGNED_INT,0);
  glBindVertexArray(0);
#endif
  // NOTE(Jovan): End SPHERE
  // -----------------------
  
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
      // TODO(Jovan): Remove local persist, use game state
      local_persist real32 Angle = 50.0f;
      Angle += dt * +5.0f;
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
