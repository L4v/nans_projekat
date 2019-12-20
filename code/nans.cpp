#include "nans.h"

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
  if(State->Camera.Position.y < 0.5)
    {
      State->Camera.Position.y = 0.5;
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

internal phys_return
RotationFunction(glm::vec3 RotVelocity, glm::vec3 SummedTorque, real32 MomentOfInertia)
{
  phys_return Result = {};
  Result.X = RotVelocity;
  Result.Y = (real32)(1.0f / MomentOfInertia) * (SummedTorque - GLOBAL_FRICTION * RotVelocity);
  return Result;
}

internal void
CubeAddForce(sdl_state* State, int32 CubeIndex, glm::vec3 Force)
{
  State->Cubes[CubeIndex].Forces += Force;
};

internal void
CubeAddTorque(sdl_state* State, int32 CubeIndex, glm::vec3 Torque)
{
  State->Cubes[CubeIndex].Torque += Torque;
};

internal void
CubeClearForces(sdl_state* State, int32 CubeIndex)
{
  State->Cubes[CubeIndex].Forces = glm::vec3(0.0);
  State->Cubes[CubeIndex].Torque = glm::vec3(0.0);
};


internal void
HandleInput(sdl_state* State, sdl_input* Input, real32 dt)
{

  if(Input->KeyboardController.MoveForward.EndedDown)
    {
      State->Camera.Position += State->Camera.Speed * State->Camera.Front;
    }
  if(Input->KeyboardController.MoveLeft.EndedDown)
    {
      State->Camera.Position -= glm::normalize(glm::cross(State->Camera.Front, State->Camera.Up)) * State->Camera.Speed;
    }
  if(Input->KeyboardController.MoveBack.EndedDown)
    {
      State->Camera.Position -= State->Camera.Speed * State->Camera.Front;
    }
  if(Input->KeyboardController.MoveRight.EndedDown)
    {
      State->Camera.Position += glm::normalize(glm::cross(State->Camera.Front, State->Camera.Up)) * State->Camera.Speed;
    }
  if(Input->KeyboardController.ShootAction.EndedDown)
    {
      CubeAddForce(State, 0, 20.0f * State->Camera.Front);
    }
  if(Input->KeyboardController.DebugLeft.EndedDown)
    {
      State->Cubes[0].Position += dt * glm::vec3(1.0, 0.0, 0.0);
    }
  if(Input->KeyboardController.DebugRight.EndedDown)
    {
      State->Cubes[0].Position += dt * glm::vec3(-1.0, 0.0, 0.0);
    }
  if(Input->KeyboardController.DebugUp.EndedDown)
    {
      State->Cubes[0].Position += dt * glm::vec3(0.0, 1.0, 0.0);
    }
  if(Input->KeyboardController.DebugDown.EndedDown)
    {
      State->Cubes[0].Position += dt * glm::vec3(0.0, -1.0, 0.0);
    }
  if(Input->KeyboardController.DebugForward.EndedDown)
    {
      State->Cubes[0].Position += dt * glm::vec3(0.0, 0.0, 1.0);
    }
  if(Input->KeyboardController.DebugBack.EndedDown)
    {
      State->Cubes[0].Position += dt * glm::vec3(0.0, 0.0, -1.0);
    }
  // TODO(Jovan): For debugging
  if(Input->KeyboardController.DebugReset.EndedDown)
    {
      State->Cubes[0].Position = glm::vec3(2.0, 2.6, 2.0);
      State->Cubes[0].V = glm::vec3(0.0f);
      State->Cubes[0].W = glm::vec3(0.0f);
      State->Cubes[0].Angles = glm::vec3(0.0f);
      
      State->Cubes[1].Position = glm::vec3(2.0, 1.5, 2.0);
      State->Cubes[1].V = glm::vec3(0.0f);
      State->Cubes[1].W = glm::vec3(0.0f);
      State->Cubes[1].Angles = glm::vec3(0.0f);
      
      State->Cubes[2].Position = glm::vec3(2.0, 4.0, 2.0);
      State->Cubes[2].V = glm::vec3(0.0f);
      State->Cubes[2].W = glm::vec3(0.0f);
      State->Cubes[2].Angles = glm::vec3(0.0f);
    }
}

internal void
RemoveVertex(simplex* Simplex, uint32 VertexIndex)
{
  Assert(VertexIndex < Simplex->Count);
  for(uint32 i = VertexIndex;
      i < Simplex->Count - 1;
      ++i)
    {
      Simplex->Vertices[i] = Simplex->Vertices[i + 1];
    }
  Simplex->Count--;
}
internal void
PushVertex(simplex* Simplex, vertex V)
{
  Simplex->Vertices[Simplex->Count] = V;
  Simplex->Count++;
}

internal inline void
ClearVertices(simplex* Simplex)
{
  Simplex->Count = 0;
}

internal void
RemoveEdge(edge* Edge, uint32 EdgeIndex)
{
  Assert(EdgeIndex < Edge->Count);
  for(uint32 i = EdgeIndex;
      i < Edge->Count - 1;
      ++i)
    {
      Edge->A[i] = Edge->A[i + 1];
      Edge->B[i] = Edge->B[i + 1];
    }
  Edge->Count--;
}

internal void
PushEdge(edge* Edge, vertex A, vertex B)
{
  // NOTE(Jovan): Checks whether an edge of opposite
  // winding already exists on the list.
  // If it does, remove it and don't add the new one
  for(uint32 EdgeIndex = 0;
      EdgeIndex < Edge->Count;
      ++EdgeIndex)
    {
      if(Edge->A[EdgeIndex] == B && Edge->B[EdgeIndex] == A)
	{
	  RemoveEdge(Edge, EdgeIndex);
	  return;
	}
    }
  // NOTE(Jovan): The edge didn't exist, so it's being added
  Edge->A[Edge->Count] = A;
  Edge->B[Edge->Count] = B;
  Edge->Count++;
}

internal void
ClearEdges(edge* Edge)
{
  Edge->Count = 0;
}

internal void
RemoveTriangle(triangle* Triangle, uint32 TriangleIndex)
{
  Assert(TriangleIndex < Triangle->Count);
  for(uint32 i = TriangleIndex;
      i < Triangle->Count - 1;
      ++i)
    {
      Triangle->A[i] = Triangle->A[i+1];
      Triangle->B[i] = Triangle->B[i+1];
      Triangle->C[i] = Triangle->C[i+1];
      Triangle->N[i] = Triangle->N[i+1];
    }
  Triangle->Count--;
}

// NOTE(Jovan): Push triangle using vertices
internal void
PushTriangle(triangle* Triangle, vertex A, vertex B, vertex C)
{
  Triangle->A[Triangle->Count] = A;
  Triangle->B[Triangle->Count] = B;
  Triangle->C[Triangle->Count] = C;
  
  // NOTE(Jovan): Calculate normal and make sure it's
  // pointing away from the origin
  Triangle->N[Triangle->Count] = glm::cross(B.P - A.P, C.P - A.P);
  // NOTE(Jovan): Normalize vector
  //real32 Len = sqrt(glm::dot(Triangle->N[Triangle->Count], Triangle->N[Triangle->Count]));
  Triangle->N[Triangle->Count] = glm::normalize(Triangle->N[Triangle->Count]);
  
  if(glm::dot(Triangle->A[Triangle->Count].P, Triangle->N[Triangle->Count]) < 0)
    {
      Triangle->N[Triangle->Count] *= -1.0f;
    }
  
  Triangle->Count++;
}

internal inline void
ClearTriangles(triangle* Triangle)
{
  Triangle->Count = 0;
}

internal void
FloorUpdateVertices(sdl_state* State)
{
  glm::mat4 Model = State->Floor.Model;
  
  State->Floor.Vertices[0] = glm::vec3(Model * glm::vec4(0.5, 0.5, 0.5, 1.0));
  State->Floor.Vertices[1] = glm::vec3(Model * glm::vec4(0.5, 0.5, -0.5, 1.0));
  State->Floor.Vertices[2] = glm::vec3(Model * glm::vec4(-0.5, 0.5, 0.5, 1.0));
  State->Floor.Vertices[3] = glm::vec3(Model * glm::vec4(-0.5, 0.5, -0.5, 1.0));
  State->Floor.Vertices[4] = glm::vec3(Model * glm::vec4(0.5, -0.5, 0.5, 1.0));
  State->Floor.Vertices[5] = glm::vec3(Model * glm::vec4(0.5, -0.5, -0.5, 1.0));
  State->Floor.Vertices[6] = glm::vec3(Model * glm::vec4(-0.5, -0.5, 0.5, 1.0));
  State->Floor.Vertices[7] = glm::vec3(Model * glm::vec4(-0.5, -0.5, -0.5, 1.0));
}

internal void
UpdateVertices(sdl_state* State, int32 CubeIndex)
{
  glm::mat4 Model = State->Cubes[CubeIndex].Model;
  
  State->Cubes[CubeIndex].Vertices[0] = glm::vec3(Model * glm::vec4(0.5, 0.5, 0.5, 1.0));
  State->Cubes[CubeIndex].Vertices[1] = glm::vec3(Model * glm::vec4(0.5, 0.5, -0.5, 1.0));
  State->Cubes[CubeIndex].Vertices[2] = glm::vec3(Model * glm::vec4(-0.5, 0.5, 0.5, 1.0));
  State->Cubes[CubeIndex].Vertices[3] = glm::vec3(Model * glm::vec4(-0.5, 0.5, -0.5, 1.0));
  State->Cubes[CubeIndex].Vertices[4] = glm::vec3(Model * glm::vec4(0.5, -0.5, 0.5, 1.0));
  State->Cubes[CubeIndex].Vertices[5] = glm::vec3(Model * glm::vec4(0.5, -0.5, -0.5, 1.0));
  State->Cubes[CubeIndex].Vertices[6] = glm::vec3(Model * glm::vec4(-0.5, -0.5, 0.5, 1.0));
  State->Cubes[CubeIndex].Vertices[7] = glm::vec3(Model * glm::vec4(-0.5, -0.5, -0.5, 1.0));
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
GetSphereSupport(sdl_state* State, int32 SphereIndex, glm::vec3 Direction)
{
  real32 MaxDistance = -FLT_MAX;
  glm::vec3 Result = glm::vec3(0.0f);
  sphere* Sphere = &State->Spheres[SphereIndex];
  Result = Sphere->Position + Sphere->Radius * Direction;
  
  return Result;
}

internal glm::vec3
GetFloorSupport(sdl_state* State, glm::vec3 Direction)
{
  real32 MaxDistance = -FLT_MAX;
  glm::vec3 Result = glm::vec3(0.0f);

  // NOTE(Jovan) 8 vertices in a cube
  for(uint32 VertexIndex = 0;
      VertexIndex < 8;
      ++VertexIndex)
    {
      glm::vec3 CurrVertex = State->Floor.Vertices[VertexIndex];
      real32 Distance = glm::dot(CurrVertex, Direction);
      if(Distance > MaxDistance)
	{
	  MaxDistance = Distance;
	  Result = CurrVertex;
	}
    }

  return Result;
}

internal vertex
CalculateSupport(sdl_state* State, glm::vec3 Direction)
{
  vertex Result = {};
  glm::vec3 SupportA, SupportB;
  switch(State->CurrentCollisionType)
    {
    case CC:
      {
	SupportA = GetCubeSupport(State, State->IndexA, Direction);
	SupportB = GetCubeSupport(State, State->IndexB, -1.0f * Direction);
      }break;
    case CS:
      {
	SupportA = GetCubeSupport(State, State->IndexA, Direction);
	SupportB = GetSphereSupport(State, State->IndexB, Direction);
      }break;
    case SC:
      {
	SupportA = GetSphereSupport(State, State->IndexA, Direction);	
	SupportB = GetCubeSupport(State, State->IndexB, Direction);
      }break;
    case CF:
      {
	SupportA = GetCubeSupport(State, State->IndexA, Direction);
	SupportB = GetFloorSupport(State, Direction);
	State->Spheres[0].Position = SupportA;
	State->Spheres[1].Position = SupportB;
      }break;
    case FC:
      {
	SupportA = GetFloorSupport(State, Direction);
	SupportB = GetCubeSupport(State, State->IndexB, Direction);
      }break;
    case SS:
      {
	SupportA = GetSphereSupport(State, State->IndexA, Direction);
	SupportB = GetSphereSupport(State, State->IndexB, Direction);
      }break;
    case SF:
      {
	SupportA = GetSphereSupport(State, State->IndexA, Direction);
	SupportB = GetFloorSupport(State, Direction);
      }break;
    case FS:
      {
	SupportA = GetFloorSupport(State, Direction);
	SupportB = GetSphereSupport(State, State->IndexB, Direction);
      }break;
    }
  Result.P = SupportA - SupportB;
  Result.SupA = SupportA;
  PushVertex(State->Simplex, Result);

  return Result;
}

internal bool32
AddSupport(sdl_state* State, glm::vec3 Direction)
{
  bool32 Result = 0;
  vertex NewVertex = {};
  NewVertex = CalculateSupport(State, Direction);
  
  if(glm::dot(Direction, NewVertex.P) >= 0)
    {
      Result = 1;
    }
  else
    {
      Result = 0;
    }
  return Result;
}

internal glm::vec3
ClosestPointOnLine(glm::vec3 A, glm::vec3 B, glm::vec3 Target,
		   real32* U, real32* V)
{
  glm::vec3 Result = glm::vec3(0.0f);
  glm::vec3 LineSegment = glm::normalize(B - A);
  real32 Length = glm::length(LineSegment);
  *V = glm::dot(-A, LineSegment) / Length;
  *U = glm::dot(B, LineSegment) / Length;
  if(*U <= 0.0f)
    {
      Result = B;
    }
  else if(*V <= 0.0f)
    {
      Result = A;
    }
  else
    {
      Result = (*U) * A + (*V) * B;
    }

  return Result;
}

internal evolve_result
EvolveSimplex(sdl_state* State, glm::vec3 PositionA, glm::vec3 PositionB)
{
  // IMPORTANT TODO(Jovan): REDO GJK
  evolve_result Result = NoIntersection;
  glm::vec3 Direction = glm::normalize(glm::vec3(1.0f));
  simplex* Simplex = State->Simplex;
  uint32 NoVertices = Simplex->Count;
  switch(NoVertices)
    {
    case 0:
      {
	// NOTE(Jovan): Best presumption of dir is the relative direction
	// of the 2 shapes
	
	//Direction = PositionB - PositionA;
      }break;
    case 1:
      {
	// NOTE(Jovan): Flip the direction to point to the origin
	Direction *= -1;
      }break;
    case 2:
      {
	// // NOTE(Jovan): Form line from first 2 vertices
	glm::vec3 AB = Simplex->Vertices[1].P - Simplex->Vertices[0].P;
	// // NOTE(Jovan): Form line from origin to A
	glm::vec3 A0 = -1.0f * Simplex->Vertices[0].P;

	// // NOTE(Jovan): Direction perpendicular to AB in the direction
	// // of the origin
	//Direction = glm::cross(glm::cross(AB, A0), AB);
	// /* STUDY TODO IMPORTANT(Jovan):
	//  * Possible location for collision error
	//  * Direction vector will be a 0 vector if the origin is
	//  * on the edge of the Miknowski sum
	// */
	// if(Direction == glm::vec3(0.0))
	//   {
	//     // TODO(Jovan): How to handle???
	//   }
	real32 U = 0.0f;
	real32 V = 0.0f;
	glm::vec3 Origin = glm::vec3(0.0);
	glm::vec3 ClosestPoint = ClosestPointOnLine(Simplex->Vertices[0].P,
						    Simplex->Vertices[1].P,
						    Origin, &U, &V);
	if(V <= 0.0f)
	  {
	    RemoveVertex(Simplex, 1);
	    Direction = -ClosestPoint;
	  }
	else if(U <= 0.0f)
	  {
	    RemoveVertex(Simplex, 0);
	    Direction = -ClosestPoint;
	  }
      	else
      	  {
      	    Direction = -ClosestPoint;
      	  }
      }break;
    case 3:
      {
	glm::vec3 AC = Simplex->Vertices[2].P - Simplex->Vertices[0].P;
	glm::vec3 AB = Simplex->Vertices[1].P - Simplex->Vertices[0].P;
	Direction = glm::cross(AC, AB);

	// NOTE(Jovan): Ensure that Direction points to the origin
	glm::vec3 A0 = -1.0f * Simplex->Vertices[0].P;
	if(glm::dot(Direction, A0) < 0)
	  {
	    Direction *= -1.0f;
	  }
      }break;
    case 4:
      {
	// NOTE(Jovan): 3 edges of interest
	glm::vec3 DA = Simplex->Vertices[0].P - Simplex->Vertices[3].P;
	glm::vec3 DB = Simplex->Vertices[1].P - Simplex->Vertices[3].P;
	glm::vec3 DC = Simplex->Vertices[2].P - Simplex->Vertices[3].P;
	
	// NOTE(Jovan): Dir to the origin
	glm::vec3 D0 = -1.0f * Simplex->Vertices[3].P;

	// NOTE(Jovan): Check triangles ABD, BCD, CAD
	glm::vec3 ABDNorm = glm::cross(DA, DB);
	glm::vec3 BCDNorm = glm::cross(DB, DC);
	glm::vec3 CADNorm = glm::cross(DC, DA);
	
	if(glm::dot(ABDNorm, D0) > 0.0f)
	  {
	    // NOTE(Jovan): Origin outside of ABD -> remove C
	    RemoveVertex(Simplex, 2);
	    Direction = ABDNorm;
	  }
	else if(glm::dot(BCDNorm, D0) > 0.0f)
	  {
	    // NOTE(Jovan): Origin outside of BCD -> remove A
	    RemoveVertex(Simplex, 0);
	    Direction = BCDNorm;
	  }
	else if(glm::dot(CADNorm, D0) > 0.0f)
	  {
	    // NOTE(Jovan): Origin outside of CAD -> remove B
	    RemoveVertex(Simplex, 1);
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
  if(glm::length(Direction) <= 0.0001f)
    {
      Result = NoIntersection;
      return Result;
    }
  if(AddSupport(State, Direction))
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
DetectCollisions(sdl_state* State)
{
  bool32 Result = 0;
  evolve_result EvolutionResult = StillEvolving;
  int32 AIndex = State->IndexA;
  int32 BIndex = State->IndexB;
  while((EvolutionResult == StillEvolving) && (State->GJKIteration++ <= MAX_GJK_ITERATIONS))
    {
      switch(State->CurrentCollisionType)
	{
	case CC:
	  {
	    EvolutionResult = EvolveSimplex(State, State->Cubes[AIndex].Position,
					    State->Cubes[BIndex].Position);
	  }break;
	case CS:
	  {
	    EvolutionResult = EvolveSimplex(State, State->Cubes[AIndex].Position,
					    State->Spheres[BIndex].Position);
	  }break;
	case SC:
	  {
	    EvolutionResult = EvolveSimplex(State, State->Spheres[AIndex].Position,
					    State->Cubes[BIndex].Position);
	  }break;
	case CF:
	  {
	    EvolutionResult = EvolveSimplex(State, State->Cubes[AIndex].Position,
					    State->Floor.Position);
	  }break;
	case FC:
	  {
	    EvolutionResult = EvolveSimplex(State, State->Floor.Position,
					    State->Cubes[BIndex].Position);
	  }break;
	case SS:
	  {
	    EvolutionResult = EvolveSimplex(State, State->Spheres[AIndex].Position,
					    State->Spheres[BIndex].Position);
	  }break;
	case SF:
	  {
	    EvolutionResult = EvolveSimplex(State, State->Spheres[AIndex].Position,
					    State->Floor.Position);
	  }break;
	case FS:
	  {
	    EvolutionResult = EvolveSimplex(State, State->Floor.Position,
					    State->Spheres[BIndex].Position);
	  }break;
	}
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

internal void
Barycentric(glm::vec3 P, glm::vec3 A, glm::vec3 B, glm::vec3 C,
	    real32* U, real32* V, real32* W)
{
  glm::vec3 V0 = B - A, V1 = C - A, V2 = P - A;
  real32 D00 = glm::dot(V0, V0);
  real32 D01 = glm::dot(V0, V1);
  real32 D11 = glm::dot(V1, V1);
  real32 D20 = glm::dot(V2, V0);
  real32 D21 = glm::dot(V2, V1);
  real32 Denom = D00 * D11 - D01 * D01;
  *V = (D11 * D20 - D01 * D21) / Denom;
  *W = (D00 * D21 - D01 * D20) / Denom;
  *U = 1.0f - *V - *W;
}

internal void
DrawTriangles(sdl_state* State, sdl_render* Render)
{
  for(uint32 TriangleIndex = 0;
      TriangleIndex < State->Triangle->Count;
      ++TriangleIndex)
    {
      real32 vertices[] =
	{
	 State->Triangle->A[TriangleIndex].P.x, State->Triangle->A[TriangleIndex].P.y, State->Triangle->A[TriangleIndex].P.z,    
	 State->Triangle->B[TriangleIndex].P.x, State->Triangle->B[TriangleIndex].P.y, State->Triangle->B[TriangleIndex].P.z,
	 State->Triangle->C[TriangleIndex].P.x, State->Triangle->C[TriangleIndex].P.y, State->Triangle->C[TriangleIndex].P.z
	};
      glBindBuffer(GL_ARRAY_BUFFER, Render->VAOs[3]);
      glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_DYNAMIC_DRAW);
      glm::mat4 Model = glm::mat4(1.0);
      SetUniformM4(Render->Shaders[2], "Model", Model);
      glm::vec3 LineColor = glm::vec3(1.0, (real32)(TriangleIndex % 3)/2.0, (real32)(TriangleIndex % 4)/3.0);
      SetUniformV3(Render->Shaders[2], "LineColor", LineColor);
      glBindVertexArray(Render->VAOs[3]);
      glDrawArrays(GL_TRIANGLES, 0, 3);
    }
}

internal int32
ResolveCollision(sdl_state* State, sdl_input* Input, contact_info* Info)
{
  int32 CurrIter = 0;
  simplex* Simplex = State->Simplex;
  edge* Edge = State->Edge;
  triangle* Triangle = State->Triangle;
  
  // NOTE(Jovan): Take over points from GJK and construct a tetrahedron
  vertex A = Simplex->Vertices[0];
  vertex B = Simplex->Vertices[1];
  vertex C = Simplex->Vertices[2];
  vertex D = Simplex->Vertices[3];
  PushTriangle(Triangle, A, B, C); // ABC
  PushTriangle(Triangle, A, C, D); // ACD
  PushTriangle(Triangle, A, D, B); // ADB
  PushTriangle(Triangle, B, D, C); // BDC
  //Assert(Simplex->Count >= 4);
  while(CurrIter++ <= MAX_EPA_ITERATIONS)
    {
      // NOTE(Jovan): Find the closest triangle
      uint32 ClosestIndex = 0;
      real32 CurrentDistance = glm::abs(glm::dot(Triangle->N[0], Triangle->A[0].P));
      for(uint32 TriangleIndex = 0;
	  TriangleIndex < Triangle->Count;
	  ++TriangleIndex)
	{
	  // TODO(Jovan): Sanity check
	  glm::vec3 AB = Triangle->B[TriangleIndex].P - Triangle->A[TriangleIndex].P;
	  glm::vec3 AC = Triangle->C[TriangleIndex].P - Triangle->A[TriangleIndex].P;
	  glm::vec3 Normal = glm::normalize(glm::cross(AB, AC));
	  real32 Distance = glm::abs(glm::dot(Normal, Triangle->A[TriangleIndex].P));
	  if(Distance < CurrentDistance)
	    {
	      CurrentDistance = Distance;
	      ClosestIndex = TriangleIndex;
	    }
	}
      
      //CurrentDistance = CurrentDistance == 0 ? 0.01 : CurrentDistance;
      
      glm::vec3 Direction = Triangle->N[ClosestIndex];
      vertex NewSupport = CalculateSupport(State, Direction);


      // NOTE(Jovan): Calculate collision point and normal as linear combination
      // of barycentric pointss
      if(glm::dot(Triangle->N[ClosestIndex], NewSupport.P) - CurrentDistance < MAX_EPA_ERROR)
	{
	  real32 BaryU, BaryV, BaryW;
	  Barycentric(Triangle->N[ClosestIndex] * CurrentDistance,
	  	      Triangle->A[ClosestIndex].P,
	  	      Triangle->B[ClosestIndex].P,
	  	      Triangle->C[ClosestIndex].P,
	  	      &BaryU, &BaryV, &BaryW);
	  glm::vec3 CollisionPoint = ((BaryU * Triangle->A[ClosestIndex].SupA) +
				      (BaryV * Triangle->B[ClosestIndex].SupA) +
				      (BaryW * Triangle->C[ClosestIndex].SupA));
	  glm::vec3 CollisionNormal = -1.0f * Triangle->N[ClosestIndex];
	  real32 Depth = CurrentDistance;
	  State->Depth = Depth;

	  Info->Point = CollisionPoint;
	  // TODO(Jovan): Probably not needed
	  Info->Normal = CollisionNormal;
	  return ClosestIndex;
	}
      // NOTE(Jovan): Removing triangle that can be "seen" by the new point
      for(uint32 TriangleIndex = 0;
	  TriangleIndex < Triangle->Count;)
	{
	  // TODO(Jovan): Sanity check
	  glm::vec3 Temp = NewSupport.P - Triangle->A[TriangleIndex].P;
	  glm::vec3 AB = Triangle->B[TriangleIndex].P - Triangle->A[TriangleIndex].P;
	  glm::vec3 AC = Triangle->C[TriangleIndex].P - Triangle->A[TriangleIndex].P;
	  glm::vec3 Normal = glm::normalize(glm::cross(AB, AC));
	  if(glm::dot(Normal, Triangle->A[TriangleIndex].P) < 0)
	    {
	      Normal *= -1.0f;
	    }
	  if(glm::dot(Normal, Temp) > 0)
	    {
	      // NOTE(Jovan): "Disolve" the removed triangle into it's edges
	      // and push them onto the edge list
	      PushEdge(Edge, Triangle->A[TriangleIndex], Triangle->B[TriangleIndex]); // AB
	      PushEdge(Edge, Triangle->B[TriangleIndex], Triangle->C[TriangleIndex]); // BC
	      PushEdge(Edge, Triangle->C[TriangleIndex], Triangle->A[TriangleIndex]); // CA
	      RemoveTriangle(Triangle, TriangleIndex);
	      continue;
	    }
	  ++TriangleIndex;
	}

      // NOTE(Jovan): Construct new triangles
      for(uint32 EdgeIndex = 0;
	  EdgeIndex < Edge->Count;
	  ++EdgeIndex)
	{
	  PushTriangle(Triangle, NewSupport, Edge->A[EdgeIndex], Edge->B[EdgeIndex]);
	}
      // NOTE(Jovan): Clear the edge list
      ClearEdges(Edge);
    }
  return -1;
}

internal void
IntegrateForces(sdl_state* State, uint32 CubeIndex, real32 dt)
{
  phys_return Y0 = {};
  phys_return Y = {};
  // NOTE(Jovan): Linear
  Y0.X = State->Cubes[CubeIndex].Position;
  Y0.Y = State->Cubes[CubeIndex].V;
  Y = Euler(MovementFunction, dt, Y0, State->Cubes[CubeIndex].Forces, State->Cubes[CubeIndex].Mass);
  State->Cubes[CubeIndex].V = Y.Y;

  // NOTE(Jovan): Rotational
  Y0.X = State->Cubes[CubeIndex].Angles;
  Y0.Y = State->Cubes[CubeIndex].W;
  Y = Euler(RotationFunction, dt, Y0, State->Cubes[CubeIndex].Torque, State->Cubes[CubeIndex].MOI);
  State->Cubes[CubeIndex].W = Y.Y;
  
  // NOTE(Jovan): Update states and clear the forces/torques
  UpdateVertices(State, CubeIndex);
  CubeClearForces(State, CubeIndex);
}

internal void
Constraint(sdl_state* State, glm::vec3 N,
	   contact_info* InfoA, contact_info* InfoB, real32 dt)
{
  if(N == glm::vec3(0.0))
    {
      return;
    }

  // NOTE(Jovan): Calculating J(M^-1)(J^T)
  // ------------------------------------
  glm::vec3 R1, R2;
  glm::vec3 PosA, PosB;
  glm::vec3 V1, W1, V2, W2;
  // NOTE(Jovan): I^-1
  real32 InvI1, InvI2;
  // NOTE(Jovan): M^-1
  real32 InvM1, InvM2;
  switch(State->CurrentCollisionType)
    {
    case CC:
      {
	PosA = State->Cubes[State->IndexA].Position;
        PosB = State->Cubes[State->IndexB].Position;
	InvI1 = 1.0f / State->Cubes[State->IndexA].MOI;
	InvI2 = 1.0f / State->Cubes[State->IndexB].MOI;
	InvM1 = 1.0f / State->Cubes[State->IndexA].Mass;
	InvM2 = 1.0f /State->Cubes[State->IndexB].Mass;
	V1 = State->Cubes[State->IndexA].V;
	V2 = State->Cubes[State->IndexB].V;
	W1 = State->Cubes[State->IndexA].W;
	W2 = State->Cubes[State->IndexB].W;
      }break;
    case CS:
      {
	PosA = State->Cubes[State->IndexA].Position;
        PosB = State->Spheres[State->IndexB].Position;
	InvI1 = 1.0f / State->Cubes[State->IndexA].MOI;
	InvI2 = 1.0f / State->Spheres[State->IndexB].MOI;
	InvM1 = 1.0f / State->Cubes[State->IndexA].Mass;
	InvM2 = 1.0f / State->Spheres[State->IndexB].Mass;
	V1 = State->Cubes[State->IndexA].V;
	V2 = State->Spheres[State->IndexB].V;
	W1 = State->Cubes[State->IndexA].W;
	W2 = State->Spheres[State->IndexB].W;
      }break;
    case SC:
      {
	PosA = State->Spheres[State->IndexA].Position;
        PosB = State->Cubes[State->IndexB].Position;
	InvI1 = 1.0f / State->Spheres[State->IndexA].MOI;
	InvI2 = 1.0f / State->Cubes[State->IndexB].MOI;
	InvM1 = 1.0f / State->Spheres[State->IndexA].Mass;
	InvM2 = 1.0f / State->Cubes[State->IndexB].Mass;
	V1 = State->Spheres[State->IndexA].V;
	V2 = State->Cubes[State->IndexB].V;
	W1 = State->Spheres[State->IndexA].W;
	W2 = State->Cubes[State->IndexB].W;
      }break;
    case CF:
      {
	PosA = State->Cubes[State->IndexA].Position;
	PosB = State->Floor.Position;
	InvI1 = 1.0f / State->Cubes[State->IndexA].MOI;
	InvI2 = 1.0f / State->Floor.MOI;
	InvM1 = 1.0f / State->Cubes[State->IndexA].Mass;
	InvM2 = 1.0f / State->Floor.Mass;
	V1 = State->Cubes[State->IndexA].V;
	V2 = State->Floor.V;
	W1 = State->Cubes[State->IndexA].W;
	W2 = State->Floor.W;
      }break;
    case FC:
      {
	PosA = State->Floor.Position;
	PosB = State->Cubes[State->IndexB].Position;
	InvI1 = 1.0f / State->Floor.MOI;
	InvI2 = 1.0f / State->Cubes[State->IndexB].MOI;
	InvM1 = 1.0f / State->Floor.Mass;
	InvM2 = 1.0f / State->Cubes[State->IndexB].Mass;
	V1 = State->Floor.V;
	V2 = State->Cubes[State->IndexB].V;
	W1 = State->Floor.W;
	W2 = State->Cubes[State->IndexB].W;
      }break;
    case SS:
      {
	PosA = State->Spheres[State->IndexA].Position;
	PosB = State->Spheres[State->IndexB].Position;
	InvI1 = 1.0f / State->Spheres[State->IndexA].MOI;
	InvI2 = 1.0f / State->Spheres[State->IndexB].MOI;
	InvM1 = 1.0f / State->Spheres[State->IndexA].Mass;
	InvM2 = 1.0f / State->Spheres[State->IndexB].Mass;
	V1 = State->Spheres[State->IndexA].V;
	V2 = State->Spheres[State->IndexB].V;
	W1 = State->Spheres[State->IndexA].W;
	W2 = State->Spheres[State->IndexB].W;
      }break;
    case SF:
      {
	PosA = State->Spheres[State->IndexA].Position;
	PosB = State->Floor.Position;
	InvI1 = 1.0f / State->Spheres[State->IndexA].MOI;
	InvI2 = 1.0f / State->Floor.MOI;
	InvM1 = 1.0f / State->Spheres[State->IndexA].Mass;
	InvM2 = 1.0f / State->Floor.Mass;
	V1 = State->Spheres[State->IndexA].V;
	V2 = State->Floor.V;
	W1 = State->Spheres[State->IndexA].W;
	W2 = State->Floor.W;
      }break;
    case FS:
      {
	PosA = State->Floor.Position;
	PosB = State->Spheres[State->IndexB].Position;
	InvI1 = 1.0f / State->Floor.MOI;
	InvI2 = 1.0f / State->Spheres[State->IndexB].MOI;
	InvM1 = 1.0f / State->Floor.Mass;
	InvM2 = 1.0f / State->Spheres[State->IndexB].Mass;
	V1 = State->Floor.V;
	V2 = State->Spheres[State->IndexB].V;
	W1 = State->Floor.W;
	W2 = State->Spheres[State->IndexB].W;
      }break;
    default:
      {
	printf("ERROR::CONSTRAINT::Unknown collision type\n");
	return;
      }break;
    }
  R1 = InfoA->Point - PosA;
  R2 = InfoB->Point - PosB;
  
  // NOTE(Jovan): Penetration depth
  real32 Depth = glm::dot(((PosA + R1) - (PosB + R2)), N);
  
  // NOTE(Jovan): R x N
  glm::vec3 RN1 = glm::cross(R1, N);
  glm::vec3 RN2 = glm::cross(R2, N);
  // NOTE(Jovan): JMJ
  real32 JMJ = (InvM1) + (InvM2);
  JMJ += InvI1 * (glm::dot(RN1, RN1)) - InvI2 * (glm::dot(-RN2, -RN2));
  // NOTE(Jovan): 1 / (JMJ) (Effective mass)
  JMJ = 1.0f / JMJ;
  
  // NOTE(Jovan): Baumgarte
  real32 Beta = 0.2f;
  // NOTE(Jovan): Bias
  real32 B = Beta/dt * Depth;
  // TODO(Jovan): Restitution bias
  // TODO(Jovan): Friction?
  int32 Iter = 32;
  real32 Impulse = 0;
  while(Iter--)
    {
      // real32 Projection = glm::dot((V1 + glm::cross(W1, R1) - V2 - glm::cross(W2, R2)), N);
      // if(Projection < 0)
      // 	{
      // 	  break;
      // 	}
      // NOTE(Jovan): Derivative of V
      glm::vec3 dV = V1 + glm::cross(W1, N) - V2 - glm::cross(W2, N);
      real32 JdV = glm::dot(dV, N);
      real32 Lambda = -(JdV + B) * JMJ;

      // NOTE(Jovan): Accumulating and clamping impulse
      real32 OldAccumI = State->AccumI;
      State->AccumI += Lambda;
      if(State->AccumI < 0)
	{
	  State->AccumI = 0.0f;
	}

      Impulse = State->AccumI - OldAccumI;

    }
  // NOTE(Jovan): Calculating linear impulse
  glm::vec3 LinearImpulse = N * Impulse;
  // NOTE(Jovan): Calculating angular impulses
  glm::vec3 AngularI1 = RN1 * Impulse;
  glm::vec3 AngularI2 = RN2 * Impulse;

  // NOTE(Jovan): Applying linear and angular impulses
  switch(State->CurrentCollisionType)
    {
    case CC:
      {
	State->Cubes[State->IndexA].V += InvM1 * LinearImpulse;
	State->Cubes[State->IndexB].V -= InvM2 * LinearImpulse;
	State->Cubes[State->IndexA].W += InvI1 * AngularI1;
	State->Cubes[State->IndexB].W -= InvI2 * AngularI2;
      }break;
    case CS:
      {
	State->Cubes[State->IndexA].V += InvM1 * LinearImpulse;
	State->Spheres[State->IndexB].V -= InvM2 * LinearImpulse;
	State->Cubes[State->IndexA].W += InvI1 * AngularI1;
	State->Spheres[State->IndexB].W -= InvI2 * AngularI2;
      }break;
    case SC:
      {
	State->Spheres[State->IndexA].V += InvM1 * LinearImpulse;
	State->Cubes[State->IndexB].V -= InvM2 * LinearImpulse;
	State->Spheres[State->IndexA].W += InvI1 * AngularI1;
	State->Cubes[State->IndexB].W -= InvI2 * AngularI2;
      }break;
    case CF:
      {
	// NOTE(Jovan): Not updating floor velocities
	State->Cubes[State->IndexA].V += InvM1 * LinearImpulse;
	State->Cubes[State->IndexA].W += InvI1 * AngularI1;
      }break;
    case FC:
      {
	// NOTE(Jovan): Not updating floor velocities
	State->Cubes[State->IndexB].V -= InvM2 * LinearImpulse;
	State->Cubes[State->IndexB].W -= InvI2 * AngularI2;
      }break;
    case SS:
      {
	State->Spheres[State->IndexA].V += InvM1 * LinearImpulse;
	State->Spheres[State->IndexB].V -= InvM2 * LinearImpulse;
	State->Spheres[State->IndexA].W += InvI1 * AngularI1;
	State->Spheres[State->IndexB].W -= InvI2 * AngularI2;
      }break;
    case SF:
      {
	// NOTE(Jovan): Not updating floor velocities
	State->Spheres[State->IndexA].V += InvM1 * LinearImpulse;
	State->Spheres[State->IndexA].W += InvI1 * AngularI1;
      }break;
    case FS:
      {
	// NOTE(Jovan): Not updating floor velocities
	State->Cubes[State->IndexB].V -= InvM2 * LinearImpulse;
	State->Cubes[State->IndexB].W -= InvI2 * AngularI2;
      }break;
    default:
      {
	printf("ERROR::CONSTRAINT::Unknown collision type\n");
	return;
      }break;
    }
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
      InitializeArena(&SimState->SimplexArena, Memory->PermanentStorageSize - sizeof(sdl_state),
		      (uint8*)Memory->PermanentStorage + sizeof(sdl_state));
      SimState->Simplex = PushStruct(&SimState->SimplexArena, simplex);
      simplex* Simplex = SimState->Simplex;
      Simplex->Count = 0;
      Simplex->Vertices = PushArray(&SimState->SimplexArena, 64, vertex);
      SimState->GJKIteration = 0;

      InitializeArena(&SimState->EdgeArena, Memory->TransientStorageSize / 2,
		      (uint8*)Memory->TransientStorage);
      SimState->Edge = PushStruct(&SimState->EdgeArena, edge);
      edge* Edge = SimState->Edge;
      Edge->Count = 0;
      Edge->A = PushArray(&SimState->EdgeArena, 32, vertex);
      Edge->B = PushArray(&SimState->EdgeArena, 32, vertex);

      InitializeArena(&SimState->TriangleArena, Memory->TransientStorageSize / 2,
		      (uint8*)Memory->TransientStorage + Memory->TransientStorageSize / 2);
      SimState->Triangle = PushStruct(&SimState->TriangleArena, triangle);
      triangle* Triangle = SimState->Triangle;
      Triangle->Count = 0;
      Triangle->A = PushArray(&SimState->TriangleArena, 128, vertex);
      Triangle->B = PushArray(&SimState->TriangleArena, 128, vertex);
      Triangle->C = PushArray(&SimState->TriangleArena, 128, vertex);
      Triangle->N = PushArray(&SimState->TriangleArena, 128, glm::vec3);
      
      // NOTE(Jovan): Camera init
      SimState->Camera.FOV = 45.0f; 
      SimState->Camera.Pitch = 0.0f;
      SimState->Camera.Yaw = -90.0f;
      SimState->Camera.Speed = 0.05f;
      
      SimState->Camera.Position = glm::vec3(0.0f, 0.0f, 3.0f);
      glm::vec3 Up = glm::vec3(0.0f, 1.0f, 0.0f);
      SimState->Camera.Target = glm::vec3(0.0f, 0.0f, 0.0f);
      SimState->Camera.Direction = glm::normalize(SimState->Camera.Position - SimState->Camera.Target);
      SimState->Camera.Front = glm::vec3(0.0f, 0.0f, -1.0f);
      SimState->Camera.Right = glm::normalize(glm::cross(Up, SimState->Camera.Direction));
      SimState->Camera.Up = glm::cross(SimState->Camera.Direction, SimState->Camera.Right);
      
      // NOTE(Jovan): Cube init
      SimState->CubeCount = 3;
      SimState->Cubes[0].Model = glm::mat4(1.0);
      UpdateVertices(SimState, 0);
      SimState->Cubes[0].Position = glm::vec3(2.0, 4.5f, 2.1);//(2.1, 2.5 ,2.1);
      SimState->Cubes[0].V = glm::vec3(0.0);
      SimState->Cubes[0].Forces = glm::vec3(0.0);
      SimState->Cubes[0].Angles = glm::vec3(0.0);
      SimState->Cubes[0].W = glm::vec3(0.0);
      SimState->Cubes[0].Torque = glm::vec3(0.0);
      SimState->Cubes[0].Size = 1.0f;
      SimState->Cubes[0].Mass = 1.0f;
      SimState->Cubes[0].MOI = (SimState->Cubes[0].Mass / 12.0f) *
      	(2.0f * SimState->Cubes[0].Size * SimState->Cubes[0].Size);

      SimState->Cubes[1].Model = glm::mat4(1.0);
      UpdateVertices(SimState, 1);
      SimState->Cubes[1].Position = glm::vec3(2.1, 3.5, 2.0);
      SimState->Cubes[1].V = glm::vec3(0.0);
      SimState->Cubes[1].Forces = glm::vec3(0.0);
      SimState->Cubes[1].Angles = glm::vec3(0.0);
      SimState->Cubes[1].W = glm::vec3(0.0);
      SimState->Cubes[1].Torque = glm::vec3(0.0);
      SimState->Cubes[1].Size = 1.0f;
      SimState->Cubes[1].Mass = 10.0f;
      SimState->Cubes[1].MOI = (SimState->Cubes[1].Mass / 12.0f) *
	(2.0f * SimState->Cubes[1].Size * SimState->Cubes[1].Size);
      
      SimState->Cubes[2].Model = glm::mat4(1.0);
      UpdateVertices(SimState, 1);
      SimState->Cubes[2].Position = glm::vec3(2.0, 1.5, 2.0);
      SimState->Cubes[2].V = glm::vec3(0.0);
      SimState->Cubes[2].Forces = glm::vec3(0.0);
      SimState->Cubes[2].Angles = glm::vec3(0.0);
      SimState->Cubes[2].W = glm::vec3(0.0);
      SimState->Cubes[2].Torque = glm::vec3(0.0);
      SimState->Cubes[2].Size = 1.0f;
      SimState->Cubes[2].Mass = 10.0f;
      SimState->Cubes[2].MOI = (SimState->Cubes[2].Mass / 12.0f) *
	(2.0f * SimState->Cubes[2].Size * SimState->Cubes[2].Size);

      // NOTE(Jovan): Sphere init
      SimState->SphereCount = 2;
      SimState->Spheres[0].Model = glm::mat4(1.0f);
      SimState->Spheres[0].Position = glm::vec3(3.0f, 4.0f, 2.0f);
      SimState->Spheres[0].V = glm::vec3(0.0f);
      SimState->Spheres[0].Forces = glm::vec3(0.0f);
      SimState->Spheres[0].Angles = glm::vec3(0.0f);
      SimState->Spheres[0].W = glm::vec3(0.0f);
      SimState->Spheres[0].Torque = glm::vec3(0.0f);
  
      SimState->Spheres[0].Radius = 0.025f;
      SimState->Spheres[0].Mass = 10.0f;
      SimState->Spheres[0].MOI = (2.0f / 5.0f) * SimState->Spheres[0].Mass *
	pow(SimState->Spheres[0].Radius, 2);
      
      SimState->Spheres[1].Model = glm::mat4(1.0f);
      SimState->Spheres[1].Position = glm::vec3(3.0f, 4.0f, 2.0f);
      SimState->Spheres[1].V = glm::vec3(0.0f);
      SimState->Spheres[1].Forces = glm::vec3(0.0f);
      SimState->Spheres[1].Angles = glm::vec3(0.0f);
      SimState->Spheres[1].W = glm::vec3(0.0f);
      SimState->Spheres[1].Torque = glm::vec3(0.0f);
  
      SimState->Spheres[1].Radius = 0.025f;
      SimState->Spheres[1].Mass = 10.0f;
      SimState->Spheres[1].MOI = (2.0f / 5.0f) * SimState->Spheres[1].Mass *
	pow(SimState->Spheres[1].Radius, 2);

      // NOTE(Jovan): Floor init
      SimState->Floor.Model = glm::mat4(1.0);
      FloorUpdateVertices(SimState);
      SimState->Floor.Position = glm::vec3(0.0f, -0.5f, 0.0f);
      SimState->Floor.V = glm::vec3(0.0);
      SimState->Floor.Forces = glm::vec3(0.0);
      SimState->Floor.Angles = glm::vec3(0.0);
      SimState->Floor.W = glm::vec3(0.0);
      SimState->Floor.Torque = glm::vec3(0.0);
      SimState->Floor.Size = 1.0f;
      SimState->Floor.Mass = 1.0f;
      SimState->Floor.MOI = (SimState->Floor.Mass / 12.0f) *
      	(2.0f * SimState->Floor.Size * SimState->Floor.Size);
      
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
  
  // NOTE(Jovan): Input
  // ------------------
  HandleInput(SimState, Input, dt);
  UpdateCamera(SimState, Input);
  // NOTE(Jovan) End input
  // ---------------------

  // NOTE(Jovan): Physics stuff
  // --------------------------
  int32 Closest = -1;
  contact_info A = {};
  contact_info B = {};
  glm::vec3 Normal;
  bool32 CollisionHappened = 0;
  // TODO IMPORTANT (Jovan): Generalize collision AND STABILIZE PLS
  for(uint32 CubeIndex1 = 0;
      CubeIndex1 < SimState->CubeCount;
      ++CubeIndex1)
    {
      IntegrateForces(SimState, CubeIndex1, dt);
      for(uint32 CubeIndex2 = CubeIndex1;
	  CubeIndex2 < SimState->CubeCount;
	  ++CubeIndex2)
	{
	  if(CubeIndex1 == CubeIndex2)
	    {
	      continue;
	    }
	  SimState->IndexA = CubeIndex2;
	  SimState->IndexB = CubeIndex1;
	  SimState->CurrentCollisionType = CC;
	  CollisionHappened = DetectCollisions(SimState);
	  Closest = -1;
	  A = {};
	  B = {};
	  IntegrateForces(SimState, CubeIndex2, dt);
	  if(CollisionHappened)
	    {
	      // NOTE(Jovan): For B to A
	      Closest = ResolveCollision(SimState, Input, &B);

	      // NOTE(Jovan): Clear B to A collision info so GJK can start over
	      // for A to B
	      ClearTriangles(SimState->Triangle);
	      ClearEdges(SimState->Edge);
	      ClearVertices(SimState->Simplex);

	      // NOTE(Jovan): For A to B
	      SimState->IndexA = CubeIndex1;
	      SimState->IndexB = CubeIndex2;
	      DetectCollisions(SimState);
	      Closest = ResolveCollision(SimState, Input, &A);
	      ClearTriangles(SimState->Triangle);
	      ClearEdges(SimState->Edge);
	      ClearVertices(SimState->Simplex);
	      
	      Normal = -SimState->Triangle->N[Closest];
	      Constraint(SimState, Normal, &A, &B, dt);
	      
	    }
	}
    }
#if 0
  ClearTriangles(SimState->Triangle);
  ClearEdges(SimState->Edge);
  ClearVertices(SimState->Simplex);
  // NOTE(Jovan): Floor collisions
  FloorUpdateVertices(SimState);
  SimState->IndexA = 0;
  SimState->IndexB = 0;
  SimState->CurrentCollisionType = FC;
  A = {};
  B = {};
  Normal = glm::vec3(0.0);
  bool32 CollisionFloor = DetectCollisions(SimState);
  if(CollisionFloor)
    {
      printf("FLOOR COLLISION!\n");
      ResolveCollision(SimState, Input, &B);
      ClearTriangles(SimState->Triangle);
      ClearEdges(SimState->Edge);
      ClearVertices(SimState->Simplex);

      SimState->IndexA = 0;
      SimState->IndexB = 0;
      SimState->CurrentCollisionType = CF;
      DetectCollisions(SimState);
      Closest = ResolveCollision(SimState, Input, &A);
      ClearTriangles(SimState->Triangle);
      ClearEdges(SimState->Edge);
      ClearVertices(SimState->Simplex);
      
      Normal = -SimState->Triangle->N[Closest];
      printf("Happening!\n");
      Constraint(SimState, Normal, &A, &B, dt);
      
    }
#endif
  // CubeAddForce(SimState, 0, SimState->Cubes[0].Mass * GRAVITY_ACCEL * glm::vec3(0.0, -1.0, 0.0));
  // CubeAddForce(SimState, 1, SimState->Cubes[1].Mass * GRAVITY_ACCEL * glm::vec3(0.0, -1.0, 0.0));
  //CubeAddTorque(SimState,0, glm::vec3(1.0));
  
  // TODO(Jovan): Integrate velocities here?
  for(uint32 CubeIndex = 0;
      CubeIndex < SimState->CubeCount;
      ++CubeIndex)
    {
      SimState->Cubes[CubeIndex].Position += dt * SimState->Cubes[CubeIndex].V;
      SimState->Cubes[CubeIndex].Angles += dt * SimState->Cubes[CubeIndex].W;
      UpdateVertices(SimState, CubeIndex);
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

#if DRAW_WIRE
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
#else
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
#endif
  // NOTE(Jovan): Collision depth drawing
  glUseProgram(Render->Shaders[2]);
  SetUniformM4(Render->Shaders[2], "View", View);
  SetUniformM4(Render->Shaders[2], "Projection", Projection);
  if(Closest != -1)
    {
      real32 PenDepth = glm::dot(((SimState->Cubes[0].Position + A.Point) -
				    (SimState->Cubes[1].Position + B.Point)),
				   Normal);
      glm::vec3 LineColor = glm::vec3(1.0, 0.0, 0.0);
      real32 Vertices[] =
	{
	 A.Point.x, A.Point.y, A.Point.z,
	 B.Point.x, B.Point.y, B.Point.z 
	};
      glBindBuffer(GL_ARRAY_BUFFER, Render->VAOs[3]);
      glBufferData(GL_ARRAY_BUFFER, sizeof(Vertices), Vertices, GL_DYNAMIC_DRAW);
      Model = glm::mat4(1.0);
      SetUniformM4(Render->Shaders[2], "Model", Model);
      SetUniformV3(Render->Shaders[2], "LineColor", LineColor);
      glBindVertexArray(Render->VAOs[3]);
      glDrawArrays(GL_LINE_STRIP, 0, 2);
    }
  
  // NOTE(Jovan): Coordinate drawing
  // -------------------------------
  glUseProgram(Render->Shaders[2]);
  SetUniformM4(Render->Shaders[2], "View", View);
  SetUniformM4(Render->Shaders[2], "Projection", Projection);
  real32 LineLength = 100.0f;
  glm::vec3 LineColor = glm::vec3(0.0);

  // NOTE(Jovan): Drawing polytope
  // ----------------
#if DRAW_EPA
  DrawTriangles(SimState, Render);

  if(Closest != -1)
    {
      real32 X = -SimState->Triangle->N[Closest].x;
      real32 Y = -SimState->Triangle->N[Closest].y;
      real32 Z = -SimState->Triangle->N[Closest].z;
      real32 vertices[] =
	{
	 0, 0, 0,
	 X, Y, Z
	};
      glBindBuffer(GL_ARRAY_BUFFER, Render->VAOs[3]);
      glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_DYNAMIC_DRAW);
      Model = glm::mat4(1.0);
      LineColor = glm::vec3(1.0f, 0.0f, 1.0f);
      SetUniformM4(Render->Shaders[2], "Model", Model);
      SetUniformV3(Render->Shaders[2], "LineColor", LineColor);
      glBindVertexArray(Render->VAOs[3]);
      glDrawArrays(GL_LINE_STRIP, 0, 2);
    }
  // NOTE(Jovan): End test
  // ---------------------
#endif
#if DRAW_COORDINATES
  for(uint32 CubeIndex = 0;
      CubeIndex < SimState->CubeCount;
      ++CubeIndex)
    {
      for(uint32 i = 0;
	  i < 8;
	  ++i)
	{
	  
	  // NOTE(Jovan): X
	  Model = glm::mat4(1.0);
	  Model = glm::translate(Model, SimState->Cubes[CubeIndex].Vertices[i]);
	  Model = glm::scale(Model, glm::vec3(LineLength));
	  SetUniformM4(Render->Shaders[2], "Model", Model);
	  LineColor = glm::vec3(1.0f, 0.0f, 0.0f);
	  SetUniformV3(Render->Shaders[2], "LineColor", LineColor);
	  glBindVertexArray(Render->VAOs[3]);
	  glDrawArrays(GL_LINES, 0, 6);
	  // NOTE(Jovan): Y
	  Model = glm::mat4(1.0);
	  Model = glm::translate(Model, SimState->Cubes[CubeIndex].Vertices[i]);
	  Model = glm::rotate(Model, glm::radians(90.0f), glm::vec3(0.0f, 0.0f, 1.0f));
	  Model = glm::scale(Model, glm::vec3(LineLength));
	  SetUniformM4(Render->Shaders[2], "Model", Model);
	  LineColor = glm::vec3(0.0f, 1.0f, 0.0f);
	  SetUniformV3(Render->Shaders[2], "LineColor", LineColor);
	  glBindVertexArray(Render->VAOs[3]);
	  glDrawArrays(GL_LINES, 0, 6);
	  // NOTE(Jovan): Z
	  Model = glm::mat4(1.0);
	  Model = glm::translate(Model, SimState->Cubes[CubeIndex].Vertices[i]);
	  Model = glm::rotate(Model, glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
	  Model = glm::scale(Model, glm::vec3(LineLength));
	  SetUniformM4(Render->Shaders[2], "Model", Model);
	  LineColor = glm::vec3(0.0f, 0.0f, 1.0f);
	  SetUniformV3(Render->Shaders[2], "LineColor", LineColor);
	  glBindVertexArray(Render->VAOs[3]);
	  glDrawArrays(GL_LINES, 0, 6);
	}
    }
#endif
  glBindVertexArray(0);

#if DRAW_FLOOR
  // NOTE(Jovan): Floor drawing
  // --------------------------

  glUseProgram(Render->Shaders[0]);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, Render->Textures[2]);
  SetUniformM4(Render->Shaders[0], "View", View);
  SetUniformM4(Render->Shaders[0], "Projection", Projection);
  SimState->Floor.Model = glm::mat4(1.0);
  SimState->Floor.Model = glm::translate(SimState->Floor.Model, SimState->Floor.Position);
  SimState->Floor.Model = glm::rotate(SimState->Floor.Model, glm::radians(SimState->Floor.Angles.x),
				      glm::vec3(1.0f, 0.0f, 0.0f));
  SimState->Floor.Model = glm::rotate(SimState->Floor.Model, glm::radians(SimState->Floor.Angles.y),
				      glm::vec3(0.0f, 1.0f, 0.0f));
  SimState->Floor.Model = glm::rotate(SimState->Floor.Model, glm::radians(SimState->Floor.Angles.z),
				      glm::vec3(0.0f, 0.0f, 1.0f));
  SimState->Floor.Model = glm::scale(SimState->Floor.Model,
				     glm::vec3(SimState->Floor.Size));
  SetUniformM4(Render->Shaders[0], "Model", SimState->Floor.Model);
  FloorUpdateVertices(SimState);
  glBindVertexArray(Render->VAOs[0]);
  glDrawArrays(GL_TRIANGLES, 0, 36);
			 
  // NOTE(Jovan): End floor drawing
  // ------------------------------
  SetUniformM4(Render->Shaders[0], "View", View);
  SetUniformM4(Render->Shaders[0], "Projection", Projection);
#endif 

  // NOTE(Jovan): Minkowski sum drawing
  // ---------------------------------

#if DRAW_MINKOWSKI
   SimState->IndexA = 0;
   SimState->IndexB = 0;
   SimState->CurrentCollisionType = CF;
  vertex v1 = CalculateSupport(SimState, glm::vec3(1.0, 1.0, 1.0));
  vertex v2 = CalculateSupport(SimState, glm::vec3(-1.0, 1.0, 1.0));
  vertex v3 = CalculateSupport(SimState, glm::vec3(1.0, -1.0, 1.0));
  vertex v4 = CalculateSupport(SimState, glm::vec3(-1.0, -1.0, 1.0));
  vertex v5 = CalculateSupport(SimState, glm::vec3(1.0, 1.0, -1.0));
  vertex v6 = CalculateSupport(SimState, glm::vec3(-1.0, 1.0, -1.0));
  vertex v7 = CalculateSupport(SimState, glm::vec3(1.0, -1.0, -1.0));
  vertex v8 = CalculateSupport(SimState, glm::vec3(-1.0, -1.0, -1.0));
  real32 MinkowskiVertices[] = {
				v1.P.x, v1.P.y, v1.P.z,
				v2.P.x, v2.P.y, v2.P.z,
				
				v1.P.x, v1.P.y, v1.P.z,
				v3.P.x, v3.P.y, v3.P.z,

				v3.P.x, v3.P.y, v3.P.z,
				v4.P.x, v4.P.y, v4.P.z,
				
				v5.P.x, v5.P.y, v5.P.z,
				v6.P.x, v6.P.y, v6.P.z,
				
				v5.P.x, v5.P.y, v5.P.z,
				v7.P.x, v7.P.y, v7.P.z,
				
				v7.P.x, v7.P.y, v7.P.z,
				v8.P.x, v8.P.y, v8.P.z,
  };
  glBindBuffer(GL_ARRAY_BUFFER, Render->VAOs[3]);
  glBufferData(GL_ARRAY_BUFFER, sizeof(MinkowskiVertices), MinkowskiVertices, GL_DYNAMIC_DRAW);
  Model = glm::mat4(1.0);
  SetUniformM4(Render->Shaders[2], "Model", Model);
  LineColor = glm::vec3(1.0f);
  SetUniformV3(Render->Shaders[2], "LineColor", LineColor);
  glBindVertexArray(Render->VAOs[3]);
  glDrawArrays(GL_LINES, 0, 36);  
#endif
  
  // NOTE(Jovan): End Minkowski sum drawing
  // --------------------------------------
  
  // NOTE(Jovan): Cube drawing
  // -------------------------
  glUseProgram(Render->Shaders[0]);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, Render->Textures[0]);

  SetUniformM4(Render->Shaders[0], "View", View);
  SetUniformM4(Render->Shaders[0], "Projection", Projection);
#if DRAW_CUBES
  for(uint32 CubeIndex = 0;
      CubeIndex < SimState->CubeCount;
      ++CubeIndex)
    {
      Model = glm::mat4(1.0f);
      Model = glm::translate(Model, SimState->Cubes[CubeIndex].Position);
      Model = glm::rotate(Model, glm::radians(SimState->Cubes[CubeIndex].Angles.x),
			  glm::vec3(1.0f, 0.0f, 0.0f));
      Model = glm::rotate(Model, glm::radians(SimState->Cubes[CubeIndex].Angles.y),
			  glm::vec3(0.0f, 1.0f, 0.0f));
      Model = glm::rotate(Model, glm::radians(SimState->Cubes[CubeIndex].Angles.z),
			  glm::vec3(0.0f, 0.0f, 1.0f));
      Model = glm::scale(Model, glm::vec3(SimState->Cubes[CubeIndex].Size,
					  SimState->Cubes[CubeIndex].Size,
					  SimState->Cubes[CubeIndex].Size));
      SimState->Cubes[CubeIndex].Model = Model;
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
#if DRAW_SPHERES
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
      Model = glm::rotate(Model, glm::radians(SimState->Spheres[SphereIndex].Angles.x),
			  glm::vec3(1.0f, 0.0f, 0.0f));
      Model = glm::rotate(Model, glm::radians(SimState->Spheres[SphereIndex].Angles.y),
			  glm::vec3(0.0f, 1.0f, 0.0f));
      Model = glm::rotate(Model, glm::radians(SimState->Spheres[SphereIndex].Angles.z),
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
#if LOGGING
  // TODO(Jovan): Make better logging
  printf("Floor coords:");
  PrintVector(SimState->Floor.Position);
  printf("Collision: %d\n", CollisionHappened);
  printf("Triangles: %d\n", SimState->Triangle->Count);
  printf("Normals:\n");
  PrintVector(SimState->Triangle->N[0]);
  printf("Cube 0 position:");
  PrintVector(SimState->Cubes[0].Position);
#endif
  // NOTE(Jovan): End logging
  // ------------------------
  ClearVertices(SimState->Simplex);
  ClearTriangles(SimState->Triangle);
  ClearEdges(SimState->Edge);
  // NOTE(Jovan): Clear accumulated impulse
  //SimState->AccumI = 0.0f;
}
