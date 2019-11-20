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
  // Triangle->N[Triangle->Count] = glm::cross(Triangle->A[Triangle->Count].P,
  // 					    Triangle->B[Triangle->Count].P);
  Triangle->N[Triangle->Count] = glm::triangleNormal(A.P, B.P, C.P);
  // NOTE(Jovan): Normalize vector
  real32 Len = sqrt(glm::dot(Triangle->N[Triangle->Count], Triangle->N[Triangle->Count]));
  Triangle->N[Triangle->Count] *= 1.0f/Len;
  glm::vec3 A0 = -1.0f * Triangle->A[Triangle->Count].P;
  if(glm::dot(A0, Triangle->N[Triangle->Count]) < 0)
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

internal vertex
CalculateSupport(sdl_state* State, int32 AIndex, int32 BIndex, glm::vec3 Direction, collision_type Type)
{
  vertex Result = {};
  switch(Type)
    {
    case CC:
      {
	glm::vec3 SupportA = GetCubeSupport(State, AIndex, Direction);
	glm::vec3 SupportB = GetCubeSupport(State, BIndex, -1.0f * Direction);
	Result.P = SupportA - SupportB;
	Result.SupA = SupportA;
	PushVertex(State->Simplex, Result);
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

internal bool32
AddSupport(sdl_state* State, int32 AIndex, int32 BIndex, glm::vec3 Direction)
{
  bool32 Result = 0;
  vertex NewVertex = CalculateSupport(State, AIndex, BIndex, Direction, CC);
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

internal evolve_result
EvolveSimplex(sdl_state* State, int32 AIndex, int32 BIndex)
{
  // TODO(Jovan): For now, under the assumption that it's 2 cubes
  evolve_result Result = NoIntersection;
  glm::vec3 Direction = glm::vec3(0.0f);
  cube* ShapeA = &State->Cubes[AIndex];
  cube* ShapeB = &State->Cubes[BIndex];
  simplex* Simplex = State->Simplex;
  uint32 NoVertices = Simplex->Count;
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
	glm::vec3 AB = Simplex->Vertices[1].P - Simplex->Vertices[0].P;
	// NOTE(Jovan): Form line from origin to A
	glm::vec3 A0 = -1.0f * Simplex->Vertices[0].P;

	// NOTE(Jovan): Direction perpendicular to AB in the direction
	// of the origin
	Direction = glm::cross(glm::cross(AB, A0), AB);
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
	// glm::vec3 DA = Simplex->Vertices[3] - Simplex->Vertices[0];
	// glm::vec3 DB = Simplex->Vertices[3] - Simplex->Vertices[1];
	// glm::vec3 DC = Simplex->Vertices[3] - Simplex->Vertices[2];

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

internal void
Barycentric(glm::vec3 P, glm::vec3 A, glm::vec3 B, glm::vec3 C,
	    real32* U, real32* V, real32* W)
{
  // TODO(Jovan): Use State?
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
ResolveCollision(sdl_state* State, int32 AIndex, int32 BIndex, collision_type Type)
{
  int32 CurrIter = 0;
  while(CurrIter++ <= MAX_EPA_ITERATIONS)
    {
      printf("Iter:%d\n", CurrIter);
      simplex* Simplex = State->Simplex;
      edge* Edge = State->Edge;
      triangle* Triangle = State->Triangle;

      Assert(Simplex->Count >= 4);
      // TODO(Jovan): Maybe use already loaded points for simplex
      // and just push triangles from their indices?
      vertex A = Simplex->Vertices[0];
      vertex B = Simplex->Vertices[1];
      vertex C = Simplex->Vertices[2];
      vertex D = Simplex->Vertices[3];
      PushTriangle(Triangle, A, B, C); // ABC
      PushTriangle(Triangle, A, C, D); // ACD
      PushTriangle(Triangle, A, D, B); // ADB
      PushTriangle(Triangle, B, D, C); // BDC

      // NOTE(Jovan): Generate new point in the direction
      // of the closest triangle
      uint32 ClosestIndex = 0;
      real32 CurrentDistance = FLT_MAX;
      for(uint32 TriangleIndex = 0;
	  TriangleIndex < Triangle->Count;
	  ++TriangleIndex)
	{
	  real32 Distance = glm::dot(Triangle->N[TriangleIndex], Triangle->A[TriangleIndex].P);
	  if(Distance < CurrentDistance)
	    {
	      CurrentDistance = Distance;
	      ClosestIndex = TriangleIndex;
	    }
	}
      printf("Dist: %f\n", CurrentDistance);
      
      glm::vec3 Direction = Triangle->N[ClosestIndex];
      vertex NewSupport = CalculateSupport(State, AIndex, BIndex, Direction, Type);
      // NOTE(Jovan): Removing triangle that can be "seen" by the new point
      for(uint32 TriangleIndex = 0;
	  TriangleIndex < Triangle->Count;)
	{
	  if(glm::dot(Triangle->N[TriangleIndex], NewSupport.P - Triangle->A[TriangleIndex].P) > 0)
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
	  glm::vec3 CollisionPoint = (BaryU * Triangle->A[ClosestIndex].SupA +
				      BaryV * Triangle->B[ClosestIndex].SupA +
				      BaryW * Triangle->C[ClosestIndex].SupA);
	  glm::vec3 CollisionNormal = -1.0f * Triangle->N[ClosestIndex];
	  real32 Depth = CurrentDistance;

	  printf("Point = %f, %f, %f | Depth:%f | Iteration: %d\nU: %f, V: %f, W: %f | Triangles: %d\n",
		 CollisionPoint.x,
		 CollisionPoint.y,
		 CollisionPoint.z,
		 Depth,
		 CurrIter,
		 BaryU,
		 BaryV,
		 BaryW,
		 Triangle->Count);
	  State->Spheres[0].Position = CollisionPoint;
	  //State->Cubes[0].Position += Depth * CollisionNormal;
	  break;
	}
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
      SimState->Camera.Position = glm::vec3(0.0f, 0.0f, 3.0f);
      glm::vec3 Up = glm::vec3(0.0f, 1.0f, 0.0f);
      SimState->Camera.Target = glm::vec3(0.0f, 0.0f, 0.0f);
      SimState->Camera.Direction = glm::normalize(SimState->Camera.Position - SimState->Camera.Target);
      SimState->Camera.Front = glm::vec3(0.0f, 0.0f, -1.0f);
      SimState->Camera.Right = glm::normalize(glm::cross(Up, SimState->Camera.Direction));
      SimState->Camera.Up = glm::cross(SimState->Camera.Direction, SimState->Camera.Right);
      
      // NOTE(Jovan): Cube init
      SimState->CubeCount = 2;
      SimState->Cubes[0].Position = glm::vec3(2.1,
					      0.5,
					      2.0);
      SimState->Cubes[0].Velocity = glm::vec3(0.0);
      SimState->Cubes[0].Forces = glm::vec3(0.0);
      SimState->Cubes[0].XAngle = 0.0f;
      SimState->Cubes[0].YAngle = 0.0f;
      SimState->Cubes[0].ZAngle = 0.0f;
      SimState->Cubes[0].Size = 1.0f;
      SimState->Cubes[0].Mass = 1.0f;
      UpdateVertices(SimState, 0);

      
      SimState->Cubes[1].Position = glm::vec3(2.0,
					      1.0,
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
  bool32 CollisionHappened = DetectCollisions(SimState, 0, 1, CC);
  if(CollisionHappened)
    {
      ResolveCollision(SimState, 0, 1, CC);
    }
  ClearVertices(SimState->Simplex);
  ClearTriangles(SimState->Triangle);
  ClearEdges(SimState->Edge);

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
  glm::vec3 LineColor = glm::vec3(0.0);
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
