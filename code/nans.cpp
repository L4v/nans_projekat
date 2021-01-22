#include "nans.h"
#include "utilities.cpp"
// NOTE(Jovan): Checks if a float is NaN or inf
static bool32
IsValid(real32 f)
{
    bool32 Result = 0;
    if (isnan(f) || isinf(f))
    {
        Result = 0;
    }
    else
    {
        Result = 1;
    }
    return Result;
}

static void
UpdateCamera(sdl_state *State, sdl_input *Input)
{
    State->Camera.Yaw += Input->MouseController.XRel * Input->MouseController.Sensitivity;
    State->Camera.Pitch += -Input->MouseController.YRel * Input->MouseController.Sensitivity;
    Input->MouseController.XRel = 0.0f;
    Input->MouseController.YRel = 0.0f;
    if (State->Camera.Pitch > 89.0f)
    {
        State->Camera.Pitch = 89.0f;
    }
    if (State->Camera.Pitch < -89.0f)
    {
        State->Camera.Pitch = -89.0f;
    }
    if (State->Camera.Position.y < 0.5)
    {
        State->Camera.Position.y = 0.5;
    }
}

static glm::vec3
Euler(glm::vec3 (*F)(glm::vec3, glm::vec3, real32), real32 dt,
      glm::vec3 Y0, glm::vec3 SummedForce, real32 Mass)
{
    glm::vec3 Result = {};
    glm::vec3 Tmp = F(Y0, SummedForce, Mass);
    Result = Y0 + dt * Tmp;
    return Result;
}

static glm::vec3
RK4(glm::vec3 (*F)(glm::vec3, glm::vec3, real32), real32 dt,
    glm::vec3 Y0, glm::vec3 SummedForce, real32 Mass)
{
    glm::vec3 Result = {};
    glm::vec3 K1 = dt * F(Y0, SummedForce, Mass);
    glm::vec3 K2 = dt * F(Y0 + (K1 / 2.0f), SummedForce, Mass);
    glm::vec3 K3 = dt * F(Y0 + (K2 / 2.0f), SummedForce, Mass);
    glm::vec3 K4 = dt * F(Y0 + K3, SummedForce, Mass);
    Result = Y0 + (1.0f / 6.0f) * (K1 + 2.0f * K2 + 2.0f * K3 + K4);
    return Result;
}

static glm::vec3
MovementFunction(glm::vec3 Velocity, glm::vec3 SummedForces, real32 Mass)
{
    // NOTE(Jovan): Not using "x" (position) so it is omitted
    glm::vec3 Result = {};
    Result = (real32)(1.0 / Mass) * (SummedForces + Mass * 9.81f * glm::vec3(0.0f, -1.0f, 0.0f) - GLOBAL_FRICTION * Velocity);
    return Result;
}

static glm::vec3
RotationFunction(glm::vec3 RotVelocity, glm::vec3 SummedTorque, real32 MomentOfInertia)
{
    glm::vec3 Result = {};
    Result = (real32)(1.0f / MomentOfInertia) * (SummedTorque - GLOBAL_FRICTION * RotVelocity);
    return Result;
}

static void
CubeAddForce(sdl_state *State, int32 CubeIndex, glm::vec3 Force)
{
    State->Cubes[CubeIndex].Forces += Force;
};

static void
CubeAddTorque(sdl_state *State, int32 CubeIndex, glm::vec3 Torque)
{
    State->Cubes[CubeIndex].Torque += Torque;
};

static void
CubeClearForces(sdl_state *State, int32 CubeIndex)
{
    State->Cubes[CubeIndex].Forces = glm::vec3(0.0);
    State->Cubes[CubeIndex].Torque = glm::vec3(0.0);
};

static void
SphereAddForce(sdl_state *State, int32 SphereIndex, glm::vec3 Force)
{
    State->Spheres[SphereIndex].Forces += Force;
}

static void
SphereClearForces(sdl_state *State, int32 SphereIndex)
{
    State->Spheres[SphereIndex].Forces = glm::vec3(0.0);
    State->Spheres[SphereIndex].Torque = glm::vec3(0.0);
}

static void
ShootSphere(sdl_state *State)
{
    SphereClearForces(State, 0);
    State->Spheres[0].V = glm::vec3(0.0f);
    State->Spheres[0].W = glm::vec3(0.0f);
    State->Spheres[0].Position = State->Camera.Position;
    SphereAddForce(State, 0, SHOOT_FORCE * State->Camera.Front);
}

static void
HandleInput(sdl_state *State, sdl_input *Input, real32 dt)
{

    if (Input->KeyboardController.MoveForward.EndedDown)
    {
        State->Camera.Position += State->Camera.Speed * State->Camera.Front;
    }
    if (Input->KeyboardController.MoveLeft.EndedDown)
    {
        State->Camera.Position -= glm::normalize(glm::cross(State->Camera.Front, State->Camera.Up)) * State->Camera.Speed;
    }
    if (Input->KeyboardController.MoveBack.EndedDown)
    {
        State->Camera.Position -= State->Camera.Speed * State->Camera.Front;
    }
    if (Input->KeyboardController.MoveRight.EndedDown)
    {
        State->Camera.Position += glm::normalize(glm::cross(State->Camera.Front, State->Camera.Up)) * State->Camera.Speed;
    }
    if (Input->KeyboardController.ShootAction.EndedDown)
    {
        ShootSphere(State); //CubeAddForce(State, 0, 20.0f * State->Camera.Front);
    }
    if (Input->KeyboardController.DebugLeft.EndedDown)
    {
        State->Cubes[0].Position += dt * glm::vec3(1.0, 0.0, 0.0);
    }
    if (Input->KeyboardController.DebugRight.EndedDown)
    {
        State->Cubes[0].Position += dt * glm::vec3(-1.0, 0.0, 0.0);
    }
    if (Input->KeyboardController.DebugUp.EndedDown)
    {
        State->Cubes[0].Position += dt * glm::vec3(0.0, 1.0, 0.0);
    }
    if (Input->KeyboardController.DebugDown.EndedDown)
    {
        State->Cubes[0].Position += dt * glm::vec3(0.0, -1.0, 0.0);
    }
    if (Input->KeyboardController.DebugForward.EndedDown)
    {
        State->Cubes[0].Position += dt * glm::vec3(0.0, 0.0, 1.0);
    }
    if (Input->KeyboardController.DebugBack.EndedDown)
    {
        State->Cubes[0].Position += dt * glm::vec3(0.0, 0.0, -1.0);
    }
    // TODO(Jovan): For debugging
    if (Input->KeyboardController.DebugReset.EndedDown)
    {
        State->Cubes[0].Position = glm::vec3(2.0, 3.5, 2.0);
        State->Cubes[0].V = glm::vec3(0.0f);
        State->Cubes[0].W = glm::vec3(0.0f);

        State->Cubes[1].Position = glm::vec3(2.0, 1.0, 2.0);
        State->Cubes[1].V = glm::vec3(0.0f);
        State->Cubes[1].W = glm::vec3(0.0f);

        State->Cubes[2].Position = glm::vec3(2.0, 4.5, 2.0);
        State->Cubes[2].V = glm::vec3(0.0f);
        State->Cubes[2].W = glm::vec3(0.0f);
    }
}

static void
RemoveVertex(simplex *Simplex, uint32 VertexIndex)
{
    Assert(VertexIndex < Simplex->Count);
    for (uint32 i = VertexIndex;
         i < Simplex->Count - 1;
         ++i)
    {
        Simplex->Vertices[i] = Simplex->Vertices[i + 1];
    }
    Simplex->Count--;
}
static void
PushVertex(simplex *Simplex, vertex V)
{
    Simplex->Vertices[Simplex->Count] = V;
    Simplex->Count++;
}

static inline void
ClearVertices(simplex *Simplex)
{
    Simplex->Count = 0;
}

static void
RemoveEdge(std::vector<edge> &Edge, uint32 EdgeIndex) //edge* Edge, uint32 EdgeIndex)
{
    // Assert(EdgeIndex < Edge->Count);
    // for(uint32 i = EdgeIndex;
    //     i < Edge->Count - 1;
    //     ++i)
    //   {
    //     Edge->A[i] = Edge->A[i + 1];
    //     Edge->B[i] = Edge->B[i + 1];
    //   }
    // Edge->Count--;
    Edge.erase(Edge.begin() + EdgeIndex);
}

static void
PushEdge(std::vector<edge> &Edge, vertex A, vertex B) //edge* Edge, vertex A, vertex B)
{
    // // NOTE(Jovan): Checks whether an edge of opposite
    // // winding already exists on the list.
    // // If it does, remove it and don't add the new one
    // for(uint32 EdgeIndex = 0;
    //     EdgeIndex < Edge->Count;
    //     ++EdgeIndex)
    //   {
    //     if(Edge->A[EdgeIndex] == B && Edge->B[EdgeIndex] == A)
    // 	{
    // 	  RemoveEdge(Edge, EdgeIndex);
    // 	  return;
    // 	}
    //   }
    // // NOTE(Jovan): The edge didn't exist, so it's being added
    // Edge->A[Edge->Count] = A;
    // Edge->B[Edge->Count] = B;
    // Edge->Count++;
    for (std::vector<edge>::iterator it = Edge.begin();
         it != Edge.end();
         ++it)
    {
        if (it->A == B && it->B == A)
        {
            Edge.erase(it);
            return;
        }
    }
    edge Tmp = {};
    Tmp.A = A;
    Tmp.B = B;
    Edge.push_back(Tmp);
}

static void
ClearEdges(edge *Edge)
{
    Edge->Count = 0;
}

static void
RemoveTriangle(std::vector<triangle> &Triangle, uint32 TriangleIndex) //triangle* Triangle, uint32 TriangleIndex)
{
    // Assert(TriangleIndex < Triangle->Count);
    // for(uint32 i = TriangleIndex;
    //     i < Triangle->Count - 1;
    //     ++i)
    //   {
    //     Triangle->A[i] = Triangle->A[i+1];
    //     Triangle->B[i] = Triangle->B[i+1];
    //     Triangle->C[i] = Triangle->C[i+1];
    //     Triangle->N[i] = Triangle->N[i+1];
    //   }
    // Triangle->Count--;
    Triangle.erase(Triangle.begin() + TriangleIndex);
}

// NOTE(Jovan): Push triangle using vertices
static void
PushTriangle(std::vector<triangle> &Triangle, vertex A, vertex B, vertex C) //triangle* Triangle, vertex A, vertex B, vertex C)
{
    // Triangle->A[Triangle->Count] = A;
    // Triangle->B[Triangle->Count] = B;
    // Triangle->C[Triangle->Count] = C;

    // // NOTE(Jovan): Calculate normal and make sure it's
    // // pointing away from the origin
    // Triangle->N[Triangle->Count] = glm::cross(B.P - A.P, C.P - A.P);
    // // NOTE(Jovan): Normalize vector
    // //real32 Len = sqrt(glm::dot(Triangle->N[Triangle->Count], Triangle->N[Triangle->Count]));
    // Triangle->N[Triangle->Count] = glm::normalize(Triangle->N[Triangle->Count]);

    // if(glm::dot(Triangle->A[Triangle->Count].P, Triangle->N[Triangle->Count]) < 0)
    //   {
    //     Triangle->N[Triangle->Count] *= -1.0f;
    //   }

    // Triangle->Count++;
    triangle Tmp = {};
    Tmp.A = A;
    Tmp.B = B;
    Tmp.C = C;
    Tmp.N = glm::normalize(glm::cross(B.P - A.P, C.P - A.P));
    if (glm::dot(Tmp.A.P, Tmp.N) < 0)
    {
        Tmp.N *= -1.0f;
    }
    Triangle.push_back(Tmp);
}

static inline void
ClearTriangles(triangle *Triangle)
{
    Triangle->Count = 0;
}

// NOTE(Jovan): Checks whether or not a contact pair already
// exists in the list
static bool32
PairExists(sdl_state *State, contact_pair *Pair)
{
    bool32 Result = 0;

    // for(uint32 PairIndex = 0;
    //     PairIndex < State->PairCount;
    //     ++PairIndex)
    //   {
    //     if(State->Pairs[PairIndex].IndexA == Pair->IndexA &&
    // 	 State->Pairs[PairIndex].IndexB == Pair->IndexB)
    // 	{
    // 	  Result = 1;
    // 	  break;
    // 	}
    //   }

    return Result;
}

// NOTE(Jovan): Pushes a contact pair into the list if it
// doesn't exist already
static void
PushPair(sdl_state *State, contact_pair Pair)
{
    // if(PairExists(State, &Pair))
    //   {
    //     return;
    //   }
    // State->Pairs[State->PairCount] = Pair;

    // State->PairCount++;
}

static void
RemovePair(sdl_state *State, int32 Index)
{
    // for(uint32 PairIndex = Index;
    //     PairIndex < State->PairCount - 1;
    //     ++PairIndex)
    //   {
    //     State->Pairs[PairIndex] = State->Pairs[PairIndex + 1];
    //   }

    // State->PairCount--;
}

static void
FloorUpdateVertices(sdl_state *State)
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

static void
UpdateVertices(sdl_state *State, int32 CubeIndex)
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

static glm::vec3
GetCubeSupport(sdl_state *State, int32 CubeIndex, glm::vec3 Direction)
{
    real32 MaxDistance = -FLT_MAX;
    glm::vec3 Result = glm::vec3(0.0f);

    // NOTE(Jovan) 8 vertices in a cube
    for (uint32 VertexIndex = 0;
         VertexIndex < 8;
         ++VertexIndex)
    {
        glm::vec3 CurrVertex = State->Cubes[CubeIndex].Vertices[VertexIndex];
        real32 Distance = glm::dot(CurrVertex, Direction);
        if (Distance > MaxDistance)
        {
            MaxDistance = Distance;
            Result = CurrVertex;
        }
    }

    return Result;
}

static glm::vec3
GetSphereSupport(sdl_state *State, int32 SphereIndex, glm::vec3 Direction)
{
    sphere *Sphere = &State->Spheres[SphereIndex];
    glm::vec3 Result = Sphere->Position + Sphere->Radius * glm::normalize(Direction);
    return Result;
}

static glm::vec3
GetFloorSupport(sdl_state *State, glm::vec3 Direction)
{
    real32 MaxDistance = -FLT_MAX;
    glm::vec3 Result = glm::vec3(0.0f);

    // NOTE(Jovan) 8 vertices in a cube
    for (uint32 VertexIndex = 0;
         VertexIndex < 8;
         ++VertexIndex)
    {
        glm::vec3 CurrVertex = State->Floor.Vertices[VertexIndex];
        real32 Distance = glm::dot(CurrVertex, Direction);
        if (Distance > MaxDistance)
        {
            MaxDistance = Distance;
            Result = CurrVertex;
        }
    }

    return Result;
}

static vertex
CalculateSupport(sdl_state *State, contact_pair *Pair, glm::vec3 Direction, std::vector<vertex> &Simplex)
{
    vertex Result = {};
    glm::vec3 SupportA, SupportB;
    switch (Pair->Type)
    {
    case CC:
    {
        glm::vec3 SupportA = GetCubeSupport(State, Pair->IndexA, Direction);
        glm::vec3 SupportB = GetCubeSupport(State, Pair->IndexB, -1.0f * Direction);
        Result.P = SupportA - SupportB;
        Result.SupA = SupportA;
        Result.SupB = SupportB;
    }
    break;
    case CS:
    {
        glm::vec3 SupportA = GetCubeSupport(State, Pair->IndexA, Direction);
        glm::vec3 SupportB = GetSphereSupport(State, Pair->IndexB, -1.0f * Direction);
        Result.P = SupportA - SupportB;
        Result.SupA = SupportA;
        Result.SupB = SupportB;
    }
    break;
    case CF:
    {
        glm::vec3 SupportA = GetCubeSupport(State, Pair->IndexA, Direction);
        glm::vec3 SupportB = GetFloorSupport(State, -1.0f * Direction);
        Result.P = SupportA - SupportB;
        Result.SupA = SupportA;
        Result.SupB = SupportB;
    }
    break;
    case SS:
    {
        glm::vec3 SupportA = GetSphereSupport(State, Pair->IndexA, Direction);
        glm::vec3 SupportB = GetSphereSupport(State, Pair->IndexB, -1.0f * Direction);
        Result.P = SupportA - SupportB;
        Result.SupA = SupportA;
        Result.SupB = SupportB;
    }
    break;
    case SF:
    {
        glm::vec3 SupportA = GetSphereSupport(State, Pair->IndexA, Direction);
        glm::vec3 SupportB = GetFloorSupport(State, -1.0f * Direction);
        Result.P = SupportA - SupportB;
        Result.SupA = SupportA;
        Result.SupB = SupportB;
    }
    break;
    }
    Simplex.push_back(Result);

    return Result;
}

static bool32
AddSupport(sdl_state *State, contact_pair *Pair, glm::vec3 Direction, std::vector<vertex> &Simplex)
{
    bool32 Result = 0;
    vertex NewVertex = {};
    NewVertex = CalculateSupport(State, Pair, Direction, Simplex);

    if (glm::dot(Direction, NewVertex.P) >= 0)
    {
        Result = 1;
    }
    else
    {
        Result = 0;
    }
    return Result;
}

static glm::vec3
ClosestPointOnLine(glm::vec3 A, glm::vec3 B, glm::vec3 Target,
                   real32 *U, real32 *V)
{
    glm::vec3 Result = glm::vec3(0.0f);
    glm::vec3 LineSegment = glm::normalize(B - A);
    real32 Length = glm::length(LineSegment);
    *V = glm::dot(-A, LineSegment) / Length;
    *U = glm::dot(B, LineSegment) / Length;
    if (*U <= 0.0f)
    {
        Result = B;
    }
    else if (*V <= 0.0f)
    {
        Result = A;
    }
    else
    {
        Result = (*U) * A + (*V) * B;
    }

    return Result;
}

static glm::vec3
TripleCross(glm::vec3 A, glm::vec3 B, glm::vec3 C)
{
    // NOTE(Jovan): (A x B) x C = B * (C . A) - A * (C . B)
    return (B * glm::dot(C, A)) - (A * glm::dot(C, B));
}

static evolve_result
EvolveSimplex(sdl_state *State, contact_pair *Pair, glm::vec3 PositionA, glm::vec3 PositionB, std::vector<vertex> &Simplex)
{
    evolve_result Result = StillEvolving;
    glm::vec3 Direction = glm::normalize(PositionB - PositionA); //glm::normalize(glm::vec3(1.0f));
    uint32 NoVertices = Simplex.size();
    switch (NoVertices)
    {
    case 0:
    {
        // NOTE(Jovan): Best presumption of dir is the relative direction
        // of the 2 shapes

        //Direction = PositionB - PositionA;
    }
    break;
    case 1:
    {
        // NOTE(Jovan): Flip the direction to point to the origin
        Direction *= -1;
    }
    break;
    case 2:
    {
        // // NOTE(Jovan): Form line from first 2 vertices
        glm::vec3 AB = Simplex[1].P - Simplex[0].P;
        // // NOTE(Jovan): Form line from origin to A
        glm::vec3 A0 = -1.0f * Simplex[0].P;

        real32 U = 0.0f;
        real32 V = 0.0f;
        glm::vec3 Origin = glm::vec3(0.0);
        glm::vec3 ClosestPoint = ClosestPointOnLine(Simplex[0].P,
                                                    Simplex[1].P,
                                                    Origin, &U, &V);
        if (V <= 0.0f)
        {
            Simplex.erase(Simplex.begin() + 1);
            Direction = -ClosestPoint;
        }
        else if (U <= 0.0f)
        {
            Simplex.erase(Simplex.begin());
            Direction = -ClosestPoint;
        }
        else
        {
            Direction = -ClosestPoint;
        }
    }
    break;
    case 3:
    {
#if 0
	glm::vec3 AC = Simplex[2].P - Simplex[0].P;
	glm::vec3 AB = Simplex[1].P - Simplex[0].P;
	Direction = glm::cross(AC, AB);

	// NOTE(Jovan): Ensure that Direction points to the origin
	glm::vec3 A0 = -1.0f * Simplex[0].P;
	if(glm::dot(Direction, A0) < 0)
	  {
	    Direction *= -1.0f;
	  }
#endif
        glm::vec3 NewToOrigin = -Simplex[0].P;         //A0
        glm::vec3 Edge1 = Simplex[1].P - Simplex[0].P; //AB
        glm::vec3 Edge2 = Simplex[2].P - Simplex[0].P; //AC

        glm::vec3 TriangleNormal = glm::cross(Edge1, Edge2);
        glm::vec3 Edge1Normal = glm::cross(Edge1, TriangleNormal);
        glm::vec3 Edge2Normal = glm::cross(TriangleNormal, Edge2);

        // NOTE(Jovan): If closer along second edge normal
        if (glm::dot(Edge2Normal, NewToOrigin) > 0.0f)
        {
            if (glm::dot(Edge2, NewToOrigin) > 0.0f)
            {
                Direction = TripleCross(Edge2, NewToOrigin, Edge2);
                Simplex.erase(Simplex.begin() + 1);
                break;
            }
            else
            {
                // NOTE(jovan): Star case
                if (glm::dot(Edge1, NewToOrigin) > 0.0f)
                {
                    Direction = TripleCross(Edge1, NewToOrigin, Edge1);
                    Simplex.erase(Simplex.begin() + 2);
                    break;
                }
                else
                {
                    Direction = NewToOrigin;
                    Simplex.erase(Simplex.begin() + 2);
                    Simplex.erase(Simplex.begin() + 1);
                    break;
                }
            }
        }
        else
        {
            // NOTE(Jovan): Star case
            if (glm::dot(Edge1Normal, NewToOrigin) > 0.0f)
            {
                if (glm::dot(Edge1, NewToOrigin) > 0.0f)
                {
                    Direction = TripleCross(Edge1, NewToOrigin, Edge1);
                    Simplex.erase(Simplex.begin() + 2);
                    break;
                }
                else
                {
                    Direction = NewToOrigin;
                    Simplex.erase(Simplex.begin() + 2);
                    Simplex.erase(Simplex.begin() + 1);
                    break;
                }
            }
            else
            {
                if (glm::dot(TriangleNormal, NewToOrigin) > 0.0f)
                {
                    Direction = TriangleNormal;
                    break;
                }
                else
                {
                    Direction = -TriangleNormal;
                    vertex Tmp = Simplex[1];
                    Simplex[1] = Simplex[2];
                    Simplex[2] = Tmp;
                    break;
                }
            }
        }
    }
    break;
    case 4:
    {
        // NOTE(Jovan): 3 edges of interest
        glm::vec3 DA = Simplex[0].P - Simplex[3].P;
        glm::vec3 DB = Simplex[1].P - Simplex[3].P;
        glm::vec3 DC = Simplex[2].P - Simplex[3].P;

        // NOTE(Jovan): Dir to the origin
        glm::vec3 D0 = -1.0f * Simplex[3].P;

        // NOTE(Jovan): Check triangles ABD, BCD, CAD
        glm::vec3 ABDNorm = glm::cross(DA, DB);
        glm::vec3 BCDNorm = glm::cross(DB, DC);
        glm::vec3 CADNorm = glm::cross(DC, DA);

        if (glm::dot(ABDNorm, D0) > 0.0f)
        {
            // NOTE(Jovan): Origin outside of ABD -> remove C
            Simplex.erase(Simplex.begin() + 2);
            Direction = ABDNorm;
        }
        else if (glm::dot(BCDNorm, D0) > 0.0f)
        {
            // NOTE(Jovan): Origin outside of BCD -> remove A
            Simplex.erase(Simplex.begin());
            Direction = BCDNorm;
        }
        else if (glm::dot(CADNorm, D0) > 0.0f)
        {
            // NOTE(Jovan): Origin outside of CAD -> remove B
            Simplex.erase(Simplex.begin() + 1);
            Direction = CADNorm;
        }
        else
        {
            // NOTE(Jovan): Origin is inside of the tetrahedron
            Result = FoundIntersection;
            return Result;
        }
    }
    break;
    default:
    {
        printf("ERROR::There can't be more than 4 vertices\n");
    }
    }
    if (glm::length(Direction) <= 0.0001f)
    {
        Result = NoIntersection;
        return Result;
    }
    if (AddSupport(State, Pair, Direction, Simplex))
    {
        Result = StillEvolving;
    }
    else
    {
        Result = NoIntersection;
    }
    return Result;
}

static void
Barycentric(glm::vec3 P, glm::vec3 A, glm::vec3 B, glm::vec3 C,
            real32 *U, real32 *V, real32 *W)
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

static bool32
ResolveCollision(sdl_state *State, contact_pair *Pair, std::vector<vertex> &Simplex)
{
    int32 CurrIter = 0;
    std::vector<edge> Edge;
    std::vector<triangle> Triangle;
    glm::vec3 Result;
    // NOTE(Jovan): Take over points from GJK and construct a tetrahedron
    vertex A = Simplex[0];
    vertex B = Simplex[1];
    vertex C = Simplex[2];
    vertex D = Simplex[3];
    PushTriangle(Triangle, A, B, C); // ABC
    PushTriangle(Triangle, A, C, D); // ACD
    PushTriangle(Triangle, A, D, B); // ADB
    PushTriangle(Triangle, B, D, C); // BDC

    while (CurrIter++ <= MAX_EPA_ITERATIONS)
    {
        // NOTE(Jovan): Find the closest triangle
        real32 CurrentDistance = glm::abs(glm::dot(Triangle[0].N, Triangle[0].A.P));
        std::vector<triangle>::iterator ClosestIt = Triangle.begin();
        for (std::vector<triangle>::iterator it = Triangle.begin();
             it != Triangle.end();
             ++it)
        {
            glm::vec3 AB = it->B.P - it->A.P;
            glm::vec3 AC = it->C.P - it->A.P;
            glm::vec3 Normal = glm::normalize(glm::cross(AB, AC));
            real32 Distance = glm::abs(glm::dot(Normal, it->A.P));
            if (Distance < CurrentDistance)
            {
                CurrentDistance = Distance;
                ClosestIt = it;
            }
        }

        glm::vec3 Direction = ClosestIt->N; //Triangle->N[ClosestIndex];
        vertex NewSupport = CalculateSupport(State, Pair, Direction, Simplex);

        // NOTE(Jovan): Extrapolate contact data
        if (glm::dot(ClosestIt->N, NewSupport.P) - CurrentDistance < MAX_EPA_ERROR)
        {
            real32 BaryU, BaryV, BaryW;
            Barycentric(ClosestIt->N * CurrentDistance,
                        ClosestIt->A.P,
                        ClosestIt->B.P,
                        ClosestIt->C.P,
                        &BaryU, &BaryV, &BaryW);
            if (fabs(BaryU) > 1.0f || fabs(BaryV) > 1.0f || fabs(BaryW) > 1.0f)
            {
                // NOTE(Jovan): No collision
                return 0;
            }

            if (!IsValid(BaryU) || !IsValid(BaryV) || !IsValid(BaryW))
            {
                // NOTE(Jovan): No collision
                return 0;
            }
            // NOTE(Jovan): Contact info for body A
            glm::vec3 CollisionPoint = ((BaryU * ClosestIt->A.SupA) +
                                        (BaryV * ClosestIt->B.SupA) +
                                        (BaryW * ClosestIt->C.SupA));
            glm::vec3 CollisionNormal = -1.0f * ClosestIt->N;
            real32 Depth = CurrentDistance;

            // NOTE(Jovan): Return the collision Point
            Pair->PointA = CollisionPoint;
            Pair->N = CollisionNormal;

            // NOTE(Jovan): Return the closest index for indexing the collision normal
            // from the triangle list

            // NOTE(Jovan): Contact info for body B
            CollisionPoint = ((BaryU * ClosestIt->A.SupB) +
                              (BaryV * ClosestIt->B.SupB) +
                              (BaryW * ClosestIt->C.SupB));
            Pair->PointB = CollisionPoint;
            return 1;
        }
        // NOTE(Jovan): Removing triangle that can be "seen" by the new point
        for (std::vector<triangle>::iterator it = Triangle.begin();
             it != Triangle.end();)
        {
            glm::vec3 Temp = NewSupport.P - it->A.P;
            glm::vec3 AB = it->B.P - it->A.P;
            glm::vec3 AC = it->C.P - it->A.P;
            glm::vec3 Normal = glm::normalize(glm::cross(AB, AC));
            if (glm::dot(Normal, it->A.P) < 0)
            {
                Normal *= -1.0f;
            }
            if (glm::dot(Normal, Temp) > 0)
            {
                // NOTE(Jovan): "Disolve" the removed triangle into it's edges
                // and push them onto the edge list
                PushEdge(Edge, it->A, it->B);
                PushEdge(Edge, it->B, it->C);
                PushEdge(Edge, it->C, it->A);
                Triangle.erase(it);
                continue;
            }
            ++it;
        }

        // NOTE(Jovan): Construct new triangles
        for (std::vector<edge>::iterator it = Edge.begin();
             it != Edge.end();
             ++it)
        {
            PushTriangle(Triangle, NewSupport, it->A, it->B);
        }
        // NOTE(Jovan): Clear the edge list
        Edge.clear();
    }
    return 0;
}

static bool32
CheckCollision(sdl_state *State, contact_pair *Pair, std::vector<vertex> &Simplex)
{
    bool32 Result = 0;
    evolve_result EvolutionResult = StillEvolving;
    int32 IndexA = Pair->IndexA;
    int32 IndexB = Pair->IndexB;
    uint32 GJKIteration = 0;
    Simplex.clear();
    while ((EvolutionResult == StillEvolving) && (GJKIteration++ <= MAX_GJK_ITERATIONS))
    {
        switch (Pair->Type)
        {
        case CC:
        {
            EvolutionResult = EvolveSimplex(State, Pair,
                                            State->Cubes[IndexA].Position,
                                            State->Cubes[IndexB].Position, Simplex);
        }
        break;
        case CS:
        {
            EvolutionResult = EvolveSimplex(State, Pair,
                                            State->Cubes[IndexA].Position,
                                            State->Spheres[IndexB].Position, Simplex);
        }
        break;
        case CF:
        {
            EvolutionResult = EvolveSimplex(State, Pair,
                                            State->Cubes[Pair->IndexA].Position,
                                            State->Floor.Position, Simplex);
        }
        break;
        case SS:
        {
            EvolutionResult = EvolveSimplex(State, Pair,
                                            State->Spheres[IndexA].Position,
                                            State->Spheres[IndexB].Position, Simplex);
        }
        break;
        case SF:
        {
            EvolutionResult = EvolveSimplex(State, Pair,
                                            State->Spheres[IndexA].Position,
                                            State->Floor.Position, Simplex);
        }
        break;
        }
    }

    if (EvolutionResult == FoundIntersection)
    {
        Result = ResolveCollision(State, Pair, Simplex); //1;
    }
    if (EvolutionResult == NoIntersection)
    {
        Result = 0;
    }
    return Result;
}

static void
DrawTriangles(sdl_state *State, sdl_render *Render)
{
    // TODO(Jovan): Implement
}

static void
IntegrateForces(sdl_state *State, real32 dt)
{
    glm::vec3 Y0 = {};
    glm::vec3 Y = {};
    for (uint32 CubeIndex = 0;
         CubeIndex < State->CubeCount;
         ++CubeIndex)
    {
        // NOTE(Jovan): Linear

        Y0 = State->Cubes[CubeIndex].V;
        Y = RK4(MovementFunction, dt, Y0, State->Cubes[CubeIndex].Forces, State->Cubes[CubeIndex].Mass);
        State->Cubes[CubeIndex].V = Y;

        // NOTE(Jovan): Rotational
        Y0 = State->Cubes[CubeIndex].W;
        Y = RK4(RotationFunction, dt, Y0, State->Cubes[CubeIndex].Torque, State->Cubes[CubeIndex].MOI);
        State->Cubes[CubeIndex].W = Y;

        // NOTE(Jovan): Update states and clear the forces/torques
        UpdateVertices(State, CubeIndex);
        CubeClearForces(State, CubeIndex);
    }

    for (uint32 SphereIndex = 0;
         SphereIndex < State->SphereCount;
         ++SphereIndex)
    {
        // NOTE(Jovan): Linear
        Y0 = State->Spheres[SphereIndex].V;
        Y = RK4(MovementFunction, dt, Y0, State->Spheres[SphereIndex].Forces, State->Spheres[SphereIndex].Mass);
        State->Spheres[SphereIndex].V = Y;

        // // NOTE(Jovan): Rotational
        Y0 = State->Spheres[SphereIndex].W;
        Y = RK4(RotationFunction, dt, Y0, State->Spheres[SphereIndex].Torque, State->Spheres[SphereIndex].MOI);
        State->Spheres[SphereIndex].W = Y;

        // NOTE(Jovan): Update states and clear the forces/torques
        SphereClearForces(State, SphereIndex);
    }

    FloorUpdateVertices(State);
}

static void
Constraint(sdl_state *State, contact_pair *Pair, real32 dt)
{
    glm::vec3 N = glm::normalize(Pair->N);
    glm::vec3 T1, T2;

    // NOTE(Jovan): Calculating J(M^-1)(J^T)
    // ------------------------------------
    glm::vec3 R1, R2;
    glm::vec3 PosA, PosB;
    glm::vec3 V1, W1, V2, W2;
    // NOTE(Jovan): I^-1
    real32 InvI1, InvI2;
    // NOTE(Jovan): M^-1
    real32 InvM1, InvM2;
    switch (Pair->Type)
    {
    case CC:
    {
        PosA = State->Cubes[Pair->IndexA].Position;
        PosB = State->Cubes[Pair->IndexB].Position;
        InvI1 = 1.0f / State->Cubes[Pair->IndexA].MOI;
        InvI2 = 1.0f / State->Cubes[Pair->IndexB].MOI;
        InvM1 = 1.0f / State->Cubes[Pair->IndexA].Mass;
        InvM2 = 1.0f / State->Cubes[Pair->IndexB].Mass;
        V1 = State->Cubes[Pair->IndexA].V;
        V2 = State->Cubes[Pair->IndexB].V;
        W1 = State->Cubes[Pair->IndexA].W;
        W2 = State->Cubes[Pair->IndexB].W;
    }
    break;
    case CS:
    {
        PosA = State->Cubes[Pair->IndexA].Position;
        PosB = State->Spheres[Pair->IndexB].Position;
        InvI1 = 1.0f / State->Cubes[Pair->IndexA].MOI;
        InvI2 = 1.0f / State->Spheres[Pair->IndexB].MOI;
        InvM1 = 1.0f / State->Cubes[Pair->IndexA].Mass;
        InvM2 = 1.0f / State->Spheres[Pair->IndexB].Mass;
        V1 = State->Cubes[Pair->IndexA].V;
        V2 = State->Spheres[Pair->IndexB].V;
        W1 = State->Cubes[Pair->IndexA].W;
        W2 = State->Spheres[Pair->IndexB].W;
    }
    break;
    case CF:
    {
        PosA = State->Cubes[Pair->IndexA].Position;
        PosB = State->Floor.Position;
        InvI1 = 1.0f / State->Cubes[Pair->IndexA].MOI;
        InvI2 = 1.0f / State->Floor.MOI;
        InvM1 = 1.0f / State->Cubes[Pair->IndexA].Mass;
        InvM2 = 1.0f / State->Floor.Mass;
        V1 = State->Cubes[Pair->IndexA].V;
        V2 = State->Floor.V;
        W1 = State->Cubes[Pair->IndexA].W;
        W2 = State->Floor.W;
    }
    break;
    case SS:
    {
        PosA = State->Spheres[Pair->IndexA].Position;
        PosB = State->Spheres[Pair->IndexB].Position;
        InvI1 = 1.0f / State->Spheres[Pair->IndexA].MOI;
        InvI2 = 1.0f / State->Spheres[Pair->IndexB].MOI;
        InvM1 = 1.0f / State->Spheres[Pair->IndexA].Mass;
        InvM2 = 1.0f / State->Spheres[Pair->IndexB].Mass;
        V1 = State->Spheres[Pair->IndexA].V;
        V2 = State->Spheres[Pair->IndexB].V;
        W1 = State->Spheres[Pair->IndexA].W;
        W2 = State->Spheres[Pair->IndexB].W;
    }
    break;
    case SF:
    {
        PosA = State->Spheres[Pair->IndexA].Position;
        PosB = State->Floor.Position;
        InvI1 = 1.0f / State->Spheres[Pair->IndexA].MOI;
        InvI2 = 1.0f / State->Floor.MOI;
        InvM1 = 1.0f / State->Spheres[Pair->IndexA].Mass;
        InvM2 = 1.0f / State->Floor.Mass;
        V1 = State->Spheres[Pair->IndexA].V;
        V2 = State->Floor.V;
        W1 = State->Spheres[Pair->IndexA].W;
        W2 = State->Floor.W;
    }
    break;
    default:
    {
        printf("ERROR::CONSTRAINT::Unknown collision type\n");
        return;
    }
    break;
    }

    if (N == glm::vec3(0.0))
    {
        N = glm::normalize(PosB - PosA); //return;
        printf("It happened!\n");
    }
    R1 = Pair->PointA - PosA;
    R2 = Pair->PointB - PosB;

    // NOTE(Jovan): Calculating tangents (allen chou)
    // sqrt(1/3) = 0.57735 (http://box2d.org/2014/02/computing-a-basis/)
    if (N.x >= 0.57735f)
    {
        T1 = glm::vec3(N.y, -N.x, 0.0f);
    }
    else
    {
        T1 = glm::vec3(0.0f, N.z, -N.y);
    }
    glm::normalize(T1);
    T2 = glm::cross(N, T1);
    Pair->T1 = T1;
    Pair->T2 = T2;

    // NOTE(Jovan): Penetration depth
    real32 Depth = glm::dot(((PosA + R1) - (PosB + R2)), N);
    // NOTE(Jovan): R x N
    glm::vec3 RN1 = glm::cross(R1, N);
    glm::vec3 RN2 = glm::cross(R2, N);
    // NOTE(Jovan): JMJ normal
    real32 JMJn = (InvM1) + (InvM2);
    JMJn += InvI1 * (glm::dot(RN1, RN1)) - InvI2 * (glm::dot(-RN2, -RN2));
    // NOTE(Jovan): 1 / (JMJn) (Effective mass)
    JMJn = 1.0f / JMJn;
    // NOTE(Jovan): Derivative of V
    glm::vec3 dVn = V1 + glm::cross(W1, N) - V2 - glm::cross(W2, N);
    real32 JdVn = glm::dot(dVn, N);
    // NOTE(Jovan): Baumgarte
    real32 Beta = 0.3f;
    // NOTE(Jovan): Restitution
    real32 Cr = 0.1;
    // NOTE(Jovan): Bias
    real32 B = -Beta / dt * Depth + Cr * JdVn;

    // NOTE(Jovan): Friction coefficient
    real32 Cf = 0.1f;
    // NOTE(Jovan): Friction tangents
    glm::vec3 R1T1 = glm::cross(R1, T1);
    glm::vec3 R2T1 = glm::cross(R2, T1);
    glm::vec3 R1T2 = glm::cross(R1, T2);
    glm::vec3 R2T2 = glm::cross(R2, T2);
    real32 JMJt1 = (InvM1) + (InvM2);
    JMJt1 += InvI1 * (glm::dot(R1T1, R1T1)) - InvI2 * (glm::dot(-R2T1, -R2T1));
    JMJt1 = 1.0f / JMJt1;
    real32 JMJt2 = (InvM1) + (InvM2);
    JMJt2 += InvI1 * (glm::dot(R1T2, R1T2)) - InvI2 * (glm::dot(-R2T2, -R2T2));
    JMJt2 = 1.0f / JMJt2;
    glm::vec3 dVt1 = V1 + glm::cross(W1, T1) - V2 - glm::cross(W2, T1);
    real32 JdVt1 = glm::dot(dVt1, T1);
    glm::vec3 dVt2 = V1 + glm::cross(W1, T2) - V2 - glm::cross(W2, T2);
    real32 JdVt2 = glm::dot(dVt2, T2);

    int32 Iter = 70;
    while (Iter--)
    {
        // NOTE(Jovan): Calculating normal lambda
        // --------------------------------------
        real32 LambdaN = (-JdVn + B) * JMJn;
        // NOTE(Jovan): Accumulating and clamping impulse
        real32 OldAccumI = Pair->DLNormalSum;
        Pair->DLNormalSum += LambdaN;
        if (Pair->DLNormalSum < 0)
        {
            Pair->DLNormalSum = 0.0f;
        }
        Pair->DLNormal = Pair->DLNormalSum - OldAccumI;
        // NOTE(Jovan): End of normal
        // ----------------------------

        // NOTE(Jovan): Calculating tangent lambda
        // ---------------------------------------
        real32 LambdaT1 = (-JdVt1) * JMJt1;
        real32 OldAccumT1 = Pair->DLTangent1Sum;
        Pair->DLTangent1Sum += LambdaT1;
        // NOTE(Jovan): For clamping friction lambda to [-Cf*Lambda, Cf*Lambda]
        real32 MaxLT1 = sqrt(2) * Cf * Pair->DLNormalSum;
        if (Pair->DLTangent1Sum < -MaxLT1)
        {
            Pair->DLTangent1Sum = -MaxLT1;
        }
        if (Pair->DLTangent1Sum > MaxLT1)
        {
            Pair->DLTangent1Sum = MaxLT1;
        }
        Pair->DLTangent1 = Pair->DLTangent1Sum - OldAccumT1;

        real32 LambdaT2 = (-JdVt2) * JMJt2;
        real32 OldAccumT2 = Pair->DLTangent2Sum;
        Pair->DLTangent2Sum += LambdaT2;
        // NOTE(Jovan): For clamping friction lambda to [-Cf*Lambda, Cf*Lambda]
        real32 MaxLT2 = sqrt(2) * Cf * Pair->DLNormalSum;
        if (Pair->DLTangent2Sum < -MaxLT2)
        {
            Pair->DLTangent2Sum = -MaxLT2;
        }
        if (Pair->DLTangent2Sum > MaxLT2)
        {
            Pair->DLTangent2Sum = MaxLT2;
        }
        Pair->DLTangent2 = Pair->DLTangent2Sum - OldAccumT2;

        // NOTE(Jovan): End of tangents
        // ----------------------------
    }
    // NOTE(Jovan): Calculating linear impulse
    glm::vec3 LinearImpulse = N * Pair->DLNormal;
    glm::vec3 LinearImpulseT1 = T1 * Pair->DLTangent1;
    glm::vec3 LinearImpulseT2 = T2 * Pair->DLTangent2;
    // NOTE(Jovan): Calculating angular impulses
    glm::vec3 AngularI1 = RN1 * Pair->DLNormal;
    glm::vec3 AngularI2 = RN2 * Pair->DLNormal;
    glm::vec3 AngularI1T1 = R1T1 * Pair->DLTangent1;
    glm::vec3 AngularI2T1 = R2T1 * Pair->DLTangent1;
    glm::vec3 AngularI1T2 = R1T2 * Pair->DLTangent2;
    glm::vec3 AngularI2T2 = R2T2 * Pair->DLTangent2;
    // NOTE(Jovan): Applying linear and angular impulses
    switch (Pair->Type)
    {
    case CC:
    {
        State->Cubes[Pair->IndexA].V += InvM1 * LinearImpulse;
        State->Cubes[Pair->IndexB].V -= InvM2 * LinearImpulse;
        State->Cubes[Pair->IndexA].W += InvI1 * AngularI1;
        State->Cubes[Pair->IndexB].W -= InvI2 * AngularI2;

        State->Cubes[Pair->IndexA].V += InvM1 * LinearImpulseT1;
        State->Cubes[Pair->IndexB].V -= InvM2 * LinearImpulseT1;
        State->Cubes[Pair->IndexA].W += InvI1 * AngularI1T1;
        State->Cubes[Pair->IndexB].W -= InvI2 * AngularI2T1;

        State->Cubes[Pair->IndexA].V += InvM1 * LinearImpulseT2;
        State->Cubes[Pair->IndexB].V -= InvM2 * LinearImpulseT2;
        State->Cubes[Pair->IndexA].W += InvI1 * AngularI1T2;
        State->Cubes[Pair->IndexB].W -= InvI2 * AngularI2T2;
    }
    break;
    case CS:
    {
        State->Cubes[Pair->IndexA].V += InvM1 * LinearImpulse;
        State->Spheres[Pair->IndexB].V -= InvM2 * LinearImpulse;
        State->Cubes[Pair->IndexA].W += InvI1 * AngularI1;
        State->Spheres[Pair->IndexB].W -= InvI2 * AngularI2;

        State->Cubes[Pair->IndexA].V += InvM1 * LinearImpulseT1;
        State->Spheres[Pair->IndexB].V -= InvM2 * LinearImpulseT1;
        State->Cubes[Pair->IndexA].W += InvI1 * AngularI1T1;
        State->Spheres[Pair->IndexB].W -= InvI2 * AngularI2T1;

        State->Cubes[Pair->IndexA].V += InvM1 * LinearImpulseT2;
        State->Spheres[Pair->IndexB].V -= InvM2 * LinearImpulseT2;
        State->Cubes[Pair->IndexA].W += InvI1 * AngularI1T2;
        State->Spheres[Pair->IndexB].W -= InvI2 * AngularI2T2;
    }
    break;
    case CF:
    {
        // NOTE(Jovan): Not updating floor velocities
        State->Cubes[Pair->IndexA].V += InvM1 * LinearImpulse;
        State->Cubes[Pair->IndexA].W += InvI1 * AngularI1;

        State->Cubes[Pair->IndexA].V += InvM1 * LinearImpulseT1;
        State->Cubes[Pair->IndexA].W += InvI1 * AngularI1T1;

        State->Cubes[Pair->IndexA].V += InvM1 * LinearImpulseT2;
        State->Cubes[Pair->IndexA].W += InvI1 * AngularI1T2;
    }
    break;
    case SS:
    {
        State->Spheres[Pair->IndexA].V += InvM1 * LinearImpulse;
        State->Spheres[Pair->IndexB].V -= InvM2 * LinearImpulse;
        State->Spheres[Pair->IndexA].W += InvI1 * AngularI1;
        State->Spheres[Pair->IndexB].W -= InvI2 * AngularI2;

        State->Spheres[Pair->IndexA].V += InvM1 * LinearImpulseT1;
        State->Spheres[Pair->IndexB].V -= InvM2 * LinearImpulseT1;
        State->Spheres[Pair->IndexA].W += InvI1 * AngularI1T1;
        State->Spheres[Pair->IndexB].W -= InvI2 * AngularI2T1;

        State->Spheres[Pair->IndexA].V += InvM1 * LinearImpulseT2;
        State->Spheres[Pair->IndexB].V -= InvM2 * LinearImpulseT2;
        State->Spheres[Pair->IndexA].W += InvI1 * AngularI1T2;
        State->Spheres[Pair->IndexB].W -= InvI2 * AngularI2T2;
    }
    break;
    case SF:
    {
        // NOTE(Jovan): Not updating floor velocities
        State->Spheres[Pair->IndexA].V += InvM1 * LinearImpulse;
        State->Spheres[Pair->IndexA].W += InvI1 * AngularI1;

        State->Spheres[Pair->IndexA].V += InvM1 * LinearImpulseT1;
        State->Spheres[Pair->IndexA].W += InvI1 * AngularI1T1;

        State->Spheres[Pair->IndexA].V += InvM1 * LinearImpulseT2;
        State->Spheres[Pair->IndexA].W += InvI1 * AngularI1T2;
    }
    break;
    default:
    {
        printf("ERROR::CONSTRAINT::Unknown collision type\n");
        return;
    }
    break;
    }
}

static void
IntegrateVelocities(sdl_state *State, real32 dt)
{
    for (uint32 CubeIndex = 0;
         CubeIndex < State->CubeCount;
         ++CubeIndex)
    {
        State->Cubes[CubeIndex].Position += dt * State->Cubes[CubeIndex].V;
        State->Cubes[CubeIndex].Angles += dt * State->Cubes[CubeIndex].W;
        UpdateVertices(State, CubeIndex);
    }
    for (uint32 SphereIndex = 0;
         SphereIndex < State->SphereCount;
         ++SphereIndex)
    {
        State->Spheres[SphereIndex].Position += dt * State->Spheres[SphereIndex].V;
        State->Spheres[SphereIndex].Angles += dt * State->Spheres[SphereIndex].W;
    }
}

static void
DetectCollisions(sdl_state *State, sdl_input *Input, sdl_render *Render, real32 dt,
                 std::vector<contact_pair> &Pairs)
{
#if CC_COLLISIONS
    // NOTE(Jovan): CC
    for (uint32 Cube1 = 0;
         Cube1 < State->CubeCount;
         ++Cube1)
    {
        for (uint32 Cube2 = Cube1 + 1;
             Cube2 < State->CubeCount;
             ++Cube2)
        {
            contact_pair Pair = {};
            Pair.Type = CC;
            Pair.IndexA = Cube1;
            Pair.IndexB = Cube2;
            bool32 CollisionHappened = 0;
            std::vector<vertex> Simplex;
            CollisionHappened = CheckCollision(State, &Pair, Simplex);
            if (CollisionHappened)
            {
                bool32 Exists = 0;
                for (std::vector<contact_pair>::iterator it = Pairs.begin();
                     it != Pairs.end();
                     ++it)
                {
                    if ((it->IndexA == Pair.IndexA && it->IndexB == Pair.IndexB && it->Type == Pair.Type) ||
                        (it->IndexA == Pair.IndexB && it->IndexB == Pair.IndexA && it->Type == Pair.Type))
                    {
                        printf("Cache hit!\n");
                        Exists = 1;
                    }
                }
                if (Exists)
                {
                    // TODO(Jovan): Warm starting?
                }
                else
                {
                    Pairs.push_back(Pair);
                }
            }
        }
    }
#endif
#if FLOOR_COLLISIONS
    // NOTE(Jovan): With floor
    for (uint32 CubeIndex = 0;
         CubeIndex < State->CubeCount;
         ++CubeIndex)
    {
        contact_pair Pair = {};
        Pair.Type = CF;
        Pair.IndexA = CubeIndex;
        Pair.IndexB = CubeIndex;
        bool32 CollisionHappened = 0;
        std::vector<vertex> Simplex;
        CollisionHappened = CheckCollision(State, &Pair, Simplex);
        if (CollisionHappened)
        {
            bool32 Exists = 0;
            for (std::vector<contact_pair>::iterator it = Pairs.begin();
                 it != Pairs.end();
                 ++it)
            {
                if ((it->IndexA == Pair.IndexA && it->IndexB == Pair.IndexB) ||
                    (it->IndexA == Pair.IndexB && it->IndexB == Pair.IndexA))
                {
                    Exists = 1;
                    // NOTE(Jovan): Warm starting?
                }
            }
            if (Exists == 0)
            {
                Pairs.push_back(Pair);
            }
        }
    }
    for (uint32 SphereIndex = 0;
         SphereIndex < State->SphereCount;
         ++SphereIndex)
    {
        contact_pair Pair = {};
        std::vector<vertex> Simplex;
        Pair.Type = SF;
        Pair.IndexA = SphereIndex;
        Pair.IndexB = SphereIndex;
        bool32 CollisionHappened = CheckCollision(State, &Pair, Simplex);
        if (CollisionHappened)
        {
            bool32 Exists = 0;
            for (std::vector<contact_pair>::iterator it = Pairs.begin();
                 it != Pairs.end();
                 ++it)
            {
                if ((it->IndexA == Pair.IndexA && it->IndexB == Pair.IndexB && it->Type == Pair.Type) ||
                    (it->IndexA == Pair.IndexB && it->IndexB == Pair.IndexA && it->Type == Pair.Type))
                {
                    Exists = 1;
                }
            }
            if (Exists == 0)
            {
                Pairs.push_back(Pair);
            }
        }
    }
#endif
#if CS_COLLISIONS
    for (uint32 CubeIndex = 0;
         CubeIndex < State->CubeCount;
         ++CubeIndex)
    {
        for (uint32 SphereIndex = 0;
             SphereIndex < State->SphereCount;
             ++SphereIndex)
        {
            contact_pair Pair = {};
            Pair.Type = CS;
            Pair.IndexA = CubeIndex;
            Pair.IndexB = SphereIndex;
            bool32 CollisionHappened = 0;
            std::vector<vertex> Simplex;
            CollisionHappened = CheckCollision(State, &Pair, Simplex);
            if (CollisionHappened)
            {
                bool32 Exists = 0;
                for (std::vector<contact_pair>::iterator it = Pairs.begin();
                     it != Pairs.end();
                     ++it)
                {
                    if ((it->IndexA == Pair.IndexA && it->IndexB == Pair.IndexB && it->Type == Pair.Type) ||
                        (it->IndexA == Pair.IndexB && it->IndexB == Pair.IndexA && it->Type == Pair.Type))
                    {
                        Exists = 1;
                    }
                }
                if (Exists == 0)
                {
                    Pairs.push_back(Pair);
                }
            }
        }
    }
#endif
#if SS_COLLISIONS
    // NOTE(Jovan): SS
    for (uint32 Sphere1 = 0;
         Sphere1 < State->SphereCount;
         ++Sphere1)
    {
        for (uint32 Sphere2 = Sphere1 + 1;
             Sphere2 < State->SphereCount;
             ++Sphere2)
        {
            contact_pair Pair = {};
            Pair.Type = SS;
            Pair.IndexA = Sphere1;
            Pair.IndexB = Sphere2;
            bool32 CollisionHappened = 0;
            std::vector<vertex> Simplex;
            CollisionHappened = CheckCollision(State, &Pair, Simplex);
            if (CollisionHappened)
            {
                bool32 Exists = 0;
                for (std::vector<contact_pair>::iterator it = Pairs.begin();
                     it != Pairs.end();
                     ++it)
                {
                    if ((it->IndexA == Pair.IndexA && it->IndexB == Pair.IndexB && it->Type == Pair.Type) ||
                        (it->IndexA == Pair.IndexB && it->IndexB == Pair.IndexA && it->Type == Pair.Type))
                    {
                        Exists = 1;
                    }
                }
                if (Exists == 0)
                {
                    Pairs.push_back(Pair);
                }
            }
        }
    }
#endif
}

static void
SolveConstraints(sdl_state *State, real32 dt, std::vector<contact_pair> &Pairs)
{
    for (std::vector<contact_pair>::iterator it = Pairs.begin();
         it != Pairs.end();
         ++it)
    {
        contact_pair Tmp = *it;
        Constraint(State, &Tmp, dt);
    }
}

extern "C" SIM_UPDATE_AND_RENDER(SimUpdateAndRender)
{
    sdl_state *SimState = (sdl_state *)Memory->PermanentStorage;

    // NOTE(Jovan): Init
    // -----------------
    if (!Memory->IsInitialized)
    {
        // NOTE(Jovan): Init random seed
        srand(time(0));

        // NOTE(Jovan): GJK init stuff
#if 0 
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
#endif
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
        SimState->Cubes[0].Position = glm::vec3(2.0f, 3.5f, 2.0f); //glm::vec3(2.0, 4.5f, 2.1);
        SimState->Cubes[0].V = glm::vec3(0.0);
        SimState->Cubes[0].Forces = glm::vec3(0.0);
        SimState->Cubes[0].Angles = glm::vec3(0.0);
        SimState->Cubes[0].W = glm::vec3(0.0);
        SimState->Cubes[0].Torque = glm::vec3(0.0);
        SimState->Cubes[0].Size = 1.0f;
        SimState->Cubes[0].Mass = 1.0f;
        SimState->Cubes[0].MOI = (SimState->Cubes[0].Mass / 12.0f) *
                                 (2.0f * pow(SimState->Cubes[0].Size, 2));

        SimState->Cubes[1].Model = glm::mat4(1.0);
        UpdateVertices(SimState, 1);
        SimState->Cubes[1].Position = glm::vec3(2.0, 1.0, 2.0);
        SimState->Cubes[1].V = glm::vec3(0.0);
        SimState->Cubes[1].Forces = glm::vec3(0.0);
        SimState->Cubes[1].Angles = glm::vec3(0.0);
        SimState->Cubes[1].W = glm::vec3(0.0);
        SimState->Cubes[1].Torque = glm::vec3(0.0);
        SimState->Cubes[1].Size = 1.0f;
        SimState->Cubes[1].Mass = 1.0f;
        SimState->Cubes[1].MOI = (SimState->Cubes[1].Mass / 12.0f) *
                                 (2.0f * SimState->Cubes[1].Size * SimState->Cubes[1].Size);

        SimState->Cubes[2].Model = glm::mat4(1.0);
        UpdateVertices(SimState, 2);
        SimState->Cubes[2].Position = glm::vec3(2.0, 4.5, 2.0);
        SimState->Cubes[2].V = glm::vec3(0.0);
        SimState->Cubes[2].Forces = glm::vec3(0.0);
        SimState->Cubes[2].Angles = glm::vec3(0.0);
        SimState->Cubes[2].W = glm::vec3(0.0);
        SimState->Cubes[2].Torque = glm::vec3(0.0);
        SimState->Cubes[2].Size = 1.0f;
        SimState->Cubes[2].Mass = 1.0f;
        SimState->Cubes[2].MOI = (SimState->Cubes[2].Mass / 12.0f) *
                                 (2.0f * SimState->Cubes[2].Size * SimState->Cubes[2].Size);

        // NOTE(Jovan): Sphere init
        SimState->SphereCount = 1;
        SimState->Spheres[0].Model = glm::mat4(1.0f);
        SimState->Spheres[0].Position = glm::vec3(0.1f, 1.1f, 1.1f);
        SimState->Spheres[0].V = glm::vec3(0.0f);
        SimState->Spheres[0].Forces = glm::vec3(0.0f);
        SimState->Spheres[0].Angles = glm::vec3(0.0f);
        SimState->Spheres[0].W = glm::vec3(0.0f);
        SimState->Spheres[0].Torque = glm::vec3(0.0f);
        SimState->Spheres[0].Radius = 0.25f;
        SimState->Spheres[0].Mass = 2.0f;
        SimState->Spheres[0].MOI = (2.0f / 5.0f) * SimState->Spheres[0].Mass *
                                   pow(SimState->Spheres[0].Radius, 2);

        // NOTE(Jovan): Floor init
        SimState->Floor.Model = glm::mat4(1.0);
        FloorUpdateVertices(SimState);
        SimState->Floor.Position = glm::vec3(1.2f, -0.5f, 1.0f);
        SimState->Floor.V = glm::vec3(0.0);
        SimState->Floor.Forces = glm::vec3(0.0);
        SimState->Floor.Angles = glm::vec3(0.0);
        SimState->Floor.W = glm::vec3(0.0);
        SimState->Floor.Torque = glm::vec3(0.0);
        SimState->Floor.Size = 100.0f;
        SimState->Floor.Mass = 100000.0f;
        SimState->Floor.MOI = (SimState->Floor.Mass / 12.0f) *
                              (2.0f * SimState->Floor.Size * SimState->Floor.Size);

        // NOTE(Jovan): Light init
        SimState->Light.Position = glm::vec3(3.0f);
        SimState->Light.Size = 0.2f;

        Memory->IsInitialized = 1;
    }
    // NOTE(Jovan): End init
    // ---------------------

    // NOTE(Jovan): Coordinate systems
    // -------------------------------
    Render->Projection = glm::perspective(glm::radians(45.0f),
                                          (real32)DEFAULT_WINDOW_WIDTH / (real32)DEFAULT_WINDOW_HEIGHT,
                                          0.1f,
                                          100.0f);
    glm::mat4 Model = glm::mat4(1.0f);
    Render->View = glm::mat4(1.0f);
    Render->View = glm::translate(Render->View, glm::vec3(0.0f, 0.0f, -3.0f));
    Render->View = glm::lookAt(SimState->Camera.Position, SimState->Camera.Position + SimState->Camera.Front,
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
    SimState->Pairs.clear();
    IntegrateForces(SimState, dt);
    DetectCollisions(SimState, Input, Render, dt, SimState->Pairs);
    SolveConstraints(SimState, dt, SimState->Pairs);
    IntegrateVelocities(SimState, dt);
    // NOTE(Jovan): End physics stuff
    // ------------------------------

    glm::vec3 Front;
    Front.x = cos(glm::radians(SimState->Camera.Yaw)) * cos(glm::radians(SimState->Camera.Pitch));
    Front.y = sin(glm::radians(SimState->Camera.Pitch));
    Front.z = sin(glm::radians(SimState->Camera.Yaw)) * cos(glm::radians(SimState->Camera.Pitch));
    SimState->Camera.Front = glm::normalize(Front);
    Render->View = glm::lookAt(SimState->Camera.Position, SimState->Camera.Position + SimState->Camera.Front,
                               SimState->Camera.Up);



    // NOTE(Jovan): Light drawing
    // --------------------------
    glUseProgram(Render->Shaders[CUBESH]);
    SetUniformM4(Render->Shaders[CUBESH], "View", Render->View);
    SetUniformM4(Render->Shaders[CUBESH], "Projection", Render->Projection);
    Model = glm::mat4(1.0f);
    Model = glm::translate(Model, SimState->Light.Position);
    Model = glm::scale(Model, SimState->Light.Size * glm::vec3(1.0f));
    SetUniformM4(Render->Shaders[CUBESH], "Model", Model);
    glBindVertexArray(Render->VAOs[LIGHTVAO]);
    glDrawArrays(GL_TRIANGLES, 0, 36);
    glUseProgram(0);
    glBindVertexArray(0);
    glActiveTexture(0);
    glBindTexture(GL_TEXTURE_2D, 0);
    // NOTE(Jovan): End of light drawing)

#if DRAW_WIRE
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
#else
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
#endif

    glBindVertexArray(0);

#if DRAW_FLOOR
    // NOTE(Jovan): Floor drawing
    // --------------------------

    glUseProgram(Render->Shaders[LIGHTSH]);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, Render->Textures[CHECKERBOARD]);
    SetUniformI1(Render->Shaders[LIGHTSH], "Material.Diffuse", 0);
    SetUniformF1(Render->Shaders[LIGHTSH], "Material.Shininess", 64.0f);
    SetUniformF3(Render->Shaders[LIGHTSH], "Material.Specular", 0.5f, 0.5f, 0.5f);
    SetUniformM4(Render->Shaders[LIGHTSH], "View", Render->View);
    SetUniformM4(Render->Shaders[LIGHTSH], "Projection", Render->Projection);
    SetUniformF3(Render->Shaders[LIGHTSH], "Light.Ambient", 0.2f, 0.2f, 0.2f);
    SetUniformF3(Render->Shaders[LIGHTSH], "Light.Diffuse", 0.5f, 0.5f, 0.5f);
    SetUniformF3(Render->Shaders[LIGHTSH], "Light.Specular", 1.0f, 1.0f, 1.0f);
    SetUniformF3(Render->Shaders[LIGHTSH], "ViewPos", SimState->Camera.Position.x, SimState->Camera.Position.y, SimState->Camera.Position.z);
    SetUniformF3(Render->Shaders[LIGHTSH], "Light.Position", SimState->Light.Position.x, SimState->Light.Position.y, SimState->Light.Position.z);
    SimState->Floor.Model = glm::mat4(1.0);
    SimState->Floor.Model = glm::translate(SimState->Floor.Model, SimState->Floor.Position);
    SimState->Floor.Model = glm::rotate(SimState->Floor.Model, glm::radians(SimState->Floor.Angles.x),
                                        glm::vec3(1.0f, 0.0f, 0.0f));
    SimState->Floor.Model = glm::rotate(SimState->Floor.Model, glm::radians(SimState->Floor.Angles.y),
                                        glm::vec3(0.0f, 1.0f, 0.0f));
    SimState->Floor.Model = glm::rotate(SimState->Floor.Model, glm::radians(SimState->Floor.Angles.z),
                                        glm::vec3(0.0f, 0.0f, 1.0f));
    SimState->Floor.Model = glm::scale(SimState->Floor.Model,
                                       glm::vec3(SimState->Floor.Size, 1.0, SimState->Floor.Size));
    SetUniformM4(Render->Shaders[LIGHTSH], "Model", SimState->Floor.Model);
    FloorUpdateVertices(SimState);
    glBindVertexArray(Render->VAOs[CUBEVAO]);
    glDrawArrays(GL_TRIANGLES, 0, 36);
    glBindVertexArray(0);
    glBindTexture(GL_TEXTURE_2D, 0);

    // NOTE(Jovan): End floor drawing
    // ------------------------------
#endif

    // NOTE(Jovan): Cube drawing
    // -------------------------
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, Render->Textures[CONTAINER_DIFFUSE]);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, Render->Textures[CONTAINER_SPECULAR]);
    SetUniformI1(Render->Shaders[LIGHTSH], "Material.Diffuse", 0);
    SetUniformI1(Render->Shaders[LIGHTSH], "Material.Specular", 1);
    SetUniformF1(Render->Shaders[LIGHTSH], "Material.Shininess", 64.0f);
    SetUniformF3(Render->Shaders[LIGHTSH], "Material.Specular", 0.5f, 0.5f, 0.5f);
    SetUniformM4(Render->Shaders[LIGHTSH], "View", Render->View);
    SetUniformM4(Render->Shaders[LIGHTSH], "Projection", Render->Projection);
    SetUniformF3(Render->Shaders[LIGHTSH], "Light.Ambient", 0.2f, 0.2f, 0.2f);
    SetUniformF3(Render->Shaders[LIGHTSH], "Light.Diffuse", 0.5f, 0.5f, 0.5f);
    SetUniformF3(Render->Shaders[LIGHTSH], "Light.Specular", 1.0f, 1.0f, 1.0f);
    SetUniformF3(Render->Shaders[LIGHTSH], "ViewPos", SimState->Camera.Position.x, SimState->Camera.Position.y, SimState->Camera.Position.z);
    SetUniformF3(Render->Shaders[LIGHTSH], "Light.Position", SimState->Light.Position.x, SimState->Light.Position.y, SimState->Light.Position.z);
#if DRAW_CUBES
    for (uint32 CubeIndex = 0;
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
        SetUniformM4(Render->Shaders[LIGHTSH], "Model", Model);
        glBindVertexArray(Render->VAOs[CUBEVAO]);
        glDrawArrays(GL_TRIANGLES, 0, 36);
    }
    glBindVertexArray(0);
    glBindTexture(GL_TEXTURE_2D, 0);
    // NOTE(Jovan): End cube drawing
#endif

    // NOTE(Jovan): Sphere drawing
    // ---------------------------
#if DRAW_SPHERES
    glUseProgram(Render->Shaders[LIGHTSH]);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, Render->Textures[EARTH]);
    SetUniformI1(Render->Shaders[LIGHTSH], "Material.Diffuse", 0);
    SetUniformF1(Render->Shaders[LIGHTSH], "Material.Shininess", 64.0f);
    SetUniformF3(Render->Shaders[LIGHTSH], "Material.Specular", 0.5f, 0.5f, 0.5f);
    SetUniformM4(Render->Shaders[LIGHTSH], "View", Render->View);
    SetUniformM4(Render->Shaders[LIGHTSH], "Projection", Render->Projection);
    SetUniformF3(Render->Shaders[LIGHTSH], "Light.Ambient", 0.2f, 0.2f, 0.2f);
    SetUniformF3(Render->Shaders[LIGHTSH], "Light.Diffuse", 0.5f, 0.5f, 0.5f);
    SetUniformF3(Render->Shaders[LIGHTSH], "Light.Specular", 1.0f, 1.0f, 1.0f);
    SetUniformF3(Render->Shaders[LIGHTSH], "ViewPos", SimState->Camera.Position.x, SimState->Camera.Position.y, SimState->Camera.Position.z);
    SetUniformF3(Render->Shaders[LIGHTSH], "Light.Position", SimState->Light.Position.x, SimState->Light.Position.y, SimState->Light.Position.z);
    for (uint32 SphereIndex = 0;
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
        SimState->Spheres[SphereIndex].Model = Model;
        SetUniformM4(Render->Shaders[LIGHTSH], "Model", Model);
        glBindVertexArray(Render->VAOs[SPHEREVAO]);
        glDrawElements(GL_TRIANGLES, Render->Num, GL_UNSIGNED_INT, Render->Indices);
    }
    glBindVertexArray(0);
    glBindTexture(GL_TEXTURE_2D, 0);
    // NOTE(Jovan): End of sphere drawing
#endif

    // NOTE(Jovan): Model drawing
    // --------------------------
    glUseProgram(Render->Shaders[LIGHTSH]);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, Render->Textures[AMONG_US]);
    SetUniformI1(Render->Shaders[LIGHTSH], "Material.Diffuse", 0);
    SetUniformF1(Render->Shaders[LIGHTSH], "Material.Shininess", 64.0f);
    SetUniformF3(Render->Shaders[LIGHTSH], "Material.Specular", 0.5f, 0.5f, 0.5f);
    SetUniformM4(Render->Shaders[LIGHTSH], "View", Render->View);
    SetUniformM4(Render->Shaders[LIGHTSH], "Projection", Render->Projection);
    SetUniformF3(Render->Shaders[LIGHTSH], "Light.Ambient", 0.2f, 0.2f, 0.2f);
    SetUniformF3(Render->Shaders[LIGHTSH], "Light.Diffuse", 0.5f, 0.5f, 0.5f);
    SetUniformF3(Render->Shaders[LIGHTSH], "Light.Specular", 1.0f, 1.0f, 1.0f);
    SetUniformF3(Render->Shaders[LIGHTSH], "ViewPos", SimState->Camera.Position.x, SimState->Camera.Position.y, SimState->Camera.Position.z);
    SetUniformF3(Render->Shaders[LIGHTSH], "Light.Position", SimState->Light.Position.x, SimState->Light.Position.y, SimState->Light.Position.z);
    Model = glm::mat4(1.0);
    Model = glm::translate(Model, glm::vec3(0.0f));
    Model = glm::scale(Model, glm::vec3(0.005f));
    SetUniformM4(Render->Shaders[LIGHTSH], "Model", Model);
    glBindVertexArray(Render->VAOs[MODELVAO]);
    glDrawElements(GL_TRIANGLES, Render->ModelNum, GL_UNSIGNED_INT, Render->ModelIndices);
    glBindVertexArray(0);
    glBindTexture(GL_TEXTURE_2D, 0);
    // NOTE(Jovan): End of model drawing
}
