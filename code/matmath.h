#ifndef MATMATH_H
#include <cmath>
struct vector3
{
  union
  {
    // IMPORTANT(Jovan): For some reason, keep getting error
    // when trying to initialize with (vector3){.X = x , .Y = y .Z = z}
    real32 E[3];
    struct{
      real32 X;
      real32 Y;
      real32 Z;
    };
  };
  
  vector3 operator+ (vector3 const& B)
  {
    vector3 Result = {};
    Result.X = X + B.X;
    Result.Y = Y + B.Y;
    Result.Z = Z + B.Z;
    return Result;
  };
  vector3 operator- (vector3 const& B)
  {
    vector3 Result = {};
    Result.X = X - B.X;
    Result.Y = Y - B.Y;
    Result.Z = Z - B.Z;
    return Result;
  };
  vector3 operator+= (vector3 const& B)
  {
    return *this + B;
  };
  vector3 operator-= (vector3 const& B)
  {
    return *this - B;
  };
  vector3 operator* (real32 const& B)
  {
    vector3 Result = {};
    Result.X = X * B;
    Result.Y = Y * B;
    Result.Z = Z * B;
    return Result;
  };
  real32 operator* (vector3 const& B)
  {
    real32 Result;
    Result = X*B.X + Y*B.Y + Z*B.Z;
    return Result;
  };
  real32 Norm()
  {
    real32 Result;
    Result = sqrt(X*X+Y*Y+Z*Z);
    return Result;
  };
  void Normalize()
  {
    real32 Length = this->Norm();
    *this = *this * (1.0f/Length);
  };
};

struct vector4
{
  union
  {
    real32 E[4];
    struct
    {
      real32 X;
      real32 Y;
      real32 Z;
      real32 W;
    };
  };
};

struct mat3
{
  union
  {
    real32 E[3][3];
    struct
    {
      vector3 X;
      vector3 Y;
      vector3 Z;
    };
  };
  void MakeIdent()
  {
    for(int Row = 0; Row < 3; ++Row)
      for(int Col = 0; Col < 3; ++Col)
	if(Row == Col)
	  E[Row][Col] = 1.0f;
  };
};

struct mat4
{
  union
  {
    real32 E[4][4];
    struct
    {
      vector4 X;
      vector4 Y;
      vector4 Z;
      vector4 W;
    };
  };
  void MakeIdent()
  {
    for(int Row = 0; Row < 4; ++Row)
      for(int Col = 0; Col < 4; ++Col)
	if(Row == Col)
	  E[Row][Col] = 1.0f;
  };
};

#define MATMATH_H
#endif
