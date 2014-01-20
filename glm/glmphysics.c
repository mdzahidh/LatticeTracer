// GLM - a math library for OpenGL
// Programmed by Markus Ilmola
// markus.ilmola@pp.inet.fi
// http://markus_ilmola.tripod.com
// http://personal.inet.fi/koti/markus.ilmola

#include "glm.h"
#include "glmvector.h"

// Temp needs size members
void glmSolveOdeEulerf(float *result, const float *initialValue, DWORD size, float t,
                       void (*derive)(float *, const float *, float), float *temp)
{
  int i;

  derive(temp, initialValue, 0);

  for (i=0; i<size; i++) result[i]=initialValue[i]+t*temp[i];
}

void glmSolveOdeEulerd(double *result, const double *initialValue, DWORD size, double t,
                       void (*derive)(double *, const double *, double), double *temp)
{
  int i;

  derive(temp, initialValue, 0);

  for (i=0; i<size; i++) result[i]=initialValue[i]+t*temp[i];
}

// Temp needs 2*size members
void glmSolveOdeMidpointf(float *result, const float *initialValue, DWORD size, float t,
                          void (*derive)(float *, const float *, float), float *temp)
{
  int i;
  float *temp2;

  temp2=temp+size;

  derive(temp, initialValue, 0);
  for (i=0; i<size; i++) temp2[i]=initialValue[i]+(t/2.0)*temp[i];

  derive(temp, temp2, t/2.0);

  for (i=0; i<size; i++) result[i]=initialValue[i]+t*temp[i];
}

void glmSolveOdeMidpointd(double *result, const double *initialValue, DWORD size, double t,
                          void (*derive)(double *, const double *, double), double *temp)
{
  int i;
  double *temp2;

  temp2=temp+size;

  derive(temp, initialValue, 0);
  for (i=0; i<size; i++) temp2[i]=initialValue[i]+(t/2.0)*temp[i];

  derive(temp, temp2, t/2.0);

  for (i=0; i<size; i++) result[i]=initialValue[i]+t*temp[i];
}


// temp needs 6*size members
void glmSolveOdeRungeKuttaf(float *result, const float *initialValue, DWORD size, float t,
                            void (*derive)(float *, const float *, float), float *temp)
{
  int i;
  float *temp2, *k1, *k2, *k3, *k4;

  temp2=temp+size;
  k1=temp+2*size;
  k2=temp+3*size;
  k3=temp+4*size;
  k4=temp+5*size;

  derive(temp, initialValue, 0);
  for (i=0; i<size; i++) k1[i]=t*temp[i];

  for (i=0; i<size; i++) temp2[i]=initialValue[i]+0.5*k1[i];
  derive(temp, temp2, 0.5*t);
  for (i=0; i<size; i++) k2[i]=t*temp[i];

  for (i=0; i<size; i++) temp2[i]=initialValue[i]+0.5*k2[i];
  derive(temp, temp2, 0.5*t);
  for (i=0; i<size; i++) k3[i]=t*temp[i];

  for (i=0; i<size; i++) temp2[i]=initialValue[i]+k3[i];
  derive(temp, temp2, t);
  for (i=0; i<size; i++) k4[i]=t*temp[i];

  for (i=0; i<size; i++)
    result[i]=initialValue[i]+(1.0/6.0)*k1[i]+(1.0/3.0)*k2[i]+(1.0/3.0)*k3[i]+(1.0/6.0)*k4[i];
}

// temp needs 6*size members
void glmSolveOdeRungeKuttad(double *result, const double *initialValue, DWORD size, double t,
                            void (*derive)(double *, const double *, double), double *temp)
{
  int i;
  double *temp2, *k1, *k2, *k3, *k4;

  temp2=temp+size;
  k1=temp+2*size;
  k2=temp+3*size;
  k3=temp+4*size;
  k4=temp+5*size;

  derive(temp, initialValue, 0);
  for (i=0; i<size; i++) k1[i]=t*temp[i];

  for (i=0; i<size; i++) temp2[i]=initialValue[i]+0.5*k1[i];
  derive(temp, temp2, 0.5*t);
  for (i=0; i<size; i++) k2[i]=t*temp[i];

  for (i=0; i<size; i++) temp2[i]=initialValue[i]+0.5*k2[i];
  derive(temp, temp2, 0.5*t);
  for (i=0; i<size; i++) k3[i]=t*temp[i];

  for (i=0; i<size; i++) temp2[i]=initialValue[i]+k3[i];
  derive(temp, temp2, t);
  for (i=0; i<size; i++) k4[i]=t*temp[i];

  for (i=0; i<size; i++)
    result[i]=initialValue[i]+(1.0/6.0)*k1[i]+(1.0/3.0)*k2[i]+(1.0/3.0)*k3[i]+(1.0/6.0)*k4[i];
}



float glmComputeImpulse2f(const float v1[2], const float v2[2], float invM1, float invM2,
                           const float normal[2], float e)
{
  float d;
  float v[2];

  d=glmDotProduct2f(normal, normal);
  d*=invM1+invM2;

  glmMakeVector2f(v, v1, v2);

  return -((1.0f+e)*glmDotProduct2f(v, normal))/d;
}

double glmComputeImpulse2d(const double v1[2], const double v2[2], float invM1, float invM2,
                           const double normal[2], double e)
{
  double d;
  double v[2];

  d=glmDotProduct2d(normal, normal);
  d*=invM1+invM2;

  glmMakeVector2d(v, v1, v2);

  return -((1.0f+e)*glmDotProduct2d(v, normal))/d;
}

float glmComputeImpulse3f(const float v1[3], const float v2[3], float invM1, float invM2,
                           const float normal[3], float e)
{
  float d;
  float v[3];

  d=glmDotProduct3f(normal, normal);
  d*=invM1+invM2;

  glmMakeVector3f(v, v1, v2);

  return -((1.0f+e)*glmDotProduct3f(v, normal))/d;
}

float glmComputeImpulse3d(const double v1[3], const double v2[3], double invM1, double invM2,
                           const double normal[3], double e)
{
  double d;
  double v[3];

  d=glmDotProduct3d(normal, normal);
  d*=invM1+invM2;

  glmMakeVector3d(v, v1, v2);

  return -((1.0f+e)*glmDotProduct3d(v, normal))/d;
}

void glmGenerateSkewSymmetricMatrix3f(float m[9], const float v[3])
{
  m[0]=0;     m[3]=-v[2]; m[6]=v[1];
  m[1]=v[2];  m[4]=0;     m[7]=-v[0];
  m[2]=-v[1]; m[5]=v[0];  m[8]=0;
}

void glmGenerateSkewSymmetricMatrix3d(double m[9], const double v[3])
{
  m[0]=0;     m[3]=-v[2]; m[6]=v[1];
  m[1]=v[2];  m[4]=0;     m[7]=-v[0];
  m[2]=-v[1]; m[5]=v[0];  m[8]=0;
}
