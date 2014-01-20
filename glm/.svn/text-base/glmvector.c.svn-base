// GLM - a math library for OpenGL
// Programmed by Markus Ilmola
// markus.ilmola@pp.inet.fi
// http://markus_ilmola.tripod.com
// http://personal.inet.fi/koti/markus.ilmola

// Vector stuff

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "glm.h"

const float  GLM_ORIGIN_2F[2]={0, 0};
const double GLM_ORIGIN_2D[2]={0, 0};
const float  GLM_ORIGIN_3F[3]={0, 0, 0};
const double GLM_ORIGIN_3D[3]={0, 0, 0};
const float  GLM_ORIGIN_4F[4]={0, 0, 0, 0};
const double GLM_ORIGIN_4D[4]={0, 0, 0, 0};

void glmCopyVector2f(const float source[2], float target[2])
{
  target[0]=source[0];
  target[1]=source[1];
}

void glmCopyVector2d(const double source[2], double target[2])
{
  target[0]=source[0];
  target[1]=source[1];
}

void glmCopyVector3f(const float source[3], float target[3])
{
  target[0]=source[0];
  target[1]=source[1];
  target[2]=source[2];
}

void glmCopyVector3d(const double source[3], double target[3])
{
  target[0]=source[0];
  target[1]=source[1];
  target[2]=source[2];
}

void glmCopyVector4f(const float source[4], float target[4])
{
  target[0]=source[0];
  target[1]=source[1];
  target[2]=source[2];
  target[3]=source[3];
}

void glmCopyVector4d(const double source[4], double target[4])
{
  target[0]=source[0];
  target[1]=source[1];
  target[2]=source[2];
  target[3]=source[3];
}

void glmCopyVectorf(const float *source, float *target, DWORD size)
{
  memcpy(target, source, size*sizeof(float));
}

void glmCopyVectord(const double *source, double *target, DWORD size)
{
  memcpy(target, source, size*sizeof(double));
}

void glmMakeVector2f(float result[2], const float p1[2], const float p2[2])
{
  result[0]=p2[0]-p1[0];
  result[1]=p2[1]-p1[1];
}

void glmMakeVector2d(double result[2], const double p1[2], const double p2[2])
{
  result[0]=p2[0]-p1[0];
  result[1]=p2[1]-p1[1];
}

void glmMakeVector3f(float result[3], const float p1[3], const float p2[3])
{
  result[0]=p2[0]-p1[0];
  result[1]=p2[1]-p1[1];
  result[2]=p2[2]-p1[2];
}

void glmMakeVector3d(double result[3], const double p1[3], const double p2[3])
{
  result[0]=p2[0]-p1[0];
  result[1]=p2[1]-p1[1];
  result[2]=p2[2]-p1[2];
}

void glmMakeVector4f(float result[4], const float p1[4], const float p2[4])
{
  result[0]=p2[0]-p1[0];
  result[1]=p2[1]-p1[1];
  result[2]=p2[2]-p1[2];
  result[3]=p2[3]-p1[3];
}

void glmMakeVector4d(double result[4], const double p1[4], const double p2[4])
{
  result[0]=p2[0]-p1[0];
  result[1]=p2[1]-p1[1];
  result[2]=p2[2]-p1[2];
  result[3]=p2[3]-p1[3];
}

void glmMakeVectorf(float *result, const float *p1, const float *p2, DWORD size)
{
  DWORD i;
  for (i=0; i<size; i++) result[i]=p2[i]-p1[i];
}

void glmMakeVectord(double *result, const double *p1, const double *p2, DWORD size)
{
  DWORD i;
  for (i=0; i<size; i++) result[i]=p2[i]-p1[i];
}

float glmDotProduct2f(const float v1[2], const float v2[2])
{
  return v1[0]*v2[0]+v1[1]*v2[1];
}

double glmDotProduct2d(const double v1[2], const double v2[2])
{
  return v1[0]*v2[0]+v1[1]*v2[1];
}

float glmDotProduct3f(const float v1[3], const float v2[3])
{
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

double glmDotProduct3d(const double v1[3], const double v2[3])
{
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

float glmDotProduct4f(const float v1[4], const float v2[4])
{
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3];
}

float glmDotProduct4d(const double v1[4], const double v2[4])
{
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3];
}

float glmDotProductf(const float *v1, const float *v2, DWORD size)
{
  float sum=0;
  DWORD i;
  for (i=0; i<size; i++) sum+=v1[i]*v2[i];
  return sum;
}

double glmDotProductd(const double *v1, const double *v2, DWORD size)
{
  double sum=0;
  DWORD i;
  for (i=0; i<size; i++) sum+=v1[i]*v2[i];
  return sum;
}

void glmCrossProduct2f(float result[2], const float v1[2])
{
  result[0]=v1[1];
  result[1]=-v1[0];
}

void glmCrossProduct2d(double result[2], const double v1[2])
{
  result[0]=v1[1];
  result[1]=-v1[0];
}

void glmCrossProduct3f(float result[3], const float v1[3], const float v2[3])
{
  result[0]=v1[1]*v2[2]-v1[2]*v2[1];
  result[1]=v1[2]*v2[0]-v1[0]*v2[2];
  result[2]=v1[0]*v2[1]-v1[1]*v2[0];
}

void glmCrossProduct3d(double result[3], const double v1[3], const double v2[3])
{
  result[0]=v1[1]*v2[2]-v1[2]*v2[1];
  result[1]=v1[2]*v2[0]-v1[0]*v2[2];
  result[2]=v1[0]*v2[1]-v1[1]*v2[0];
}

void glmCrossProduct4f(float result[4], const float v1[4], const float v2[4], const float v3[4])
{
  result[0]=-v1[2]*v2[1]*v3[3] + v1[1]*v2[2]*v3[3] - v1[3]*v2[2]*v3[1]
            +v1[2]*v2[3]*v3[1] + v1[3]*v2[1]*v3[2] - v1[1]*v2[3]*v3[2];
  result[1]= v1[2]*v2[0]*v3[3] - v1[0]*v2[2]*v3[3] + v1[3]*v2[2]*v3[0]
            -v1[2]*v2[3]*v3[0] - v1[3]*v2[0]*v3[2] + v1[0]*v2[3]*v3[2];
  result[2]=-v1[1]*v2[0]*v3[3] + v1[0]*v2[1]*v3[3] - v1[3]*v2[1]*v3[0]
            +v1[1]*v2[3]*v3[0] + v1[3]*v2[0]*v3[1] - v1[0]*v2[3]*v3[1];
  result[3]= v1[2]*v2[1]*v3[0] - v1[1]*v2[2]*v3[0] - v1[2]*v2[0]*v3[1]
            +v1[0]*v2[2]*v3[1] + v1[1]*v2[0]*v3[2] - v1[0]*v2[1]*v3[2];
}

void glmCrossProduct4d(double result[4], const double v1[4], const double v2[4], const double v3[4])
{
  result[0]=-v1[2]*v2[1]*v3[3] + v1[1]*v2[2]*v3[3] - v1[3]*v2[2]*v3[1]
            +v1[2]*v2[3]*v3[1] + v1[3]*v2[1]*v3[2] - v1[1]*v2[3]*v3[2];
  result[1]= v1[2]*v2[0]*v3[3] - v1[0]*v2[2]*v3[3] + v1[3]*v2[2]*v3[0]
            -v1[2]*v2[3]*v3[0] - v1[3]*v2[0]*v3[2] + v1[0]*v2[3]*v3[2];
  result[2]=-v1[1]*v2[0]*v3[3] + v1[0]*v2[1]*v3[3] - v1[3]*v2[1]*v3[0]
            +v1[1]*v2[3]*v3[0] + v1[3]*v2[0]*v3[1] - v1[0]*v2[3]*v3[1];
  result[3]= v1[2]*v2[1]*v3[0] - v1[1]*v2[2]*v3[0] - v1[2]*v2[0]*v3[1]
            +v1[0]*v2[2]*v3[1] + v1[1]*v2[0]*v3[2] - v1[0]*v2[1]*v3[2];
}

/*
void glmCrossProductf(float *result, DWORD size, float *first, ...)
{
  va_list marker;
  DWORD i;
  float *vector;
  float *m;
  DWORD i, j;

  m=(float *)malloc(size*size*sizeof(float));
  if (!m) return;

  va_start(marker, first);
  for (i=0; i<size-1; i++)
  {
    vector=va_arg(marker, float *);
    for (j=0; j<size; j++) m[i+1+size*j];
  }
  va_end(marker);

  free(m);
}
*/

float glmNormalizeVector2f(float v[2])
{
  float length;
  length=sqrt(v[0]*v[0]+v[1]*v[1]);
  if (length>0)
  {
    v[0]/=length;
    v[1]/=length;
  }
  return length;
}

double glmNormalizeVector2d(double v[2])
{
  double length;
  length=sqrt(v[0]*v[0]+v[1]*v[1]);
  if (length>0)
  {
    v[0]/=length;
    v[1]/=length;
  }
  return length;
}

float glmNormalizeVector3f(float v[3])
{
  float length;
  length=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if (length>0)
  {
    v[0]/=length;
    v[1]/=length;
    v[2]/=length;
  }
  return length;
}

double glmNormalizeVector3d(double v[3])
{
  double length;
  length=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if (length>0)
  {
    v[0]/=length;
    v[1]/=length;
    v[2]/=length;
  }
  return length;
}

float glmNormalizeVector4f(float v[4])
{
  float length;
  length=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
  if (length>0)
  {
    v[0]/=length;
    v[1]/=length;
    v[2]/=length;
    v[3]/=length;
  }
  return length;
}

double glmNormalizeVector4d(double v[4])
{
  double length;
  length=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
  if (length>0)
  {
    v[0]/=length;
    v[1]/=length;
    v[2]/=length;
    v[3]/=length;
  }
  return length;
}

float glmNormalizeVectorf(float *v, DWORD size)
{
  float length=0;
  int i;
  for (i=0; i<size; i++) length+=v[i]*v[i];
  length=sqrt(length);
  if (length>0)
  {
    for (i=0; i<size; i++) v[i]/=length;
  }
  return length;
}

double glmNormalizeVectord(double *v, DWORD size)
{
  double length=0;
  int i;
  for (i=0; i<size; i++) length+=v[i]*v[i];
  length=sqrt(length);
  if (length>0)
  {
    for (i=0; i<size; i++) v[i]/=length;
  }
  return length;
}

float glmComputeVectorLength2f(const float v[2])
{
  return sqrt(v[0]*v[0]+v[1]*v[1]);
}

double glmComputeVectorLength2d(const double v[2])
{
  return sqrt(v[0]*v[0]+v[1]*v[1]);
}

float glmComputeVectorLength3f(const float v[3])
{
  return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

double glmComputeVectorLength3d(const double v[3])
{
  return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

float glmComputeVectorLength4f(const float v[4])
{
  return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
}

double glmComputeVectorLength4d(const double v[4])
{
  return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
}

float glmComputeVectorLengthf(const float *v, DWORD size)
{
  float lenght=0;
  DWORD i;
  for (i=0; i<size; i++) lenght+=v[i]*v[i];
  return sqrt(lenght);
}

double glmComputeVectorLengthd(const double *v, DWORD size)
{
  double lenght=0;
  DWORD i;
  for (i=0; i<size; i++) lenght+=v[i]*v[i];
  return sqrt(lenght);
}

float glmComputeVectorLength2fPow2(const float v[2])
{
  return v[0]*v[0]+v[1]*v[1];
}

double glmComputeVectorLength2dPow2(const double v[2])
{
  return v[0]*v[0]+v[1]*v[1];
}

float glmComputeVectorLength3fPow2(const float v[3])
{
  return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
}

double glmComputeVectorLength3dPow2(const double v[3])
{
  return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
}

float glmComputeVectorLength4fPow2(const float v[4])
{
  return v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3];
}

double glmComputeVectorLength4dPow2(const float v[4])
{
  return v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3];
}

float glmComputeVectorLengthfPow2(const float *v, DWORD size)
{
  float lenght=0;
  DWORD i;
  for (i=0; i<size; i++) lenght+=v[i]*v[i];
  return lenght;
}

double glmComputeVectorLengthdPow2(const double *v, DWORD size)
{
  double lenght=0;
  DWORD i;
  for (i=0; i<size; i++) lenght+=v[i]*v[i];
  return lenght;
}

float glmComputeNormalizedVectorAngle2fRad(const float v1[2], const float v2[2])
{
  return acos(glmDotProduct2f(v1, v2));
}

float glmComputeNormalizedVectorAngle2fDeg(const float v1[2], const float v2[2])
{
  return glmDegf(acos(glmDotProduct2f(v1, v2)));
}

double glmComputeNormalizedVectorAngle2dRad(const double v1[2], const double v2[2])
{
  return acos(glmDotProduct2d(v1, v2));
}

double glmComputeNormalizedVectorAngle2dDeg(const double v1[2], const double v2[2])
{
  return glmDegd(acos(glmDotProduct2d(v1, v2)));
}

float glmComputeNormalizedVectorAngle3fRad(const float v1[3], const float v2[3])
{
  return acos(glmDotProduct3f(v1, v2));
}

float glmComputeNormalizedVectorAngle3fDeg(const float v1[3], const float v2[3])
{
  return glmDegf(acos(glmDotProduct3f(v1, v2)));
}

double glmComputeNormalizedVectorAngle3dRad(const double v1[3], const double v2[3])
{
  return acos(glmDotProduct3d(v1, v2));
}

double glmComputeNormalizedVectorAngle3dDeg(const double v1[3], const double v2[3])
{
  return glmDegd(acos(glmDotProduct3d(v1, v2)));
}

float glmComputeNormalizedVectorAngle4fRad(const float v1[4], const float v2[4])
{
  return acos(glmDotProduct4f(v1, v2));
}

double glmComputeNormalizedVectorAngle4dRad(const double v1[4], const double v2[4])
{
  return acos(glmDotProduct4d(v1, v2));
}

float glmComputeNormalizedVectorAngle4fDeg(const float v1[4], const float v2[4])
{
  return glmDegf(acos(glmDotProduct4f(v1, v2)));
}

double glmComputeNormalizedVectorAngle4dDeg(const double v1[4], const double v2[4])
{
  return glmDegd(acos(glmDotProduct4d(v1, v2)));
}

float glmComputeNormalizedVectorAnglefRad(const float *v1, const float *v2, DWORD size)
{
  return acos(glmDotProductf(v1, v2, size));
}

float glmComputeNormalizedVectorAnglefDeg(const float *v1, const float *v2, DWORD size)
{
  return glmDegf(acos(glmDotProductf(v1, v2, size)));
}

double glmComputeNormalizedVectorAngledRad(const double *v1, const double *v2, DWORD size)
{
  return acos(glmDotProductd(v1, v2, size));
}

double glmComputeNormalizedVectorAngledDeg(const double *v1, const double *v2, DWORD size)
{
  return glmDegd(acos(glmDotProductd(v1, v2, size)));
}

float glmComputeVectorAngle2fRad(const float v1[2], const float v2[2])
{
  return acos(glmDotProduct2f(v1, v2)/
         (glmComputeVectorLength2f(v1)*glmComputeVectorLength2f(v2)));
}

float glmComputeVectorAngle2fDeg(const float v1[2], const float v2[2])
{
  return glmDegf(acos(glmDotProduct2f(v1, v2)/
         (glmComputeVectorLength2f(v1)*glmComputeVectorLength2f(v2))));
}

double glmComputeVectorAngle2dRad(const double v1[2], const double v2[2])
{
  return acos(glmDotProduct2d(v1, v2)/
         (glmComputeVectorLength2d(v1)*glmComputeVectorLength2d(v2)));
}

double glmComputeVectorAngle2dDeg(const double v1[2], const double v2[2])
{
  return glmDegd(acos(glmDotProduct2d(v1, v2)/
         (glmComputeVectorLength2d(v1)*glmComputeVectorLength2d(v2))));
}

float glmComputeVectorAngle3fRad(const float v1[3], const float v2[3])
{
  return acos(glmDotProduct3f(v1, v2)/
         (glmComputeVectorLength3f(v1)*glmComputeVectorLength3f(v2)));
}

float glmComputeVectorAngle3fDeg(const float v1[3], const float v2[3])
{
  return glmDegf(acos(glmDotProduct3f(v1, v2)/
         (glmComputeVectorLength3f(v1)*glmComputeVectorLength3f(v2))));
}

double glmComputeVectorAngle3dRad(const double v1[3], const double v2[3])
{
  return acos(glmDotProduct3d(v1, v2)/
         (glmComputeVectorLength3d(v1)*glmComputeVectorLength3d(v2)));
}

double glmComputeVectorAngle3dDeg(const double v1[3], const double v2[3])
{
  return glmDegd(acos(glmDotProduct3d(v1, v2)/
         (glmComputeVectorLength3d(v1)*glmComputeVectorLength3d(v2))));
}

float glmComputeVectorAngle4fRad(const float v1[4], const float v2[4])
{
  return acos(glmDotProduct4f(v1, v2)/
         (glmComputeVectorLength4f(v1)*glmComputeVectorLength4f(v2)));
}

float glmComputeVectorAngle4fDeg(const float v1[4], const float v2[4])
{
  return glmDegf(acos(glmDotProduct4f(v1, v2)/
         (glmComputeVectorLength4f(v1)*glmComputeVectorLength4f(v2))));
}

double glmComputeVectorAngle4dRad(const double v1[4], const double v2[4])
{
  return acos(glmDotProduct4d(v1, v2)/
         (glmComputeVectorLength4d(v1)*glmComputeVectorLength4d(v2)));
}

double glmComputeVectorAngle4dDeg(const double v1[4], const double v2[4])
{
  return glmDegd(acos(glmDotProduct4d(v1, v2)/
         (glmComputeVectorLength4d(v1)*glmComputeVectorLength4d(v2))));
}

float glmComputeVectorAnglefRad(const float *v1, const float *v2, DWORD size)
{
  return acos(glmDotProductf(v1, v2, size)/
         (glmComputeVectorLengthf(v1, size)*glmComputeVectorLengthf(v2, size)));
}

float glmComputeVectorAnglefDeg(const float *v1, const float *v2, DWORD size)
{
  return glmDegf(acos(glmDotProductf(v1, v2, size)/
         (glmComputeVectorLengthf(v1, size)*glmComputeVectorLengthf(v2, size))));
}

double glmComputeVectorAngledRad(const double *v1, const double *v2, DWORD size)
{
  return acos(glmDotProductd(v1, v2, size)/
         (glmComputeVectorLengthd(v1, size)*glmComputeVectorLengthd(v2, size)));
}

double glmComputeVectorAngledDeg(const double *v1, const double *v2, DWORD size)
{
  return glmDegd(acos(glmDotProductd(v1, v2, size)/
         (glmComputeVectorLengthd(v1, size)*glmComputeVectorLengthd(v2, size))));
}

void glmAddVector2f(float result[2], const float v1[2], const float v2[2])
{
  result[0]=v1[0]+v2[0];
  result[1]=v1[1]+v2[1];
}

void glmAddVector2d(double result[2], const double v1[2], const double v2[2])
{
  result[0]=v1[0]+v2[0];
  result[1]=v1[1]+v2[1];
}

void glmAddVector3f(float result[3], const float v1[3], const float v2[3])
{
  result[0]=v1[0]+v2[0];
  result[1]=v1[1]+v2[1];
  result[2]=v1[2]+v2[2];
}

void glmAddVector3d(double result[3], const double v1[3], const double v2[3])
{
  result[0]=v1[0]+v2[0];
  result[1]=v1[1]+v2[1];
  result[2]=v1[2]+v2[2];
}

void glmAddVector4f(float result[4], const float v1[4], const float v2[4])
{
  result[0]=v1[0]+v2[0];
  result[1]=v1[1]+v2[1];
  result[2]=v1[2]+v2[2];
  result[3]=v1[3]+v2[3];
}

void glmAddVector4d(double result[4], const double v1[4], const double v2[4])
{
  result[0]=v1[0]+v2[0];
  result[1]=v1[1]+v2[1];
  result[2]=v1[2]+v2[2];
  result[3]=v1[3]+v2[3];
}

void glmAddVectorf(float *result, const float *v1, const float *v2, DWORD size)
{
  DWORD i;
  for (i=0; i<size; i++) result[i]=v1[i]+v2[i];
}

void glmAddVectord(double *result, const double *v1, const double *v2, DWORD size)
{
  DWORD i;
  for (i=0; i<size; i++) result[i]=v1[i]+v2[i];
}

void glmSubtractVector2f(float result[2], const float v1[2], const float v2[2])
{
  result[0]=v1[0]-v2[0];
  result[1]=v1[1]-v2[1];
}

void glmSubtractVector2d(double result[2], const double v1[2], const double v2[2])
{
  result[0]=v1[0]-v2[0];
  result[1]=v1[1]-v2[1];
}

void glmSubtractVector3f(float result[3], const float v1[3], const float v2[3])
{
  result[0]=v1[0]-v2[0];
  result[1]=v1[1]-v2[1];
  result[2]=v1[2]-v2[2];
}

void glmSubtractVector3d(double result[3], const double v1[3], const double v2[3])
{
  result[0]=v1[0]-v2[0];
  result[1]=v1[1]-v2[1];
  result[2]=v1[2]-v2[2];
}

void glmSubtractVector4f(float result[4], const float v1[4], const float v2[4])
{
  result[0]=v1[0]-v2[0];
  result[1]=v1[1]-v2[1];
  result[2]=v1[2]-v2[2];
  result[3]=v1[3]-v2[3];
}

void glmSubtractVector4d(double result[4], const double v1[4], const double v2[4])
{
  result[0]=v1[0]-v2[0];
  result[1]=v1[1]-v2[1];
  result[2]=v1[2]-v2[2];
  result[3]=v1[3]-v2[3];
}

void glmSubtractVectorf(float *result, const float *v1, const float *v2, DWORD size)
{
  DWORD i;
  for (i=0; i<size; i++) result[i]=v1[i]-v2[i];
}

void glmSubtractVectord(double *result, const double *v1, const double *v2, DWORD size)
{
  DWORD i;
  for (i=0; i<size; i++) result[i]=v1[i]-v2[i];
}

void glmReflectVector2f(float result[2], const float v[2], const float normal[2])
{
  float temp;
  temp=2.0*(v[0]*normal[0]+v[1]*normal[1]);
  result[0]=v[0]-normal[0]*temp;
  result[1]=v[1]-normal[1]*temp;
}

void glmReflectVector2d(double result[2], const double v[2], const double normal[2])
{
  double temp;
  temp=2.0*(v[0]*normal[0]+v[1]*normal[1]);
  result[0]=v[0]-normal[0]*temp;
  result[1]=v[1]-normal[1]*temp;
}

void glmReflectVector3f(float result[3], const float v[3], const float normal[3])
{
  float temp;
  temp=2.0*(v[0]*normal[0]+v[1]*normal[1]+v[2]*normal[2]);
  result[0]=v[0]-normal[0]*temp;
  result[1]=v[1]-normal[1]*temp;
  result[2]=v[2]-normal[2]*temp;
}

void glmReflectVector3d(double result[3], const double v[3], const double normal[3])
{
  double temp;
  temp=2.0*(v[0]*normal[0]+v[1]*normal[1]+v[2]*normal[2]);
  result[0]=v[0]-normal[0]*temp;
  result[1]=v[1]-normal[1]*temp;
  result[2]=v[2]-normal[2]*temp;
}

void glmReflectVector4f(float result[4], const float v[4], const float normal[4])
{
  float temp;
  temp=2.0*(v[0]*normal[0]+v[1]*normal[1]+v[2]*normal[2]+v[3]*normal[3]);
  result[0]=v[0]-normal[0]*temp;
  result[1]=v[1]-normal[1]*temp;
  result[2]=v[2]-normal[2]*temp;
  result[3]=v[3]-normal[3]*temp;
}

void glmReflectVector4d(double result[4], const double v[4], const double normal[4])
{
  double temp;
  temp=2.0*(v[0]*normal[0]+v[1]*normal[1]+v[2]*normal[2]+v[3]*normal[3]);
  result[0]=v[0]-normal[0]*temp;
  result[1]=v[1]-normal[1]*temp;
  result[2]=v[2]-normal[2]*temp;
  result[3]=v[3]-normal[3]*temp;
}

void glmReflectVectorf(float *result, const float *v, const float *normal, DWORD size)
{
  float temp;
  DWORD i;
  temp=2.0*glmDotProductf(v, normal, size);
  for (i=0; i<size; i++) result[i]=v[i]-normal[i]*temp;
}

void glmReflectVectord(double *result, const double *v, const double *normal, DWORD size)
{
  float temp;
  DWORD i;
  temp=2.0*glmDotProductd(v, normal, size);
  for (i=0; i<size; i++) result[i]=v[i]-normal[i]*temp;
}

void glmGenerateRandomVector2f(float v[2])
{
  float a;
  a=glmRandomNumberf(0, GLM_2PI);
  v[0]=sin(a);
  v[1]=cos(a);
}

void glmGenerateRandomVector2d(double v[2])
{
  double a;
  a=glmRandomNumberd(0, GLM_2PI);
  v[0]=sin(a);
  v[1]=cos(a);
}

void glmGenerateRandomVector3f(float v[3])
{
  float a, r;
  v[2]=glmRandomNumberf(-1.0f, 1.0f);
  a=glmRandomNumberf(0.0f, GLM_2PI);
  r=sqrt(1.0f - v[2]*v[2]);
  v[0]=r*cos(a);
  v[1]=r*sin(a);
}

void glmGenerateRandomVector3d(double v[3])
{
  double a, r;
  v[2]=glmRandomNumberd(-1.0, 1.0);
  a=glmRandomNumberd(0.0, GLM_2PI);
  r=sqrt(1.0 - v[2]*v[2]);
  v[0]=r*cos(a);
  v[1]=r*sin(a);
}

void glmComputeLineNormalf(float normal[2], const float c1[2], const float c2[2])
{
  float v[2];

  glmMakeVector2f(v, c1, c2);
  glmCrossProduct2f(normal, v);
  glmNormalizeVector2f(normal);
}

void glmComputeLineNormald(double normal[2], const double c1[2], const double c2[2])
{
  double v[2];

  glmMakeVector2d(v, c1, c2);
  glmCrossProduct2d(normal, v);
  glmNormalizeVector2d(normal);
}

void glmComputeTriangleNormalf(float normal[3], const float c1[3], const float c2[3], const float c3[3])
{
  float v1[3], v2[3];

  glmMakeVector3f(v1, c1, c2);
  glmMakeVector3f(v2, c1, c3);
  glmCrossProduct3f(normal, v1, v2);
  glmNormalizeVector3f(normal);
}

void glmComputeTriangleNormald(double normal[3], const double c1[3], const double c2[3], const double c3[3])
{
  double v1[3], v2[3];

  glmMakeVector3d(v1, c1, c2);
  glmMakeVector3d(v2, c1, c3);
  glmCrossProduct3d(normal, v1, v2);
  glmNormalizeVector3d(normal);
}

void glmMultVectorByScalar2f(float v[2], float a)
{
  v[0]*=a;
  v[1]*=a;
}

void glmMultVectorByScalar2d(double v[2], double a)
{
  v[0]*=a;
  v[1]*=a;
}

void glmMultVectorByScalar3f(float v[3], float a)
{
  v[0]*=a;
  v[1]*=a;
  v[2]*=a;
}

void glmMultVectorByScalar4f(float v[4], float a)
{
  v[0]*=a;
  v[1]*=a;
  v[2]*=a;
  v[3]*=a;
}

void glmMultVectorByScalar3d(double v[3], double a)
{
  v[0]*=a;
  v[1]*=a;
  v[2]*=a;  
}

void glmMultVectorByScalarf(float *v, float a, DWORD size)
{
  DWORD i;
  for (i=0; i<size; i++) v[i]*=a;
}

void glmMultVectorByScalard(double *v, double a, DWORD size)
{
  DWORD i;
  for (i=0; i<size; i++) v[i]*=a;
}

void glmProjectVector2f(float result[2], const float v[2], const float normal[2])
{
  float temp;
  float tempVector[2];
  temp=glmDotProduct2f(v, normal);
  glmCopyVector2f(normal, tempVector);
  glmMultVectorByScalar2f(tempVector, temp);
  glmSubtractVector2f(result, v, tempVector);
}

void glmProjectVector2d(double result[2], const double v[2], const double normal[2])
{
  double temp;
  double tempVector[2];
  temp=glmDotProduct2d(v, normal);
  glmCopyVector2d(normal, tempVector);
  glmMultVectorByScalar2d(tempVector, temp);
  glmSubtractVector2d(result, v, tempVector);
}

void glmProjectVector3f(float result[3], const float v[3], const float normal[3])
{
  float temp;
  float tempVector[3];
  temp=glmDotProduct3f(v, normal);
  glmCopyVector3f(normal, tempVector);
  glmMultVectorByScalar3f(tempVector, temp);
  glmSubtractVector3f(result, v, tempVector);
}

void glmProjectVector3d(double result[3], const double v[3], const double normal[3])
{
  double temp;
  double tempVector[3];
  temp=glmDotProduct3d(v, normal);
  glmCopyVector3d(normal, tempVector);
  glmMultVectorByScalar3d(tempVector, temp);
  glmSubtractVector3d(result, v, tempVector);
}

void glmPrintVector2f(const float v[2])
{
  printf("[%f,%f]", v[0], v[1]);
}

void glmPrintVector2d(const double v[2])
{
  printf("[%f,%f]", v[0], v[1]);
}

void glmPrintVector3f(const float v[3])
{
  printf("[%f,%f,%f]", v[0], v[1], v[2]);
}

void glmPrintVector3d(const double v[3])
{
  printf("[%f,%f,%f]", v[0], v[1], v[2]);
}

void glmPrintVector4f(const float v[4])
{
  printf("[%f,%f,%f,%f]", v[0], v[1], v[2], v[3]);
}

void glmPrintVector4d(const double v[4])
{
  printf("[%f,%f,%f,%f]", v[0], v[1], v[2], v[3]);
}

void glmPrintVectorf(const float *v, DWORD size)
{
  DWORD i;
  printf("[");
  for (i=0; i<size; i++) printf("%f,", v[i]);
  printf("]");
}

void glmPrintVectord(const double *v, DWORD size)
{
  DWORD i;
  printf("[");
  for (i=0; i<size; i++) printf("%f,", v[i]);
  printf("]");
}

void glmGenerateZeroVector2f(float v[2])
{
  v[0]=0;
  v[1]=0;
}

void glmGenerateZeroVector2d(double v[2])
{
  v[0]=0;
  v[1]=0;
}

void glmGenerateZeroVector3f(float v[3])
{
  v[0]=0;
  v[1]=0;
  v[2]=0;
}

void glmGenerateZeroVector3d(double v[3])
{
  v[0]=0;
  v[1]=0;
  v[2]=0;
}

void glmGenerateZeroVector4f(float *v)
{
  v[0]=0;
  v[1]=0;
  v[2]=0;
  v[3]=0;
}

void glmGenerateZeroVector4d(double *v)
{
  v[0]=0;
  v[1]=0;
  v[2]=0;
  v[3]=0;
}

void glmGenerateZeroVectorf(float *v, DWORD size)
{
  memset(v, 0, size*sizeof(float));
}

void glmGenerateZeroVectord(double *v, DWORD size)
{
  memset(v, 0, size*sizeof(double));
}

void glmColorFromVectorf(float color[3], const float vector[3])
{
  color[0]=(vector[0]+1)/2.0;
  color[1]=(vector[1]+1)/2.0;
  color[2]=(vector[2]+1)/2.0;
}

void glmColorFromVectord(double color[3], const double vector[3])
{
  color[0]=(vector[0]+1)/2.0;
  color[1]=(vector[1]+1)/2.0;
  color[2]=(vector[2]+1)/2.0;
}

void glmVectorFromColorf(float vector[3], const float color[3])
{
  vector[0]=(color[0]-0.5)*2.0;
  vector[1]=(color[1]-0.5)*2.0;
  vector[2]=(color[2]-0.5)*2.0;
}

void glmVectorFromColord(double vector[3], const double color[3])
{
  vector[0]=(color[0]-0.5)*2.0;
  vector[1]=(color[1]-0.5)*2.0;
  vector[2]=(color[2]-0.5)*2.0;
}

float glmComputeDistance2f(const float p1[2], const float p2[2])
{
  float d[2];
  glmMakeVector2f(d, p1, p2);
  return glmComputeVectorLength2f(d);
}

double glmComputeDistance2d(const double p1[2], const double p2[2])
{
  double d[2];
  glmMakeVector2d(d, p1, p2);
  return glmComputeVectorLength2d(d);
}

float glmComputeDistance3f(const float p1[3], const float p2[3])
{
  float d[3];
  glmMakeVector3f(d, p1, p2);
  return glmComputeVectorLength3f(d);
}

double glmComputeDistance3d(const double p1[3], const double p2[3])
{
  double d[3];
  glmMakeVector3d(d, p1, p2);
  return glmComputeVectorLength3d(d);
}

float glmComputeDistance4f(const float p1[4], const float p2[4])
{
  float d[4];
  glmMakeVector4f(d, p1, p2);
  return glmComputeVectorLength4f(d);
}

double glmComputeDistance4d(const double p1[4], const double p2[4])
{
  double d[4];
  glmMakeVector4d(d, p1, p2);
  return glmComputeVectorLength4d(d);
}

float glmComputeDistancef(const float *p1, const float *p2, DWORD size)
{
  float sum=0;
  float temp;
  DWORD i;
  for (i=0; i<size; i++)
  {
    temp=p1[i]-p2[i];
	sum+=temp*temp;
  }
  return sqrt(sum);
}

double glmComputeDistanced(const double *p1, const double *p2, DWORD size)
{
  double sum=0;
  double temp;
  DWORD i;
  for (i=0; i<size; i++)
  {
    temp=p1[i]-p2[i];
    sum+=temp*temp;
  }
  return sqrt(sum);
}

void glmReverseVector2f(float v[2])
{
  v[0]=-v[0];
  v[1]=-v[1];
}

void glmReverseVector2d(double v[2])
{
  v[0]=-v[0];
  v[1]=-v[1];
}

void glmReverseVector3f(float v[3])
{
  v[0]=-v[0];
  v[1]=-v[1];
  v[2]=-v[2];
}

void glmReverseVector3d(double v[3])
{
  v[0]=-v[0];
  v[1]=-v[1];
  v[2]=-v[2];
}

void glmReverseVector4f(float v[4])
{
  v[0]=-v[0];
  v[1]=-v[1];
  v[2]=-v[2];
  v[3]=-v[3];
}

void glmReverseVector4d(double v[4])
{
  v[0]=-v[0];
  v[1]=-v[1];
  v[2]=-v[2];
  v[3]=-v[3];
}

void glmReverseVectorf(float *v, DWORD size)
{
  DWORD i;
  for (i=0; i<size; i++) v[i]=-v[i];
}

void glmReverseVectord(double *v, DWORD size)
{
  DWORD i;
  for (i=0; i<size; i++) v[i]=-v[i];
}

void glmGenerateReverseVector2f(float result[2], const float v[2])
{
  result[0]=-v[0];
  result[1]=-v[1];
}

void glmGenerateReverseVector2d(double result[2], const double v[2])
{
  result[0]=-v[0];
  result[1]=-v[1];
}

void glmGenerateReverseVector3f(float result[3], const float v[3])
{
  result[0]=-v[0];
  result[1]=-v[1];
  result[2]=-v[2];
}

void glmGenerateReverseVector3d(double result[3], const double v[3])
{
  result[0]=-v[0];
  result[1]=-v[1];
  result[2]=-v[2];
}

void glmGenerateReverseVector4f(float result[4], const float v[4])
{
  result[0]=-v[0];
  result[1]=-v[1];
  result[2]=-v[2];
  result[3]=-v[3];
}

void glmGenerateReverseVector4d(double result[4], const double v[4])
{
  result[0]=-v[0];
  result[1]=-v[1];
  result[2]=-v[2];
  result[3]=-v[3];
}

void glmGenerateReverseVectorf(float *result, const float *v, DWORD size)
{
  DWORD i;
  for (i=0; i<size; i++) result[i]=-v[i];
}

void glmGenerateReverseVectord(double *result, const double *v, DWORD size)
{
  DWORD i;
  for (i=0; i<size; i++) result[i]=-v[i];
}

void glmInterpolateVector2f(float result[2], const float v1[2], const float v2[2], float a)
{
  float b=1.0-a;

  result[0]=b*v1[0]+a*v2[0];
  result[1]=b*v1[1]+a*v2[1];
}

void glmInterpolateVector2d(double result[2], const double v1[2], const double v2[2], double a)
{
  float b=1.0-a;

  result[0]=b*v1[0]+a*v2[0];
  result[1]=b*v1[1]+a*v2[1];
}

void glmInterpolateVector3f(float result[3], const float v1[3], const float v2[3], float a)
{
  float b=1.0-a;

  result[0]=b*v1[0]+a*v2[0];
  result[1]=b*v1[1]+a*v2[1];
  result[2]=b*v1[2]+a*v2[2];
}

void glmInterpolateVector3d(double result[3], const double v1[3], const double v2[3], double a)
{
  float b=1.0-a;

  result[0]=b*v1[0]+a*v2[0];
  result[1]=b*v1[1]+a*v2[1];
  result[2]=b*v1[2]+a*v2[2];
}

void glmInterpolateVector4f(float result[4], const float v1[4], const float v2[4], float a)
{
  float b=1.0-a;

  result[0]=b*v1[0]+a*v2[0];
  result[1]=b*v1[1]+a*v2[1];
  result[2]=b*v1[2]+a*v2[2];
  result[3]=b*v1[3]+a*v2[3];
}

void glmInterpolateVector4d(double result[4], const double v1[4], const double v2[4], double a)
{
  float b=1.0-a;

  result[0]=b*v1[0]+a*v2[0];
  result[1]=b*v1[1]+a*v2[1];
  result[2]=b*v1[2]+a*v2[2];
  result[3]=b*v1[3]+a*v2[3];
}

void glmInterpolateVectorf(float *result, const float *v1, const float *v2, float a, DWORD size)
{
  DWORD i;
  float b=1.0-a;

  for (i=0; i<size; i++) result[i]=b*v1[i]+a*v2[i];
}

void glmInterpolateVectord(double *result, const double *v1, const double *v2, double a, DWORD size)
{
  DWORD i;
  float b=1.0-a;

  for (i=0; i<size; i++) result[i]=b*v1[i]+a*v2[i];
}

BOOL glmCompareVector2f(const float v1[2], const float v2[2])
{
  if (v1[0]!=v2[0]) return FALSE;
  if (v1[1]!=v2[1]) return FALSE;
  return TRUE;
}

BOOL glmCompareVector2d(const double v1[2], const double v2[2])
{
  if (v1[0]!=v2[0]) return FALSE;
  if (v1[1]!=v2[1]) return FALSE;
  return TRUE;
}

BOOL glmCompareVector3f(const float v1[3], const float v2[3])
{
  if (v1[0]!=v2[0]) return FALSE;
  if (v1[1]!=v2[1]) return FALSE;
  if (v1[2]!=v2[2]) return FALSE;
  return TRUE;
}

BOOL glmCompareVector3d(const double v1[3], const double v2[3])
{
  if (v1[0]!=v2[0]) return FALSE;
  if (v1[1]!=v2[1]) return FALSE;
  if (v1[2]!=v2[2]) return FALSE;
  return TRUE;
}

BOOL glmCompareVector4f(const float v1[4], const float v2[4])
{
  if (v1[0]!=v2[0]) return FALSE;
  if (v1[1]!=v2[1]) return FALSE;
  if (v1[2]!=v2[2]) return FALSE;
  if (v1[3]!=v2[3]) return FALSE;
  return TRUE;
}

BOOL glmCompareVector4d(const double v1[4], const double v2[4])
{
  if (v1[0]!=v2[0]) return FALSE;
  if (v1[1]!=v2[1]) return FALSE;
  if (v1[2]!=v2[2]) return FALSE;
  if (v1[3]!=v2[3]) return FALSE;
  return TRUE;
}

BOOL glmCompareVectorf(const float *v1, const float *v2, DWORD size)
{
  DWORD i;
  for (i=0; i<size; i++) if (v1[i]!=v2[i]) return FALSE;
  return TRUE;
}

BOOL glmCompareVectord(const double *v1, const double *v2, DWORD size)
{
  DWORD i;
  for (i=0; i<size; i++) if (v1[i]!=v2[i]) return FALSE;
  return TRUE;
}



void glmComputeTriangleTangentsf(float sTangent[3], float tTangent[3],
                                 const float c1[3], const float c2[3], const float c3[3],
                                 const float t1[2], const float t2[2], const float t3[2])
{
  float side1[3], side2[3];
  float deltaT1, deltaT2, deltaS1, deltaS2;
  //float temp[3];

  glmMakeVector3f(side1, c1, c2);
  glmMakeVector3f(side2, c1, c3);

  deltaT1=t2[1]-t1[1];
  deltaT2=t3[1]-t1[1];
  sTangent[0]=deltaT2*side1[0]-deltaT1*side2[0];
  sTangent[1]=deltaT2*side1[1]-deltaT1*side2[1];
  sTangent[2]=deltaT2*side1[2]-deltaT1*side2[2];
  glmNormalizeVector3f(sTangent);

  deltaS1=t2[0]-t1[0];
  deltaS2=t3[0]-t1[0];
  tTangent[0]=deltaS2*side1[0]-deltaS1*side2[0];
  tTangent[1]=deltaS2*side1[1]-deltaS1*side2[1];
  tTangent[2]=deltaS2*side1[2]-deltaS1*side2[2];
  glmNormalizeVector3f(tTangent);

}


