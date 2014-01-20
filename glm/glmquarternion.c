// GLM - a math library for OpenGL
// Programmed by Markus Ilmola
// markus.ilmola@pp.inet.fi
// personal.inet.fi/koti/markus.ilmola

// Quaternion stuff
// Quaternions are presented like this: [w, (x, y, z)] w is the scalar part and (x,y,z) the vector part.
// w+xi+yj+zk

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "glm.h"

void glmAddQuaternionf(float result[4], const float q1[4], const float q2[4])
{
  result[0]=q1[0]+q2[0];
  result[1]=q1[1]+q2[1];
  result[2]=q1[2]+q2[2];
  result[3]=q1[3]+q2[3];
}

void glmAddQuaterniond(double result[4], const double q1[4], const double q2[4])
{
  result[0]=q1[0]+q2[0];
  result[1]=q1[1]+q2[1];
  result[2]=q1[2]+q2[2];
  result[3]=q1[3]+q2[3];
}

float glmComputeQuaternionNormf(const float q[4])
{
  return q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3];
}

float glmComputeQuaternionNormd(const double q[4])
{
  return q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3];
}

float glmComputeQuaternionMagnitudef(float q[4])
{
  return sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
}

double glmComputeQuaternionMagnituded(double q[4])
{
  return sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
}

float glmNormalizeQuaternionf(float q[4])
{
  float m=glmComputeQuaternionMagnitudef(q);

  if (m!=0)
  {
    q[0]/=m;
    q[1]/=m;
    q[2]/=m;
    q[3]/=m;
  }

  return m;
}

double glmNormalizeQuaterniond(double q[4])
{
  double m=glmComputeQuaternionMagnituded(q);

  if (m!=0)
  {
    q[0]/=m;
    q[1]/=m;
    q[2]/=m;
    q[3]/=m;
  }

  return m;
}

void glmMultQuaternionf(float result[4], const float q1[4], const float q2[4])
{
  result[0] = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3];
  result[1] = q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2];
  result[2] = q1[0]*q2[2] + q1[2]*q2[0] + q1[3]*q2[1] - q1[1]*q2[3];
  result[3] = q1[0]*q2[3] + q1[3]*q2[0] + q1[1]*q2[2] - q1[2]*q2[1];
}

void glmMultQuaterniond(double result[4], const double q1[4], const double q2[4])
{
  result[0] = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3];
  result[1] = q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2];
  result[2] = q1[0]*q2[2] + q1[2]*q2[0] + q1[3]*q2[1] - q1[1]*q2[3];
  result[3] = q1[0]*q2[3] + q1[3]*q2[0] + q1[1]*q2[2] - q1[2]*q2[1];
}

void glmQuaternionFromAxisAndAnglefRad(float q[4], const float axis[3], float angle)
{
  float sine = sin(angle/2.0);
  q[0]=cos(angle/2.0);
  q[1]=axis[0]*sine;
  q[2]=axis[1]*sine;
  q[3]=axis[2]*sine;
}

void glmQuaternionFromAxisAndAngledRad(double q[4], const double axis[3], double angle)
{
  double sine = sin(angle/2.0);
  q[0]=cos(angle/2.0);
  q[1]=axis[0]*sine;
  q[2]=axis[1]*sine;
  q[3]=axis[2]*sine;
}

void glmQuaternionFromAxisAndAnglefDeg(float q[4], const float axis[3], float angle)
{
  float sine;
  angle=glmRadf(angle);
  sine = sin(angle/2.0);
  q[0]=cos(angle/2.0);
  q[1]=axis[0]*sine;
  q[2]=axis[1]*sine;
  q[3]=axis[2]*sine;
}

void glmQuaternionFromAxisAndAngledDeg(double q[4], const double axis[3], double angle)
{
  double sine;
  angle=glmRadd(angle);
  sine = sin(angle/2.0);
  q[0]=cos(angle/2.0);
  q[1]=axis[0]*sine;
  q[2]=axis[1]*sine;
  q[3]=axis[2]*sine;
}

void glmGenerateConjugateQuaternionf(float result[4], const float q[4])
{
  result[0]= q[0];
  result[1]=-q[1];
  result[2]=-q[2];
  result[3]=-q[3];
}

void glmGenerateConjugateQuaterniond(double result[4], const double q[4])
{
  result[0]= q[0];
  result[1]=-q[1];
  result[2]=-q[2];
  result[3]=-q[3];
}

void glmGenerateNormalizedInverseQuaternionf(float result[4], const float q[4])
{
  result[0]= q[0];
  result[1]=-q[1];
  result[2]=-q[2];
  result[3]=-q[3];
}

void glmGenerateNormalizedInverseQuaterniond(double result[4], const double q[4])
{
  result[0]= q[0];
  result[1]=-q[1];
  result[2]=-q[2];
  result[3]=-q[3];
}

void glmGenerateQuaternionMatrix3f(float m[9], const float q[4])
{
  float xx=q[1]*q[1];
  float xy=q[1]*q[2];
  float xz=q[1]*q[3];
  float xw=q[1]*q[0];
  float yy=q[2]*q[2];
  float yz=q[2]*q[3];
  float yw=q[2]*q[0];
  float zz=q[3]*q[3];
  float zw=q[3]*q[0];

  m[0] = 1 - 2 * ( yy + zz );
  m[1] =     2 * ( xy + zw );
  m[2] =     2 * ( xz - yw );
  m[3] =     2 * ( xy - zw );
  m[4] = 1 - 2 * ( xx + zz );
  m[5] =     2 * ( yz + xw );
  m[6] =     2 * ( xz + yw );
  m[7] =     2 * ( yz - xw );
  m[8] = 1 - 2 * ( xx + yy );
}

void glmGenerateQuaternionMatrix3d(double m[9], const double q[4])
{
  double xx=q[1]*q[1];
  double xy=q[1]*q[2];
  double xz=q[1]*q[3];
  double xw=q[1]*q[0];
  double yy=q[2]*q[2];
  double yz=q[2]*q[3];
  double yw=q[2]*q[0];
  double zz=q[3]*q[3];
  double zw=q[3]*q[0];

  m[0] = 1 - 2 * ( yy + zz );
  m[1] =     2 * ( xy + zw );
  m[2] =     2 * ( xz - yw );
  m[3] =     2 * ( xy - zw );
  m[4] = 1 - 2 * ( xx + zz );
  m[5] =     2 * ( yz + xw );
  m[6] =     2 * ( xz + yw );
  m[7] =     2 * ( yz - xw );
  m[8] = 1 - 2 * ( xx + yy );
}

void glmGenerateQuaternionMatrix4f(float m[16], const float q[4])
{
  float xx=q[1]*q[1];
  float xy=q[1]*q[2];
  float xz=q[1]*q[3];
  float xw=q[1]*q[0];
  float yy=q[2]*q[2];
  float yz=q[2]*q[3];
  float yw=q[2]*q[0];
  float zz=q[3]*q[3];
  float zw=q[3]*q[0];

  m[0] = 1 - 2 * ( yy + zz );
  m[1] =     2 * ( xy + zw );
  m[2] =     2 * ( xz - yw );
  m[4] =     2 * ( xy - zw );
  m[5] = 1 - 2 * ( xx + zz );
  m[6] =     2 * ( yz + xw );
  m[8] =     2 * ( xz + yw );
  m[9] =     2 * ( yz - xw );
  m[10]= 1 - 2 * ( xx + yy );
  m[3] = m[7] = m[11] = m[12] = m[13] = m[14] = 0;
  m[15]= 1;
}

void glmGenerateQuaternionMatrix4d(double m[16], const double q[4])
{
  double xx=q[1]*q[1];
  double xy=q[1]*q[2];
  double xz=q[1]*q[3];
  double xw=q[1]*q[0];
  double yy=q[2]*q[2];
  double yz=q[2]*q[3];
  double yw=q[2]*q[0];
  double zz=q[3]*q[3];
  double zw=q[3]*q[0];

  m[0] = 1 - 2 * ( yy + zz );
  m[1] =     2 * ( xy + zw );
  m[2] =     2 * ( xz - yw );
  m[4] =     2 * ( xy - zw );
  m[5] = 1 - 2 * ( xx + zz );
  m[6] =     2 * ( yz + xw );
  m[8] =     2 * ( xz + yw );
  m[9] =     2 * ( yz - xw );
  m[10]= 1 - 2 * ( xx + yy );
  m[3] = m[7] = m[11] = m[12] = m[13] = m[14] = 0;
  m[15]= 1;
}

void glmGenerateIdentityQuaternionf(float q[4])
{
  q[0]=1;
  q[1]=0;
  q[2]=0;
  q[3]=0;
}

void glmGenerateIdentityQuaterniond(double q[4])
{
  q[0]=1;
  q[1]=0;
  q[2]=0;
  q[3]=0;
}

void glmCopyQuaternionf(const float source[4], float target[4])
{
  target[0]=source[0];
  target[1]=source[1];
  target[2]=source[2];
  target[3]=source[3];
}

void glmCopyQuaterniond(const double source[4], double target[4])
{
  target[0]=source[0];
  target[1]=source[1];
  target[2]=source[2];
  target[3]=source[3];
}

void glmQuaternionFromVectorAndScalarf(float q[4], const float v[3], const float w)
{
  q[0]=w;
  q[1]=v[0];
  q[2]=v[1];
  q[3]=v[2];
}

void glmQuaternionFromVectorAndScalard(double q[4], const double v[3], const double w)
{
  q[0]=w;
  q[1]=v[0];
  q[2]=v[1];
  q[3]=v[2];
}

void glmRotateVectorByNormalizedQuaternionf(float result[3], const float v[3], const float q[4])
{
  float temp1= v[2]*q[0]+v[1]*q[1]-v[0]*q[2];
  float temp2= v[1]*q[0]-v[2]*q[1]+v[0]*q[3];
  float temp3= v[0]*q[0]+v[2]*q[2]-v[1]*q[3];
  float temp4=-v[0]*q[1]-v[1]*q[2]-v[2]*q[3];
  result[0]= q[2]*temp1-q[3]*temp2+q[0]*temp3-q[1]*temp4;
  result[1]=-q[1]*temp1+q[0]*temp2+q[3]*temp3-q[2]*temp4;
  result[2]= q[0]*temp1+q[1]*temp2-q[2]*temp3-q[3]*temp4;
}

void glmRotateVectorByNormalizedQuaterniond(double result[3], const double v[3], const double q[4])
{
  double temp1= v[2]*q[0]+v[1]*q[1]-v[0]*q[2];
  double temp2= v[1]*q[0]-v[2]*q[1]+v[0]*q[3];
  double temp3= v[0]*q[0]+v[2]*q[2]-v[1]*q[3];
  double temp4=-v[0]*q[1]-v[1]*q[2]-v[2]*q[3];
  result[0]= q[2]*temp1-q[3]*temp2+q[0]*temp3-q[1]*temp4;
  result[1]=-q[1]*temp1+q[0]*temp2+q[3]*temp3-q[2]*temp4;
  result[2]= q[0]*temp1+q[1]*temp2-q[2]*temp3-q[3]*temp4;
}

void glmPrintQuaternionf(const float q[4])
{
  printf("%f + %f i + %f j + %f k\n", q[0], q[1], q[2], q[3]);
}

void glmPrintQuaterniond(const double q[4])
{
  printf("%f + %f i + %f j + %f k\n", q[0], q[1], q[2], q[3]);
}

void glmQuaternionSlerpf(float result[4], const float from[4], const float to[4], float t)
{
  float to1[4];
  double omega, cosom, sinom, scale0, scale1;

  // calc cosine
  cosom = from[1] * to[1] + from[2] * to[2] + from[3] * to[3] + from[0] * to[0];

  // adjust signs (if necessary)
  if (cosom <0.0 )
  {
    cosom = -cosom;
	to1[0] = - to[0];
    to1[1] = - to[1];
    to1[2] = - to[2];
    to1[3] = - to[3];
  }
  else
  {
    to1[0] = to[0];
    to1[1] = to[1];
    to1[2] = to[2];
    to1[3] = to[3];
  }

  // calculate coefficients
  if ( (1.0 - cosom) > 0.00001 )
  {
    // standard case (slerp)
    omega = acos(cosom);
    sinom = sin(omega);
    scale0 = sin((1.0 - t) * omega) / sinom;
    scale1 = sin(t * omega) / sinom;
  }
  else
  {
    // "from" and "to" quaternions are very close
    //  ... so we can do a linear interpolation
    scale0 = 1.0 - t;
    scale1 = t;
  }

  // calculate final values
  result[0] = scale0 * from[0] + scale1 * to1[0];
  result[1] = scale0 * from[1] + scale1 * to1[1];
  result[2] = scale0 * from[2] + scale1 * to1[2];
  result[3] = scale0 * from[3] + scale1 * to1[3];
}

