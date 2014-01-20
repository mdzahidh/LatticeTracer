// GLM - a math library for OpenGL
// Programmed by Markus Ilmola
// markus.ilmola@pp.inet.fi
// http://markus_ilmola.tripod.com
// http://personal.inet.fi/koti/markus.ilmolarkus.ilmola

// Matrix stuff

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "glm.h"
#include "glmvector.h"

const float  GLM_IDENTITY_MATRIX_2F[4]= {1,0, 0,1};
const double GLM_IDENTITY_MATRIX_2D[4]= {1,0, 0,1};
const float  GLM_IDENTITY_MATRIX_3F[9]= {1,0,0, 0,1,0, 0,0,1};
const double GLM_IDENTITY_MATRIX_3D[9]= {1,0,0, 0,1,0, 0,0,1};
const float  GLM_IDENTITY_MATRIX_4F[16]={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
const double GLM_IDENTITY_MATRIX_4D[16]={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};

void glmMatrixFromVectors2f(float m[4], const float v1[2], const float v2[2])
{
  m[0]=v1[0]; m[2]=v2[0];
  m[1]=v1[1]; m[3]=v2[1];
}

void glmMatrixFromVectors2d(double m[4], const double v1[2], const double v2[2])
{
  m[0]=v1[0]; m[2]=v2[0];
  m[1]=v1[1]; m[3]=v2[1];
}

void glmMatrixFromVectors2fT(float m[4], const float v1[2], const float v2[2])
{
  m[0]=v1[0]; m[2]=v1[1];
  m[1]=v2[0]; m[3]=v2[1];
}

void glmMatrixFromVectors2dT(double m[4], const double v1[2], const double v2[2])
{
  m[0]=v1[0]; m[2]=v1[1];
  m[1]=v2[0]; m[3]=v2[1];
}

void glmMatrixFromVectors3f(float m[9], const float v1[3], const float v2[3], const float v3[3])
{
  m[0]=v1[0]; m[3]=v2[0]; m[6]=v3[0];
  m[1]=v1[1]; m[4]=v2[1]; m[7]=v3[1];
  m[2]=v1[2]; m[5]=v2[2]; m[8]=v3[2];
}

void glmMatrixFromVectors3d(double m[9], const double v1[3], const double v2[3], const double v3[3])
{
  m[0]=v1[0]; m[3]=v2[0]; m[6]=v3[0];
  m[1]=v1[1]; m[4]=v2[1]; m[7]=v3[1];
  m[2]=v1[2]; m[5]=v2[2]; m[8]=v3[2];
}

void glmMatrixFromVectors3fT(float m[9], const float v1[3], const float v2[3], const float v3[3])
{
  m[0]=v1[0]; m[3]=v1[1]; m[6]=v1[2];
  m[1]=v2[0]; m[4]=v2[1]; m[7]=v2[2];
  m[2]=v3[0]; m[5]=v3[1]; m[8]=v3[2];
}

void glmMatrixFromVectors3dT(double m[9], const double v1[3], const double v2[3], const double v3[3])
{
  m[0]=v1[0]; m[3]=v1[1]; m[6]=v1[2];
  m[1]=v2[0]; m[4]=v2[1]; m[7]=v2[2];
  m[2]=v3[0]; m[5]=v3[1]; m[8]=v3[2];
}

void glmMatrixFromVectors4f(float m[16], const float v1[4], const float v2[4],
                            const float v3[4], const float v4[4])
{
  m[0]=v1[0]; m[4]=v2[0]; m[8]=v3[0];  m[12]=v4[0];
  m[1]=v1[1]; m[5]=v2[1]; m[9]=v3[1];  m[13]=v4[1];
  m[2]=v1[2]; m[6]=v2[2]; m[10]=v3[2]; m[14]=v4[2];
  m[3]=v1[3]; m[7]=v2[3]; m[11]=v3[3]; m[15]=v4[3];
}

void glmMatrixFromVectors4d(double m[16], const double v1[4], const double v2[4],
                            const double v3[4], const double v4[4])
{
  m[0]=v1[0]; m[4]=v2[0]; m[8]=v3[0];  m[12]=v4[0];
  m[1]=v1[1]; m[5]=v2[1]; m[9]=v3[1];  m[13]=v4[1];
  m[2]=v1[2]; m[6]=v2[2]; m[10]=v3[2]; m[14]=v4[2];
  m[3]=v1[3]; m[7]=v2[3]; m[11]=v3[3]; m[15]=v4[3];
}

void glmMatrixFromVectors4fT(float m[16], const float v1[4], const float v2[4],
                             const float v3[4], const float v4[4])
{
  m[0]=v1[0]; m[4]=v1[1]; m[8]=v1[2];  m[12]=v1[3];
  m[1]=v2[0]; m[5]=v2[1]; m[9]=v2[2];  m[13]=v2[3];
  m[2]=v3[0]; m[6]=v3[1]; m[10]=v3[2]; m[14]=v3[3];
  m[3]=v4[0]; m[7]=v4[1]; m[11]=v4[2]; m[15]=v4[3];
}

void glmMatrixFromVectors4dT(double m[16], const double v1[4], const double v2[4],
                             const double v3[4], const double v4[4])
{
  m[0]=v1[0]; m[4]=v1[1]; m[8]=v1[2];  m[12]=v1[3];
  m[1]=v2[0]; m[5]=v2[1]; m[9]=v2[2];  m[13]=v2[3];
  m[2]=v3[0]; m[6]=v3[1]; m[10]=v3[2]; m[14]=v3[3];
  m[3]=v4[0]; m[7]=v4[1]; m[11]=v4[2]; m[15]=v4[3];
}

void glmCopyMatrix2f(const float source[4], float target[4])
{
  memcpy(target, source, 4*sizeof(float));
}

void glmCopyMatrix2d(const double source[4], double target[4])
{
  memcpy(target, source, 4*sizeof(double));
}

void glmCopyMatrix3f(const float source[9], float target[9])
{
  memcpy(target, source, 9*sizeof(float));
}

void glmCopyMatrix3d(const double source[9], double target[9])
{
  memcpy(target, source, 9*sizeof(double));
}

void glmCopyMatrix4f(const float source[16], float target[16])
{
  memcpy(target, source, 16*sizeof(float));
}

void glmCopyMatrix4d(const double source[16], double target[16])
{
  memcpy(target, source, 16*sizeof(double));
}

void glmCopyMatrix5f(const float source[25], float target[25])
{
  memcpy(target, source, 25*sizeof(float));
}

void glmCopyMatrix5d(const double source[25], double target[25])
{
  memcpy(target, source, 25*sizeof(double));
}

// Computes the determinant of colum-mojor 2x2 matrix.
float glmDeterminant2f(const float m[4])
{
  return m[0]*m[3]-m[2]*m[1];
}

// Computes the determinant of colum-mojor 2x2 matrix.
double glmDeterminant2d(const double m[4])
{
  return m[0]*m[3]-m[2]*m[1];
}

// Computes the determinant of column-mojor 3x3 matrix.
float glmDeterminant3f(const float m[9])
{
  return
    -m[6]*m[4]*m[2]+
     m[3]*m[7]*m[2]+
     m[6]*m[1]*m[5]-
     m[0]*m[7]*m[5]-
     m[3]*m[1]*m[8]+
     m[0]*m[4]*m[8];
}

// Computes the determinant of column-mojor 3x3 matrix.
double glmDeterminant3d(const double m[9])
{
  return
    -m[6]*m[4]*m[2]+
     m[3]*m[7]*m[2]+
     m[6]*m[1]*m[5]-
     m[0]*m[7]*m[5]-
     m[3]*m[1]*m[8]+
     m[0]*m[4]*m[8];
}

// Computes the determinant of OpenGL style column-mojor 4x4 matrix.
float glmDeterminant4f(const float m[16])
{
  return
    m[12]*m[9]*m[6]*m[3]-
    m[8]*m[13]*m[6]*m[3]-
    m[12]*m[5]*m[10]*m[3]+
    m[4]*m[13]*m[10]*m[3]+
    m[8]*m[5]*m[14]*m[3]-
    m[4]*m[9]*m[14]*m[3]-
    m[12]*m[9]*m[2]*m[7]+
    m[8]*m[13]*m[2]*m[7]+
    m[12]*m[1]*m[10]*m[7]-
    m[0]*m[13]*m[10]*m[7]-
    m[8]*m[1]*m[14]*m[7]+
    m[0]*m[9]*m[14]*m[7]+
    m[12]*m[5]*m[2]*m[11]-
    m[4]*m[13]*m[2]*m[11]-
    m[12]*m[1]*m[6]*m[11]+
    m[0]*m[13]*m[6]*m[11]+
    m[4]*m[1]*m[14]*m[11]-
    m[0]*m[5]*m[14]*m[11]-
    m[8]*m[5]*m[2]*m[15]+
    m[4]*m[9]*m[2]*m[15]+
    m[8]*m[1]*m[6]*m[15]-
    m[0]*m[9]*m[6]*m[15]-
    m[4]*m[1]*m[10]*m[15]+
    m[0]*m[5]*m[10]*m[15];
}

// Computes the determinant of OpenGL style column-mojor 4x4 matrix.
double glmDeterminant4d(const double m[16])
{
  return
    m[12]*m[9]*m[6]*m[3]-
    m[8]*m[13]*m[6]*m[3]-
    m[12]*m[5]*m[10]*m[3]+
    m[4]*m[13]*m[10]*m[3]+
    m[8]*m[5]*m[14]*m[3]-
    m[4]*m[9]*m[14]*m[3]-
    m[12]*m[9]*m[2]*m[7]+
    m[8]*m[13]*m[2]*m[7]+
    m[12]*m[1]*m[10]*m[7]-
    m[0]*m[13]*m[10]*m[7]-
    m[8]*m[1]*m[14]*m[7]+
    m[0]*m[9]*m[14]*m[7]+
    m[12]*m[5]*m[2]*m[11]-
    m[4]*m[13]*m[2]*m[11]-
    m[12]*m[1]*m[6]*m[11]+
    m[0]*m[13]*m[6]*m[11]+
    m[4]*m[1]*m[14]*m[11]-
    m[0]*m[5]*m[14]*m[11]-
    m[8]*m[5]*m[2]*m[15]+
    m[4]*m[9]*m[2]*m[15]+
    m[8]*m[1]*m[6]*m[15]-
    m[0]*m[9]*m[6]*m[15]-
    m[4]*m[1]*m[10]*m[15]+
    m[0]*m[5]*m[10]*m[15];
}

// Computes the inverse matrix of a 2x2 colum-mojor matrix.
BOOL glmGenerateInverseMatrix2f(float i[4], const float m[4])
{
  float x=glmDeterminant2f(m);
  if (x==0) return FALSE;

  i[0]= m[3]/x;
  i[2]=-m[2]/x;
  i[1]=-m[1]/x;
  i[3]= m[0]/x;

  return TRUE;
}

// Computes the inverse matrix of a 2x2 colum-mojor matrix.
BOOL glmGenerateInverseMatrix2d(double i[4], const double m[4])
{
  double x=glmDeterminant2d(m);
  if (x==0) return FALSE;

  i[0]= m[3]/x;
  i[2]=-m[2]/x;
  i[1]=-m[1]/x;
  i[3]= m[0]/x;

  return TRUE;
}

// Computes the inverse matrix of a 3x3 colum-mojor matrix.
BOOL glmGenerateInverseMatrix3f(float i[9], const float m[9])
{
  float x=glmDeterminant3f(m);
  if (x==0) return FALSE;

  i[0]=(-m[7]*m[5]+m[4]*m[8])/x;
  i[3]=( m[6]*m[5]-m[3]*m[8])/x;
  i[6]=(-m[6]*m[4]+m[3]*m[7])/x;
  i[1]=( m[7]*m[2]-m[1]*m[8])/x;
  i[4]=(-m[6]*m[2]+m[0]*m[8])/x;
  i[7]=( m[6]*m[1]-m[0]*m[7])/x;
  i[2]=(-m[4]*m[2]+m[1]*m[5])/x;
  i[5]=( m[3]*m[2]-m[0]*m[5])/x;
  i[8]=(-m[3]*m[1]+m[0]*m[4])/x;

  return TRUE;
}

// Computes the inverse matrix of a 3x3 column-mojor matrix.
BOOL glmGenerateInverseMatrix3d(double i[9], const double m[9])
{
  double x=glmDeterminant3d(m);
  if (x==0) return FALSE;

  i[0]=(-m[7]*m[5]+m[4]*m[8])/x;
  i[3]=( m[6]*m[5]-m[3]*m[8])/x;
  i[6]=(-m[6]*m[4]+m[3]*m[7])/x;
  i[1]=( m[7]*m[2]-m[1]*m[8])/x;
  i[4]=(-m[6]*m[2]+m[0]*m[8])/x;
  i[7]=( m[6]*m[1]-m[0]*m[7])/x;
  i[2]=(-m[4]*m[2]+m[1]*m[5])/x;
  i[5]=( m[3]*m[2]-m[0]*m[5])/x;
  i[8]=(-m[3]*m[1]+m[0]*m[4])/x;

  return TRUE;
}

// Computes the inverse matrix of a OpenGL style 4x4 column-mojor matrix.
BOOL glmGenerateInverseMatrix4f(float i[16], const float m[16])
{
  float x=glmDeterminant4f(m);
  if (x==0) return FALSE;

  i[0]= (-m[13]*m[10]*m[7] +m[9]*m[14]*m[7] +m[13]*m[6]*m[11]
         -m[5]*m[14]*m[11] -m[9]*m[6]*m[15] +m[5]*m[10]*m[15])/x;
  i[4]= ( m[12]*m[10]*m[7] -m[8]*m[14]*m[7] -m[12]*m[6]*m[11]
         +m[4]*m[14]*m[11] +m[8]*m[6]*m[15] -m[4]*m[10]*m[15])/x;
  i[8]= (-m[12]*m[9]* m[7] +m[8]*m[13]*m[7] +m[12]*m[5]*m[11]
         -m[4]*m[13]*m[11] -m[8]*m[5]*m[15] +m[4]*m[9]* m[15])/x;
  i[12]=( m[12]*m[9]* m[6] -m[8]*m[13]*m[6] -m[12]*m[5]*m[10]
         +m[4]*m[13]*m[10] +m[8]*m[5]*m[14] -m[4]*m[9]* m[14])/x;
  i[1]= ( m[13]*m[10]*m[3] -m[9]*m[14]*m[3] -m[13]*m[2]*m[11]
         +m[1]*m[14]*m[11] +m[9]*m[2]*m[15] -m[1]*m[10]*m[15])/x;
  i[5]= (-m[12]*m[10]*m[3] +m[8]*m[14]*m[3] +m[12]*m[2]*m[11]
         -m[0]*m[14]*m[11] -m[8]*m[2]*m[15] +m[0]*m[10]*m[15])/x;
  i[9]= ( m[12]*m[9]* m[3] -m[8]*m[13]*m[3] -m[12]*m[1]*m[11]
         +m[0]*m[13]*m[11] +m[8]*m[1]*m[15] -m[0]*m[9]* m[15])/x;
  i[13]=(-m[12]*m[9]* m[2] +m[8]*m[13]*m[2] +m[12]*m[1]*m[10]
         -m[0]*m[13]*m[10] -m[8]*m[1]*m[14] +m[0]*m[9]* m[14])/x;
  i[2]= (-m[13]*m[6]* m[3] +m[5]*m[14]*m[3] +m[13]*m[2]*m[7]
         -m[1]*m[14]*m[7]  -m[5]*m[2]*m[15] +m[1]*m[6]* m[15])/x;
  i[6]= ( m[12]*m[6]* m[3] -m[4]*m[14]*m[3] -m[12]*m[2]*m[7]
         +m[0]*m[14]*m[7]  +m[4]*m[2]*m[15] -m[0]*m[6]* m[15])/x;
  i[10]=(-m[12]*m[5]* m[3] +m[4]*m[13]*m[3] +m[12]*m[1]*m[7]
         -m[0]*m[13]*m[7]  -m[4]*m[1]*m[15] +m[0]*m[5]* m[15])/x;
  i[14]=( m[12]*m[5]* m[2] -m[4]*m[13]*m[2] -m[12]*m[1]*m[6]
         +m[0]*m[13]*m[6]  +m[4]*m[1]*m[14] -m[0]*m[5]* m[14])/x;
  i[3]= ( m[9]* m[6]* m[3] -m[5]*m[10]*m[3] -m[9]* m[2]*m[7]
         +m[1]*m[10]*m[7]  +m[5]*m[2]*m[11] -m[1]*m[6]* m[11])/x;
  i[7]= (-m[8]* m[6]* m[3] +m[4]*m[10]*m[3] +m[8]* m[2]*m[7]
         -m[0]*m[10]*m[7]  -m[4]*m[2]*m[11] +m[0]*m[6]* m[11])/x;
  i[11]=( m[8]* m[5]* m[3] -m[4]*m[9]* m[3] -m[8]* m[1]*m[7]
         +m[0]*m[9]* m[7]  +m[4]*m[1]*m[11] -m[0]*m[5]* m[11])/x;
  i[15]=(-m[8]* m[5]* m[2] +m[4]*m[9]* m[2] +m[8]* m[1]*m[6]
         -m[0]*m[9]* m[6]  -m[4]*m[1]*m[10] +m[0]*m[5]* m[10])/x;

  return TRUE;
}

// Computes the inverse matrix of a OpenGL style 4x4 colum-mojor matrix.
BOOL glmGenerateInverseMatrix4d(double i[16], const double m[16])
{
  double x=glmDeterminant4d(m);
  if (x==0) return FALSE;

  i[0]= (-m[13]*m[10]*m[7] +m[9]*m[14]*m[7] +m[13]*m[6]*m[11]
         -m[5]*m[14]*m[11] -m[9]*m[6]*m[15] +m[5]*m[10]*m[15])/x;
  i[4]= ( m[12]*m[10]*m[7] -m[8]*m[14]*m[7] -m[12]*m[6]*m[11]
         +m[4]*m[14]*m[11] +m[8]*m[6]*m[15] -m[4]*m[10]*m[15])/x;
  i[8]= (-m[12]*m[9]* m[7] +m[8]*m[13]*m[7] +m[12]*m[5]*m[11]
         -m[4]*m[13]*m[11] -m[8]*m[5]*m[15] +m[4]*m[9]* m[15])/x;
  i[12]=( m[12]*m[9]* m[6] -m[8]*m[13]*m[6] -m[12]*m[5]*m[10]
         +m[4]*m[13]*m[10] +m[8]*m[5]*m[14] -m[4]*m[9]* m[14])/x;
  i[1]= ( m[13]*m[10]*m[3] -m[9]*m[14]*m[3] -m[13]*m[2]*m[11]
         +m[1]*m[14]*m[11] +m[9]*m[2]*m[15] -m[1]*m[10]*m[15])/x;
  i[5]= (-m[12]*m[10]*m[3] +m[8]*m[14]*m[3] +m[12]*m[2]*m[11]
         -m[0]*m[14]*m[11] -m[8]*m[2]*m[15] +m[0]*m[10]*m[15])/x;
  i[9]= ( m[12]*m[9]* m[3] -m[8]*m[13]*m[3] -m[12]*m[1]*m[11]
         +m[0]*m[13]*m[11] +m[8]*m[1]*m[15] -m[0]*m[9]* m[15])/x;
  i[13]=(-m[12]*m[9]* m[2] +m[8]*m[13]*m[2] +m[12]*m[1]*m[10]
         -m[0]*m[13]*m[10] -m[8]*m[1]*m[14] +m[0]*m[9]* m[14])/x;
  i[2]= (-m[13]*m[6]* m[3] +m[5]*m[14]*m[3] +m[13]*m[2]*m[7]
         -m[1]*m[14]*m[7]  -m[5]*m[2]*m[15] +m[1]*m[6]* m[15])/x;
  i[6]= ( m[12]*m[6]* m[3] -m[4]*m[14]*m[3] -m[12]*m[2]*m[7]
         +m[0]*m[14]*m[7]  +m[4]*m[2]*m[15] -m[0]*m[6]* m[15])/x;
  i[10]=(-m[12]*m[5]* m[3] +m[4]*m[13]*m[3] +m[12]*m[1]*m[7]
         -m[0]*m[13]*m[7]  -m[4]*m[1]*m[15] +m[0]*m[5]* m[15])/x;
  i[14]=( m[12]*m[5]* m[2] -m[4]*m[13]*m[2] -m[12]*m[1]*m[6]
         +m[0]*m[13]*m[6]  +m[4]*m[1]*m[14] -m[0]*m[5]* m[14])/x;
  i[3]= ( m[9]* m[6]* m[3] -m[5]*m[10]*m[3] -m[9]* m[2]*m[7]
         +m[1]*m[10]*m[7]  +m[5]*m[2]*m[11] -m[1]*m[6]* m[11])/x;
  i[7]= (-m[8]* m[6]* m[3] +m[4]*m[10]*m[3] +m[8]* m[2]*m[7]
         -m[0]*m[10]*m[7]  -m[4]*m[2]*m[11] +m[0]*m[6]* m[11])/x;
  i[11]=( m[8]* m[5]* m[3] -m[4]*m[9]* m[3] -m[8]* m[1]*m[7]
         +m[0]*m[9]* m[7]  +m[4]*m[1]*m[11] -m[0]*m[5]* m[11])/x;
  i[15]=(-m[8]* m[5]* m[2] +m[4]*m[9]* m[2] +m[8]* m[1]*m[6]
         -m[0]*m[9]* m[6]  -m[4]*m[1]*m[10] +m[0]*m[5]* m[10])/x;

  return TRUE;
}

void glmGenerateTransposeMatrix2f(float t[4], const float m[4])
{
  t[0]=m[0]; t[2]=m[1];
  t[1]=m[2]; t[3]=m[3];
}

void glmGenerateTransposeMatrix2d(double t[4], const double m[4])
{
  t[0]=m[0]; t[2]=m[1];
  t[1]=m[2]; t[3]=m[3];
}

void glmGenerateTransposeMatrix3f(float t[9], const float m[9])
{
  t[0]=m[0]; t[3]=m[1]; t[6]=m[2];
  t[1]=m[3]; t[4]=m[4]; t[7]=m[5];
  t[2]=m[6]; t[5]=m[7]; t[8]=m[8];
}

void glmGenerateTransposeMatrix3d(double t[9], const double m[9])
{
  t[0]=m[0]; t[3]=m[1]; t[6]=m[2];
  t[1]=m[3]; t[4]=m[4]; t[7]=m[5];
  t[2]=m[6]; t[5]=m[7]; t[8]=m[8];
}

void glmGenerateTransposeMatrix4f(float t[16], const float m[16])
{
  t[0]=m[0];  t[4]=m[1];  t[8]=m[2];   t[12]=m[3];
  t[1]=m[4];  t[5]=m[5];  t[9]=m[6];   t[13]=m[7];
  t[2]=m[8];  t[6]=m[9];  t[10]=m[10]; t[14]=m[11];
  t[3]=m[12]; t[7]=m[13]; t[11]=m[14]; t[15]=m[15];
}

void glmGenerateTransposeMatrix4d(double t[16], const double m[16])
{
  t[0]=m[0];  t[4]=m[1];  t[8]=m[2];   t[12]=m[3];
  t[1]=m[4];  t[5]=m[5];  t[9]=m[6];   t[13]=m[7];
  t[2]=m[8];  t[6]=m[9];  t[10]=m[10]; t[14]=m[11];
  t[3]=m[12]; t[7]=m[13]; t[11]=m[14]; t[15]=m[15];
}

void glmTransposeMatrix2f(float m[4])
{
  float temp;
  temp=m[1];
  m[1]=m[2];
  m[2]=temp;
}

void glmTransposeMatrix2d(double m[4])
{
  double temp;
  temp=m[1];
  m[1]=m[2];
  m[2]=temp;
}

void glmTransposeMatrix3f(float m[9])
{
  float temp;
  DWORD i, j;
  for (i=0; i<2; i++)
  {
    for (j=i+1; j<3; j++)
    {
      temp=m[j+3*i];
      m[j+3*i]=m[3*j+i];
      m[3*j+i]=temp;
    }
  }
}

void glmTransposeMatrix3d(double m[9])
{
  double temp;
  DWORD i, j;
  for (i=0; i<2; i++)
  {
    for (j=i+1; j<3; j++)
    {
      temp=m[j+3*i];
      m[j+3*i]=m[3*j+i];
      m[3*j+i]=temp;
    }
  }
}

void glmTransposeMatrix4f(float m[16])
{
  float temp;
  DWORD i, j;
  for (i=0; i<3; i++)
  {
    for (j=i+1; j<4; j++)
    {
      temp=m[j+4*i];
      m[j+4*i]=m[4*j+i];
      m[4*j+i]=temp;
    }
  }
}

void glmTransposeMatrix4d(double m[16])
{
  double temp;
  DWORD i, j;
  for (i=0; i<3; i++)
  {
    for (j=i+1; j<4; j++)
    {
      temp=m[j+4*i];
      m[j+4*i]=m[4*j+i];
      m[4*j+i]=temp;
    }
  }
}

void glmTransposeMatrixf(float *m, DWORD size)
{
  float temp;
  DWORD i, j;
  for (i=0; i<size-1; i++)
  {
    for (j=i+1; j<size; j++)
    {
      temp=m[j+size*i];
      m[j+size*i]=m[size*j+i];
      m[size*j+i]=temp;
    }
  }
}

void glmTransposeMatrixd(double *m, DWORD size)
{
  double temp;
  DWORD i, j;
  for (i=0; i<size-1; i++)
  {
    for (j=i+1; j<size; j++)
    {
      temp=m[j+size*i];
      m[j+size*i]=m[size*j+i];
      m[size*j+i]=temp;
    }
  }
}

// Computes m1*m2
void glmMultMatrix2f(float result[4], const float m1[4], const float m2[4])
{
  float sum;
  int x, y, i;
  for (x=0; x<2; x++)
  {
    for (y=0; y<2; y++)
    {
      sum=0;
      for (i=0; i<2; i++) sum+=m1[y+2*i]*m2[2*x+i];
      result[y+2*x]=sum;
    }
  }
}

// Computes m1*m2
void glmMultMatrix2d(double result[4], const double m1[4], const double m2[4])
{
  double sum;
  int x, y, i;
  for (x=0; x<2; x++)
  {
    for (y=0; y<2; y++)
    {
      sum=0;
      for (i=0; i<2; i++) sum+=m1[y+2*i]*m2[2*x+i];
      result[y+2*x]=sum;
    }
  }
}

// Computes m1*m2
void glmMultMatrix3f(float result[9], const float m1[9], const float m2[9])
{
  float sum;
  int x, y, i;
  for (x=0; x<3; x++)
  {
    for (y=0; y<3; y++)
    {
      sum=0;
      for (i=0; i<3; i++) sum+=m1[y+3*i]*m2[3*x+i];
      result[y+3*x]=sum;
    }
  }
}

// Computes m1*m2
void glmMultMatrix3d(double result[9], const double m1[9], const double m2[9])
{
  double sum;
  int x, y, i;
  for (x=0; x<3; x++)
  {
    for (y=0; y<3; y++)
    {
      sum=0;
      for (i=0; i<3; i++) sum+=m1[y+3*i]*m2[3*x+i];
      result[y+3*x]=sum;
    }
  }
}

// Computes m1*m2
void glmMultMatrix4f(float result[16], const float m1[16], const float m2[16])
{
  float sum;
  int x, y, i;
  for (x=0; x<4; x++)
  {
    for (y=0; y<4; y++)
    {
      sum=0;
      for (i=0; i<4; i++) sum+=m1[y+4*i]*m2[4*x+i];
      result[y+4*x]=sum;
    }
  }
}

void glmMultMatrix4d(double result[16], const double m1[16], const double m2[16])
{
  double sum;
  int x, y, i;
  for (x=0; x<4; x++)
  {
    for (y=0; y<4; y++)
    {
      sum=0;
      for (i=0; i<4; i++) sum+=m1[y+4*i]*m2[4*x+i];
      result[y+4*x]=sum;
    }
  }
}

void glmMultMatrix5f(float result[25], const float m1[25], const float m2[25])
{
  float sum;
  int x, y, i;
  for (x=0; x<5; x++)
  {
    for (y=0; y<5; y++)
    {
      sum=0;
      for (i=0; i<5; i++) sum+=m1[y+5*i]*m2[5*x+i];
      result[y+5*x]=sum;
    }
  }
}

void glmMultMatrix5d(double result[25], const double m1[25], const double m2[25])
{
  double sum;
  int x, y, i;
  for (x=0; x<5; x++)
  {
    for (y=0; y<5; y++)
    {
      sum=0;
      for (i=0; i<5; i++) sum+=m1[y+5*i]*m2[5*x+i];
      result[y+5*x]=sum;
    }
  }
}

void glmAddMatrix2f(float result[4], const float m1[4], const float m2[4])
{
  int i;
  for (i=0; i<4; i++)
    result[i]=m1[i]+m2[i];
}

void glmAddMatrix2d(double result[4], const double m1[4], const double m2[4])
{
  int i;
  for (i=0; i<4; i++)
    result[i]=m1[i]+m2[i];
}

void glmAddMatrix3f(float result[9], const float m1[9], const float m2[9])
{
  int i;
  for (i=0; i<9; i++)
    result[i]=m1[i]+m2[i];
}

void glmAddMatrix3d(double result[9], const double m1[9], const double m2[9])
{
  int i;
  for (i=0; i<9; i++)
    result[i]=m1[i]+m2[i];
}

void glmAddMatrix4f(float result[16], const float m1[16], const float m2[16])
{
  int i;
  for (i=0; i<16; i++)
    result[i]=m1[i]+m2[i];
}

void glmAddMatrix4d(double result[16], const double m1[16], const double m2[16])
{
  int i;
  for (i=0; i<16; i++)
    result[i]=m1[i]+m2[i];
}

// Computes m1-m2
void glmSubtractMatrix2f(float result[4], const float m1[4], const float m2[4])
{
  int i;
  for (i=0; i<4; i++)
    result[i]=m1[i]-m2[i];
}

// Computes m1-m2
void glmSubtractMatrix2d(double result[4], const double m1[4], const double m2[4])
{
  int i;
  for (i=0; i<4; i++)
    result[i]=m1[i]-m2[i];
}

// Computes m1-m2
void glmSubtractMatrix3f(float result[9], const float m1[9], const float m2[9])
{
  int i;
  for (i=0; i<9; i++)
    result[i]=m1[i]-m2[i];
}

// Computes m1-m2
void glmSubtractMatrix3d(double result[9], const double m1[9], const double m2[9])
{
  int i;
  for (i=0; i<9; i++)
    result[i]=m1[i]-m2[i];
}

// Computes m1-m2
void glmSubtractMatrix4f(float result[16], const float m1[16], const float m2[16])
{
  int i;
  for (i=0; i<16; i++)
    result[i]=m1[i]-m2[i];
}

// Computes m1-m2
void glmSubtractMatrix4d(double result[16], const double m1[16], const double m2[16])
{
  int i;
  for (i=0; i<16; i++)
    result[i]=m1[i]-m2[i];
}

void glmGenerateRotationMatrix2fRad(float m[4], float angle)
{
  // | cos(A)  -sin(A) |
  // | sin(A)   cos(A) |
  m[0]=cos(angle);  m[2]=-sin(angle);
  m[1]=sin(angle);	m[3]=cos(angle);
}

void glmGenerateRotationMatrix2fDeg(float m[4], float angle)
{
  angle=glmRadf(angle);
  // | cos(A)  -sin(A) |
  // | sin(A)   cos(A) |
  m[0]=cos(angle);  m[2]=-sin(angle);
  m[1]=sin(angle);	m[3]=cos(angle);
}

void glmGenerateRotationMatrix2dRad(double m[4], double angle)
{
  // | cos(A)  -sin(A) |
  // | sin(A)   cos(A) |
  m[0]=cos(angle);  m[2]=-sin(angle);
  m[1]=sin(angle);	m[3]=cos(angle);
}

void glmGenerateRotationMatrix2dDeg(double m[4], double angle)
{
  angle=glmRadd(angle);
  // | cos(A)  -sin(A) |
  // | sin(A)   cos(A) |
  m[0]=cos(angle);  m[2]=-sin(angle);
  m[1]=sin(angle);	m[3]=cos(angle);
}

void glmGenerateXRotationMatrix3fRad(float m[16], float angle)
{
  //  |  1  0       0     |
  //M=|  0  cos(A) -sin(A)|
  //  |  0  sin(A)  cos(A)|
  m[0]=1; m[3]=0;          m[6]=0;
  m[1]=0; m[4]=cos(angle); m[7]=-sin(angle);
  m[2]=0; m[5]=sin(angle); m[8]=cos(angle);
}

void glmGenerateXRotationMatrix3fDeg(float m[9], float angle)
{
  angle=glmRadf(angle);
  //  |  1  0       0     |
  //M=|  0  cos(A) -sin(A)|
  //  |  0  sin(A)  cos(A)|
  m[0]=1; m[3]=0;          m[6]=0;
  m[1]=0; m[4]=cos(angle); m[7]=-sin(angle);
  m[2]=0; m[5]=sin(angle); m[8]=cos(angle);
}

void glmGenerateXRotationMatrix3dRad(double m[9], double angle)
{
  //  |  1  0       0     |
  //M=|  0  cos(A) -sin(A)|
  //  |  0  sin(A)  cos(A)|
  m[0]=1; m[3]=0;          m[6]=0;
  m[1]=0; m[4]=cos(angle); m[7]=-sin(angle);
  m[2]=0; m[5]=sin(angle); m[8]=cos(angle);
}

void glmGenerateXRotationMatrix3dDeg(double m[9], double angle)
{
  angle=glmRadd(angle);
  //  |  1  0       0     |
  //M=|  0  cos(A) -sin(A)|
  //  |  0  sin(A)  cos(A)|
  m[0]=1; m[3]=0;          m[6]=0;
  m[1]=0; m[4]=cos(angle); m[7]=-sin(angle);
  m[2]=0; m[5]=sin(angle); m[8]=cos(angle);
}

void glmGenerateXRotationMatrix4fRad(float m[16], float angle)
{
  //  |  1  0       0       0 |
  //M=|  0  cos(A) -sin(A)  0 |
  //  |  0  sin(A)  cos(A)  0 |
  //  |  0  0       0       1 |
  m[0]=1; m[4]=0;          m[8]=0;           m[12]=0;
  m[1]=0; m[5]=cos(angle); m[9]=-sin(angle); m[13]=0;
  m[2]=0; m[6]=sin(angle); m[10]=cos(angle); m[14]=0;
  m[3]=0; m[7]=0;          m[11]=0;          m[15]=1;
}

void glmGenerateXRotationMatrix4fDeg(float m[16], float angle)
{
  angle=glmRadf(angle);
  //  |  1  0       0       0 |
  //M=|  0  cos(A) -sin(A)  0 |
  //  |  0  sin(A)  cos(A)  0 |
  //  |  0  0       0       1 |
  m[0]=1; m[4]=0;          m[8]=0;           m[12]=0;
  m[1]=0; m[5]=cos(angle); m[9]=-sin(angle); m[13]=0;
  m[2]=0; m[6]=sin(angle); m[10]=cos(angle); m[14]=0;
  m[3]=0; m[7]=0;          m[11]=0;          m[15]=1;
}

void glmGenerateXRotationMatrix4dRad(double m[16], double angle)
{
  //  |  1  0       0       0 |
  //M=|  0  cos(A) -sin(A)  0 |
  //  |  0  sin(A)  cos(A)  0 |
  //  |  0  0       0       1 |
  m[0]=1; m[4]=0;          m[8]=0;           m[12]=0;
  m[1]=0; m[5]=cos(angle); m[9]=-sin(angle); m[13]=0;
  m[2]=0; m[6]=sin(angle); m[10]=cos(angle); m[14]=0;
  m[3]=0; m[7]=0;          m[11]=0;          m[15]=1;
}

void glmGenerateXRotationMatrix4dDeg(double m[16], double angle)
{
  angle=glmRadd(angle);
  //  |  1  0       0       0 |
  //M=|  0  cos(A) -sin(A)  0 |
  //  |  0  sin(A)  cos(A)  0 |
  //  |  0  0       0       1 |
  m[0]=1; m[4]=0;          m[8]=0;           m[12]=0;
  m[1]=0; m[5]=cos(angle); m[9]=-sin(angle); m[13]=0;
  m[2]=0; m[6]=sin(angle); m[10]=cos(angle); m[14]=0;
  m[3]=0; m[7]=0;          m[11]=0;          m[15]=1;
}

void glmGenerateYRotationMatrix3fRad(float m[9], float angle)
{
  //  |  cos(A)  0  sin(A) |
  //M=|  0       1  0      |
  //  | -sin(A)  0  cos(A) |
  m[0]=cos(angle);  m[3]=0; m[6]=sin(angle);
  m[1]=0;           m[4]=1; m[7]=0;
  m[2]=-sin(angle); m[5]=0; m[8]=cos(angle);
}

void glmGenerateYRotationMatrix3fDeg(float m[9], float angle)
{
  angle=glmRadf(angle);
  //  |  cos(A)  0  sin(A) |
  //M=|  0       1  0      |
  //  | -sin(A)  0  cos(A) |
  m[0]=cos(angle);  m[3]=0; m[6]=sin(angle);
  m[1]=0;           m[4]=1; m[7]=0;
  m[2]=-sin(angle); m[5]=0; m[8]=cos(angle);
}

void glmGenerateYRotationMatrix3dRad(double m[9], double angle)
{
  //  |  cos(A)  0  sin(A) |
  //M=|  0       1  0      |
  //  | -sin(A)  0  cos(A) |
  m[0]=cos(angle);  m[3]=0; m[6]=sin(angle);
  m[1]=0;           m[4]=1; m[7]=0;
  m[2]=-sin(angle); m[5]=0; m[8]=cos(angle);
}

void glmGenerateYRotationMatrix3dDeg(double m[9], double angle)
{
  angle=glmRadd(angle);
  //  |  cos(A)  0  sin(A) |
  //M=|  0       1  0      |
  //  | -sin(A)  0  cos(A) |
  m[0]=cos(angle);  m[3]=0; m[6]=sin(angle);
  m[1]=0;           m[4]=1; m[7]=0;
  m[2]=-sin(angle); m[5]=0; m[8]=cos(angle);
}

void glmGenerateYRotationMatrix4fRad(float m[16], float angle)
{
  //  |  cos(A)  0  sin(A)  0 |
  //M=|  0       1  0       0 |
  //  | -sin(A)  0  cos(A)  0 |
  //  |  0       0  0       1 |
  m[0]=cos(angle);  m[4]=0; m[8]=sin(angle);  m[12]=0;
  m[1]=0;           m[5]=1; m[9]=0;           m[13]=0;
  m[2]=-sin(angle); m[6]=0; m[10]=cos(angle); m[14]=0;
  m[3]=0;           m[7]=0; m[11]=0;          m[15]=1;
}

void glmGenerateYRotationMatrix4fDeg(float m[16], float angle)
{
  angle=glmRadf(angle);
  //  |  cos(A)  0  sin(A)  0 |
  //M=|  0       1  0       0 |
  //  | -sin(A)  0  cos(A)  0 |
  //  |  0       0  0       1 |
  m[0]=cos(angle);  m[4]=0; m[8]=sin(angle);  m[12]=0;
  m[1]=0;           m[5]=1; m[9]=0;           m[13]=0;
  m[2]=-sin(angle); m[6]=0; m[10]=cos(angle); m[14]=0;
  m[3]=0;           m[7]=0; m[11]=0;          m[15]=1;
}

void glmGenerateYRotationMatrix4dRad(double m[16], double angle)
{
  //  |  cos(A)  0  sin(A)  0 |
  //M=|  0       1  0       0 |
  //  | -sin(A)  0  cos(A)  0 |
  //  |  0       0  0       1 |
  m[0]=cos(angle);  m[4]=0; m[8]=sin(angle);  m[12]=0;
  m[1]=0;           m[5]=1; m[9]=0;           m[13]=0;
  m[2]=-sin(angle); m[6]=0; m[10]=cos(angle); m[14]=0;
  m[3]=0;           m[7]=0; m[11]=0;          m[15]=1;
}

void glmGenerateYRotationMatrix4dDeg(double m[16], double angle)
{
  angle=glmRadd(angle);
  //  |  cos(A)  0  sin(A)  0 |
  //M=|  0       1  0       0 |
  //  | -sin(A)  0  cos(A)  0 |
  //  |  0       0  0       1 |
  m[0]=cos(angle);  m[4]=0; m[8]=sin(angle);  m[12]=0;
  m[1]=0;           m[5]=1; m[9]=0;           m[13]=0;
  m[2]=-sin(angle); m[6]=0; m[10]=cos(angle); m[14]=0;
  m[3]=0;           m[7]=0; m[11]=0;          m[15]=1;
}

void glmGenerateZRotationMatrix3fRad(float m[9], float angle)
{
  //  |  cos(A)  -sin(A) 0 |
  //M=|  sin(A)   cos(A) 0 |
  //  |  0        0      1 |
  m[0]=cos(angle); m[3]=-sin(angle); m[6]=0;
  m[1]=sin(angle); m[4]=cos(angle);  m[7]=0;
  m[2]=0;          m[5]=0;           m[8]=1;
}

void glmGenerateZRotationMatrix3fDeg(float m[9], float angle)
{
  angle=glmRadf(angle);
  //  |  cos(A)  -sin(A) 0 |
  //M=|  sin(A)   cos(A) 0 |
  //  |  0        0      1 |
  m[0]=cos(angle); m[3]=-sin(angle); m[6]=0;
  m[1]=sin(angle); m[4]=cos(angle);  m[7]=0;
  m[2]=0;          m[5]=0;           m[8]=1;
}

void glmGenerateZRotationMatrix3dRad(double m[9], double angle)
{
  //  |  cos(A)  -sin(A) 0 |
  //M=|  sin(A)   cos(A) 0 |
  //  |  0        0      1 |
  m[0]=cos(angle); m[3]=-sin(angle); m[6]=0;
  m[1]=sin(angle); m[4]=cos(angle);  m[7]=0;
  m[2]=0;          m[5]=0;           m[8]=1;
}

void glmGenerateZRotationMatrix3dDeg(double m[9], double angle)
{
  angle=glmRadd(angle);
  //  |  cos(A)  -sin(A) 0 |
  //M=|  sin(A)   cos(A) 0 |
  //  |  0        0      1 |
  m[0]=cos(angle); m[3]=-sin(angle); m[6]=0;
  m[1]=sin(angle); m[4]=cos(angle);  m[7]=0;
  m[2]=0;          m[5]=0;           m[8]=1;
}

void glmGenerateZRotationMatrix4fRad(float m[16], float angle)
{
  //  |  cos(A)  -sin(A) 0  0 |
  //M=|  sin(A)   cos(A) 0  0 |
  //  |  0        0      1  0 |
  //  |  0        0      0  1 |
  m[0]=cos(angle); m[4]=-sin(angle); m[8]=0;  m[12]=0;
  m[1]=sin(angle); m[5]=cos(angle);  m[9]=0;  m[13]=0;
  m[2]=0;          m[6]=0;           m[10]=1; m[14]=0;
  m[3]=0;          m[7]=0;           m[11]=0; m[15]=1;
}

void glmGenerateZRotationMatrix4fDeg(float m[16], float angle)
{
  angle=glmRadf(angle);
  //  |  cos(A)  -sin(A) 0  0 |
  //M=|  sin(A)   cos(A) 0  0 |
  //  |  0        0      1  0 |
  //  |  0        0      0  1 |
  m[0]=cos(angle); m[4]=-sin(angle); m[8]=0;  m[12]=0;
  m[1]=sin(angle); m[5]=cos(angle);  m[9]=0;  m[13]=0;
  m[2]=0;          m[6]=0;           m[10]=1; m[14]=0;
  m[3]=0;          m[7]=0;           m[11]=0; m[15]=1;
}

void glmGenerateZRotationMatrix4dRad(double m[16], double angle)
{
  //  |  cos(A)  -sin(A) 0  0 |
  //M=|  sin(A)   cos(A) 0  0 |
  //  |  0        0      1  0 |
  //  |  0        0      0  1 |
  m[0]=cos(angle); m[4]=-sin(angle); m[8]=0;  m[12]=0;
  m[1]=sin(angle); m[5]=cos(angle);  m[9]=0;  m[13]=0;
  m[2]=0;          m[6]=0;           m[10]=1; m[14]=0;
  m[3]=0;          m[7]=0;           m[11]=0; m[15]=1;
}

void glmGenerateZRotationMatrix4dDeg(double m[16], double angle)
{
  angle=glmRadd(angle);
  //  |  cos(A)  -sin(A) 0  0 |
  //M=|  sin(A)   cos(A) 0  0 |
  //  |  0        0      1  0 |
  //  |  0        0      0  1 |
  m[0]=cos(angle); m[4]=-sin(angle); m[8]=0;  m[12]=0;
  m[1]=sin(angle); m[5]=cos(angle);  m[9]=0;  m[13]=0;
  m[2]=0;          m[6]=0;           m[10]=1; m[14]=0;
  m[3]=0;          m[7]=0;           m[11]=0; m[15]=1;
}

void glmGenerateXYRotationMatrix4fDeg(float m[16], float angle)
{
  angle=glmRadf(angle);
  m[0]=1; m[4]=0; m[8 ]=0;          m[12]=0;
  m[1]=0; m[5]=1; m[9 ]=0;          m[13]=0;
  m[2]=0; m[6]=0; m[10]=cos(angle); m[14]=-sin(angle);
  m[3]=0; m[7]=0; m[11]=sin(angle); m[15]=cos(angle);
}

void glmGenerateXYRotationMatrix4dDeg(double m[16], double angle)
{
  angle=glmRadd(angle);
  m[0]=1; m[4]=0; m[8 ]=0;          m[12]=0;
  m[1]=0; m[5]=1; m[9 ]=0;          m[13]=0;
  m[2]=0; m[6]=0; m[10]=cos(angle); m[14]=-sin(angle);
  m[3]=0; m[7]=0; m[11]=sin(angle); m[15]=cos(angle);
}

void glmGenerateYZRotationMatrix4fDeg(float m[16], float angle)
{
  angle=glmRadf(angle);
  m[0]=cos(angle);  m[4]=0; m[8 ]=0; m[12]=sin(angle);
  m[1]=0;           m[5]=1; m[9 ]=0; m[13]=0;
  m[2]=0;           m[6]=0; m[10]=1; m[14]=0;
  m[3]=-sin(angle); m[7]=0; m[11]=0; m[15]=cos(angle);
}

void glmGenerateYZRotationMatrix4dDeg(double m[16], double angle)
{
  angle=glmRadd(angle);
  m[0]=cos(angle);  m[4]=0; m[8 ]=0; m[12]=sin(angle);
  m[1]=0;           m[5]=1; m[9 ]=0; m[13]=0;
  m[2]=0;           m[6]=0; m[10]=1; m[14]=0;
  m[3]=-sin(angle); m[7]=0; m[11]=0; m[15]=cos(angle);
}

void glmGenerateXZRotationMatrix4fDeg(float m[16], float angle)
{
  angle=glmRadf(angle);
  m[0]=1; m[4]=0;           m[8 ]=0; m[12]=0;
  m[1]=0; m[5]=cos(angle);  m[9 ]=0; m[13]=-sin(angle);
  m[2]=0; m[6]=0;           m[10]=1; m[14]=0;
  m[3]=0; m[7]=sin(angle);  m[11]=0; m[15]=cos(angle);
}

void glmGenerateXZRotationMatrix4dDeg(double m[16], double angle)
{
  angle=glmRadd(angle);
  m[0]=1; m[4]=0;           m[8 ]=0; m[12]=0;
  m[1]=0; m[5]=cos(angle);  m[9 ]=0; m[13]=-sin(angle);
  m[2]=0; m[6]=0;           m[10]=1; m[14]=0;
  m[3]=0; m[7]=sin(angle);  m[11]=0; m[15]=cos(angle);
}

void glmGenerateXWRotationMatrix4fDeg(float m[16], float angle)
{
  angle=glmRadf(angle);
  m[0]=1; m[4]=0;           m[8 ]=0;          m[12]=0;
  m[1]=0; m[5]=cos(angle);  m[9 ]=sin(angle); m[13]=0;
  m[2]=0; m[6]=-sin(angle); m[10]=cos(angle); m[14]=0;
  m[3]=0; m[7]=0;           m[11]=0;          m[15]=1;
}

void glmGenerateXWRotationMatrix4dDeg(double m[16], double angle)
{
  angle=glmRadd(angle);
  m[0]=1; m[4]=0;           m[8 ]=0;          m[12]=0;
  m[1]=0; m[5]=cos(angle);  m[9 ]=sin(angle); m[13]=0;
  m[2]=0; m[6]=-sin(angle); m[10]=cos(angle); m[14]=0;
  m[3]=0; m[7]=0;           m[11]=0;          m[15]=1;
}

void glmGenerateYWRotationMatrix4fDeg(float m[16], float angle)
{
  angle=glmRadf(angle);
  m[0]=cos(angle); m[4]=0; m[8 ]=-sin(angle); m[12]=0;
  m[1]=0;          m[5]=1; m[9 ]=0;           m[13]=0;
  m[2]=sin(angle); m[6]=0; m[10]=cos(angle);  m[14]=0;
  m[3]=0;          m[7]=0; m[11]=0;           m[15]=1;
}

void glmGenerateYWRotationMatrix4dDeg(double m[16], double angle)
{
  angle=glmRadd(angle);
  m[0]=cos(angle); m[4]=0; m[8 ]=-sin(angle); m[12]=0;
  m[1]=0;          m[5]=1; m[9 ]=0;           m[13]=0;
  m[2]=sin(angle); m[6]=0; m[10]=cos(angle);  m[14]=0;
  m[3]=0;          m[7]=0; m[11]=0;           m[15]=1;
}

void glmGenerateZWRotationMatrix4fDeg(float m[16], float angle)
{
  angle=glmRadf(angle);
  m[0]=cos(angle);  m[4]=sin(angle); m[8 ]=0; m[12]=0;
  m[1]=-sin(angle); m[5]=cos(angle); m[9 ]=0; m[13]=0;
  m[2]=0;           m[6]=0;          m[10]=1; m[14]=0;
  m[3]=0;           m[7]=0;          m[11]=0; m[15]=1;
}

void glmGenerateZWRotationMatrix4dDeg(double m[16], double angle)
{
  angle=glmRadd(angle);
  m[0]=cos(angle);  m[4]=sin(angle); m[8 ]=0; m[12]=0;
  m[1]=-sin(angle); m[5]=cos(angle); m[9 ]=0; m[13]=0;
  m[2]=0;           m[6]=0;          m[10]=1; m[14]=0;
  m[3]=0;           m[7]=0;          m[11]=0; m[15]=1;
}

void glmGenerateRotationMatrixFromAngles4fDeg(float matrix[16], float a1, float a2, float a3)
{
  float	sr, sp, sy, cr, cp, cy;
  a1=glmRadf(a1);
  a2=glmRadf(a2);
  a3=glmRadf(a3);

  sy=sin(a3);
  cy=cos(a3);
  sp=sin(a2);
  cp=cos(a2);
  sr=sin(a1);
  cr=cos(a1);

  matrix[0]=cp*cy;
  matrix[1]=cp*sy;
  matrix[2]=-sp;
  matrix[3]=0;
  matrix[4]=sr*sp*cy+cr*-sy;
  matrix[5]=sr*sp*sy+cr*cy;
  matrix[6]=sr*cp;
  matrix[7]=0;
  matrix[8]=(cr*sp*cy+-sr*-sy);
  matrix[9]=(cr*sp*sy+-sr*cy);
  matrix[10]=cr*cp;
  matrix[11]=0;
  matrix[12]=0.0;
  matrix[13]=0.0;
  matrix[14]=0.0;
  matrix[15]=1;
}

void glmGenerateRotationMatrixFromAngles4fRad(float matrix[16], float a1, float a2, float a3)
{
  float	sr, sp, sy, cr, cp, cy;

  sy=sin(a3);
  cy=cos(a3);
  sp=sin(a2);
  cp=cos(a2);
  sr=sin(a1);
  cr=cos(a1);

  matrix[0]=cp*cy;
  matrix[1]=cp*sy;
  matrix[2]=-sp;
  matrix[3]=0;
  matrix[4]=sr*sp*cy+cr*-sy;
  matrix[5]=sr*sp*sy+cr*cy;
  matrix[6]=sr*cp;
  matrix[7]=0;
  matrix[8]=(cr*sp*cy+-sr*-sy);
  matrix[9]=(cr*sp*sy+-sr*cy);
  matrix[10]=cr*cp;
  matrix[11]=0;
  matrix[12]=0.0;
  matrix[13]=0.0;
  matrix[14]=0.0;
  matrix[15]=1;
}

void glmGenerateRotationMatrixFromAngles4dDeg(double matrix[16], double a1, double a2, double a3)
{
  double sr, sp, sy, cr, cp, cy;
  a1=glmRadd(a1);
  a2=glmRadd(a2);
  a3=glmRadd(a3);

  sy=sin(a3);
  cy=cos(a3);
  sp=sin(a2);
  cp=cos(a2);
  sr=sin(a1);
  cr=cos(a1);

  matrix[0]=cp*cy;
  matrix[1]=cp*sy;
  matrix[2]=-sp;
  matrix[3]=0;
  matrix[4]=sr*sp*cy+cr*-sy;
  matrix[5]=sr*sp*sy+cr*cy;
  matrix[6]=sr*cp;
  matrix[7]=0;
  matrix[8]=(cr*sp*cy+-sr*-sy);
  matrix[9]=(cr*sp*sy+-sr*cy);
  matrix[10]=cr*cp;
  matrix[11]=0;
  matrix[12]=0.0;
  matrix[13]=0.0;
  matrix[14]=0.0;
  matrix[15]=1;
}

void glmGenerateRotationMatrixFromAngles4dRad(double matrix[16], double a1, double a2, double a3)
{
  double sr, sp, sy, cr, cp, cy;

  sy=sin(a3);
  cy=cos(a3);
  sp=sin(a2);
  cp=cos(a2);
  sr=sin(a1);
  cr=cos(a1);

  matrix[0]=cp*cy;
  matrix[1]=cp*sy;
  matrix[2]=-sp;
  matrix[3]=0;
  matrix[4]=sr*sp*cy+cr*-sy;
  matrix[5]=sr*sp*sy+cr*cy;
  matrix[6]=sr*cp;
  matrix[7]=0;
  matrix[8]=(cr*sp*cy+-sr*-sy);
  matrix[9]=(cr*sp*sy+-sr*cy);
  matrix[10]=cr*cp;
  matrix[11]=0;
  matrix[12]=0.0;
  matrix[13]=0.0;
  matrix[14]=0.0;
  matrix[15]=1;
}

void glmGenerateTranslationMatrix3f(float m[9], const float p[2])
{
  m[0]=1; m[3]=0; m[6]=p[0];
  m[1]=0; m[4]=1; m[7]=p[1];
  m[2]=0; m[5]=0; m[8]=1;
}

void glmGenerateTranslationMatrix3d(double m[9], const double p[2])
{
  m[0]=1; m[3]=0; m[6]=p[0];
  m[1]=0; m[4]=1; m[7]=p[1];
  m[2]=0; m[5]=0; m[8]=1;
}

// Generates OpenGL style 4x4 column major translation matrix.
void glmGenerateTranslationMatrix4f(float m[16], const float p[3])
{
  //  | 1  0  0  X |
  //M=| 0  1  0  Y |
  //  | 0  0  1  Z |
  //  | 0  0  0  1 |
  m[0]=1; m[4]=0; m[8]=0;  m[12]=p[0];
  m[1]=0; m[5]=1; m[9]=0;  m[13]=p[1];
  m[2]=0; m[6]=0; m[10]=1; m[14]=p[2];
  m[3]=0; m[7]=0; m[11]=0; m[15]=1;
}

// Generates OpenGL style 4x4 column major translation matrix.
void glmGenerateTranslationMatrix4d(double m[16], const double p[3])
{
  //  | 1  0  0  X |
  //M=| 0  1  0  Y |
  //  | 0  0  1  Z |
  //  | 0  0  0  1 |
  m[0]=1; m[4]=0; m[8]=0;  m[12]=p[0];
  m[1]=0; m[5]=1; m[9]=0;  m[13]=p[1];
  m[2]=0; m[6]=0; m[10]=1; m[14]=p[2];
  m[3]=0; m[7]=0; m[11]=0; m[15]=1;
}

void glmGenerateTranslationMatrix5f(float m[25], const float p[4])
{
  m[0]=1; m[5]=0; m[10]=0; m[15]=0; m[20]=p[0];
  m[1]=0; m[6]=1; m[11]=0; m[16]=0; m[21]=p[1];
  m[2]=0; m[7]=0; m[12]=1; m[17]=0; m[22]=p[2];
  m[3]=0; m[8]=0; m[13]=0; m[18]=1; m[23]=p[3];
  m[4]=0; m[9]=0; m[14]=0; m[19]=0; m[24]=1;
}

void glmGenerateTranslationMatrix5d(double m[25], const double p[4])
{
  m[0]=1; m[5]=0; m[10]=0; m[15]=0; m[20]=p[0];
  m[1]=0; m[6]=1; m[11]=0; m[16]=0; m[21]=p[1];
  m[2]=0; m[7]=0; m[12]=1; m[17]=0; m[22]=p[2];
  m[3]=0; m[8]=0; m[13]=0; m[18]=1; m[23]=p[3];
  m[4]=0; m[9]=0; m[14]=0; m[19]=0; m[24]=1;
}

// Generates 2x2 identity matrix.
void glmGenerateIdentityMatrix2f(float m[4])
{
  m[0]=1; m[2]=0;
  m[1]=0; m[3]=1;
}

// Generates 2x2 identity matrix.
void glmGenerateIdentityMatrix2d(double m[4])
{
  m[0]=1; m[2]=0;
  m[1]=0; m[3]=1;
}

// Generates 3x3 identity matrix.
void glmGenerateIdentityMatrix3f(float m[9])
{
  m[0]=1; m[3]=0; m[6]=0;
  m[1]=0; m[4]=1; m[7]=0;
  m[2]=0; m[5]=0; m[8]=1;
}

// Generates 3x3 colum major identity matrix.
void glmGenerateIdentityMatrix3d(double m[9])
{
  m[0]=1; m[3]=0; m[6]=0;
  m[1]=0; m[4]=1; m[7]=0;
  m[2]=0; m[5]=0; m[8]=1;
}

// Generates 4x4 identity matrix.
void glmGenerateIdentityMatrix4f(float m[16])
{
  m[0]=1; m[4]=0; m[8]=0;  m[12]=0;
  m[1]=0; m[5]=1; m[9]=0;  m[13]=0;
  m[2]=0; m[6]=0; m[10]=1; m[14]=0;
  m[3]=0; m[7]=0; m[11]=0; m[15]=1;
}

// Generates 4x4 identity matrix.
void glmGenerateIdentityMatrix4d(double m[16])
{
  m[0]=1; m[4]=0; m[8]=0;  m[12]=0;
  m[1]=0; m[5]=1; m[9]=0;  m[13]=0;
  m[2]=0; m[6]=0; m[10]=1; m[14]=0;
  m[3]=0; m[7]=0; m[11]=0; m[15]=1;
}

// Generates a 2x2 zero matrix.
void glmGenerateZeroMatrix2f(float m[4])
{
  memset(m, 0, 4*sizeof(float));
}

// Generates a 2x2 zero matrix.
void glmGenerateZeroMatrix2d(double m[4])
{
  memset(m, 0, 4*sizeof(double));
}

// Generates 3x3 zero matrix.
void glmGenerateZeroMatrix3f(float m[9])
{
  memset(m, 0, 9*sizeof(float));
}

// Generates 3x3 zero matrix.
void glmGenerateZeroMatrix3d(double m[9])
{
  memset(m, 0, 9*sizeof(double));
}

// Generates 4x4 zero matrix.
void glmGenerateZeroMatrix4f(float m[16])
{
  memset(m, 0, 16*sizeof(float));
}

// Generates 4x4 zero matrix.
void glmGenerateZeroMatrix4d(double m[16])
{
  memset(m, 0, 16*sizeof(double));
}

void glmMultVectorByMatrix2f(float result[2], const float v[2], const float m[4])
{
  result[0]=m[0]*v[0]+m[2]*v[1];
  result[1]=m[1]*v[0]+m[3]*v[1];
}

void glmMultVectorByMatrix2d(double result[2], const double v[2], const double m[4])
{
  result[0]=m[0]*v[0]+m[2]*v[1];
  result[1]=m[1]*v[0]+m[3]*v[1];
}

void glmMultVectorByMatrix32f(float result[2], const float v[2], const float m[9])
{
  result[0]=m[0]*v[0]+m[3]*v[1]+m[6];
  result[1]=m[1]*v[0]+m[4]*v[1]+m[7];
}

void glmMultVectorByMatrix32d(double result[2], const double v[2], const double m[9])
{
  result[0]=m[0]*v[0]+m[3]*v[1]+m[6];
  result[1]=m[1]*v[0]+m[4]*v[1]+m[7];
}

void glmMultVectorByMatrix3f(float result[3], const float v[3], const float m[9])
{
  result[0]=m[0]*v[0]+m[3]*v[1]+m[6]*v[2];
  result[1]=m[1]*v[0]+m[4]*v[1]+m[7]*v[2];
  result[2]=m[2]*v[0]+m[5]*v[1]+m[8]*v[2];
}

void glmMultVectorByMatrix3d(double result[3], const double v[3], const double m[9])
{
  result[0]=m[0]*v[0]+m[3]*v[1]+m[6]*v[2];
  result[1]=m[1]*v[0]+m[4]*v[1]+m[7]*v[2];
  result[2]=m[2]*v[0]+m[5]*v[1]+m[8]*v[2];
}

void glmMultVectorByMatrix43f(float result[3], const float v[3], const float m[16])
{
  result[0]=m[0]*v[0]+m[4]*v[1]+m[8 ]*v[2]+m[12];
  result[1]=m[1]*v[0]+m[5]*v[1]+m[9 ]*v[2]+m[13];
  result[2]=m[2]*v[0]+m[6]*v[1]+m[10]*v[2]+m[14];
}

void glmMultVectorByMatrix43d(double result[3], const double v[3], const double m[16])
{
  result[0]=m[0]*v[0]+m[4]*v[1]+m[8 ]*v[2]+m[12];
  result[1]=m[1]*v[0]+m[5]*v[1]+m[9 ]*v[2]+m[13];
  result[2]=m[2]*v[0]+m[6]*v[1]+m[10]*v[2]+m[14];
}

void glmMultVectorBySubMatrix43f(float result[3], const float v[3], const float m[16])
{
  result[0]=m[0]*v[0]+m[4]*v[1]+m[8 ]*v[2];
  result[1]=m[1]*v[0]+m[5]*v[1]+m[9 ]*v[2];
  result[2]=m[2]*v[0]+m[6]*v[1]+m[10]*v[2];
}

void glmMultVectorBySubMatrix43d(double result[3], const double v[3], const double m[16])
{
  result[0]=m[0]*v[0]+m[4]*v[1]+m[8 ]*v[2];
  result[1]=m[1]*v[0]+m[5]*v[1]+m[9 ]*v[2];
  result[2]=m[2]*v[0]+m[6]*v[1]+m[10]*v[2];
}

void glmMultVectorByMatrix4f(float result[4], const float v[4], const float m[16])
{
  result[0]=m[0]*v[0]+m[4]*v[1]+m[8 ]*v[2]+m[12]*v[3];
  result[1]=m[1]*v[0]+m[5]*v[1]+m[9 ]*v[2]+m[13]*v[3];
  result[2]=m[2]*v[0]+m[6]*v[1]+m[10]*v[2]+m[14]*v[3];
  result[3]=m[3]*v[0]+m[7]*v[1]+m[11]*v[2]+m[15]*v[3];
}

void glmMultVectorByMatrix4d(double result[4], const double v[4], const double m[16])
{
  result[0]=m[0]*v[0]+m[4]*v[1]+m[8 ]*v[2]+m[12]*v[3];
  result[1]=m[1]*v[0]+m[5]*v[1]+m[9 ]*v[2]+m[13]*v[3];
  result[2]=m[2]*v[0]+m[6]*v[1]+m[10]*v[2]+m[14]*v[3];
  result[3]=m[3]*v[0]+m[7]*v[1]+m[11]*v[2]+m[15]*v[3];
}

void glmMultVectorByMatrix5f(float result[5], const float v[5], const float m[25])
{
  result[0]=m[0]*v[0]+m[5]*v[1]+m[10]*v[2]+m[15]*v[3]+m[20]*v[4];
  result[1]=m[1]*v[0]+m[6]*v[1]+m[11]*v[2]+m[16]*v[3]+m[21]*v[4];
  result[2]=m[2]*v[0]+m[7]*v[1]+m[12]*v[2]+m[17]*v[3]+m[22]*v[4];
  result[3]=m[3]*v[0]+m[8]*v[1]+m[13]*v[2]+m[18]*v[3]+m[23]*v[4];
  result[4]=m[4]*v[0]+m[9]*v[1]+m[14]*v[2]+m[19]*v[3]+m[24]*v[4];
}

void glmMultVectorByMatrix5d(double result[5], const double v[5], const double m[25])
{
  result[0]=m[0]*v[0]+m[5]*v[1]+m[10]*v[2]+m[15]*v[3]+m[20]*v[4];
  result[1]=m[1]*v[0]+m[6]*v[1]+m[11]*v[2]+m[16]*v[3]+m[21]*v[4];
  result[2]=m[2]*v[0]+m[7]*v[1]+m[12]*v[2]+m[17]*v[3]+m[22]*v[4];
  result[3]=m[3]*v[0]+m[8]*v[1]+m[13]*v[2]+m[18]*v[3]+m[23]*v[4];
  result[4]=m[4]*v[0]+m[9]*v[1]+m[14]*v[2]+m[19]*v[3]+m[24]*v[4];
}

void glmMultVectorByMatrix54f(float result[4], const float v[4], const float m[25])
{
  result[0]=m[0]*v[0]+m[5]*v[1]+m[10]*v[2]+m[15]*v[3]+m[20];
  result[1]=m[1]*v[0]+m[6]*v[1]+m[11]*v[2]+m[16]*v[3]+m[21];
  result[2]=m[2]*v[0]+m[7]*v[1]+m[12]*v[2]+m[17]*v[3]+m[22];
  result[3]=m[3]*v[0]+m[8]*v[1]+m[13]*v[2]+m[18]*v[3]+m[23];
  result[4]=m[4]*v[0]+m[9]*v[1]+m[14]*v[2]+m[19]*v[3]+m[24];
}

void glmMultVectorByMatrix54d(double result[4], const double v[4], const double m[25])
{
  result[0]=m[0]*v[0]+m[5]*v[1]+m[10]*v[2]+m[15]*v[3]+m[20];
  result[1]=m[1]*v[0]+m[6]*v[1]+m[11]*v[2]+m[16]*v[3]+m[21];
  result[2]=m[2]*v[0]+m[7]*v[1]+m[12]*v[2]+m[17]*v[3]+m[22];
  result[3]=m[3]*v[0]+m[8]*v[1]+m[13]*v[2]+m[18]*v[3]+m[23];
  result[4]=m[4]*v[0]+m[9]*v[1]+m[14]*v[2]+m[19]*v[3]+m[24];
}

void glmGenerateInertiaMatrixForBox3f(float m[9], float x, float y, float z, float mass)
{
  mass/=12.0f;
  m[0]=mass*(y*y+z*z); m[3]=0;              m[6]=0;
  m[1]=0;              m[4]=mass*(x*x+z*z); m[7]=0;
  m[2]=0;              m[5]=0;              m[8]=mass*(x*x+y*y);
}

void glmGenerateInertiaMatrixForBox3d(double m[9], double x, double y, double z, double mass)
{
  mass/=12.0f;
  m[0]=mass*(y*y+z*z); m[3]=0;              m[6]=0;
  m[1]=0;              m[4]=mass*(x*x+z*z); m[7]=0;
  m[2]=0;              m[5]=0;              m[8]=mass*(x*x+y*y);
}

void glmGenerateInertiaMatrixForBox4f(float m[16], float x, float y, float z, float mass)
{
  mass/=12.0f;
  m[0]=mass*(y*y+z*z); m[4]=0;              m[8]=0;               m[12]=0;
  m[1]=0;              m[5]=mass*(x*x+z*z); m[9]=0;               m[13]=0;
  m[2]=0;              m[6]=0;              m[10]=mass*(x*x+y*y); m[14]=0;
  m[3]=0;              m[7]=0;              m[11]=0;              m[15]=1;
}

void glmGenerateInertiaMatrixForBox4d(double m[16], double x, double y, double z, double mass)
{
  mass/=12.0;
  m[0]=mass*(y*y+z*z); m[4]=0;              m[8]=0;               m[12]=0;
  m[1]=0;              m[5]=mass*(x*x+z*z); m[9]=0;               m[13]=0;
  m[2]=0;              m[6]=0;              m[10]=mass*(x*x+y*y); m[14]=0;
  m[3]=0;              m[7]=0;              m[11]=0;              m[15]=1;
}

void glmGenerateInertiaMatrixForSphere3f(float m[9], float r, float mass)
{
  mass*=2.0f/5.0f;
  m[0]=mass*r*r; m[3]=0;        m[6]=0;
  m[1]=0;        m[4]=mass*r*r; m[7]=0;
  m[2]=0;        m[5]=0;        m[8]=mass*r*r;
}

void glmGenerateInertiaMatrixForSphere3d(double m[9], double r, double mass)
{
  mass*=2.0f/5.0f;
  m[0]=mass*r*r; m[3]=0;        m[6]=0;
  m[1]=0;        m[4]=mass*r*r; m[7]=0;
  m[2]=0;        m[5]=0;        m[8]=mass*r*r;
}

void glmGenerateInertiaMatrixForSphere4f(float m[16], float r, float mass)
{
  mass*=2.0f/5.0f;
  m[0]=mass*r*r; m[4]=0;        m[8]=0;         m[12]=0;
  m[1]=0;        m[5]=mass*r*r; m[9]=0;         m[13]=0;
  m[2]=0;        m[6]=0;        m[10]=mass*r*r; m[14]=0;
  m[3]=0;        m[7]=0;        m[11]=0;        m[15]=1;
}

void glmGenerateInertiaMatrixForSphere4d(double m[16], double r, double mass)
{
  mass*=2.0/5.0;
  m[0]=mass*r*r; m[4]=0;        m[8]=0;         m[12]=0;
  m[1]=0;        m[5]=mass*r*r; m[9]=0;         m[13]=0;
  m[2]=0;        m[6]=0;        m[10]=mass*r*r; m[14]=0;
  m[3]=0;        m[7]=0;        m[11]=0;        m[15]=1;
}

void glmGenerateInertiaMatrixForCylinder3f(float m[9], float r, float h, float mass)
{
  m[0]=(1.0f/12.0f)*mass*h*h+(1.0f/4.0f)*mass*r*r; m[3]=0;  m[6]=0;
  m[1]=0;  m[4]=(1.0f/12.0f)*mass*h*h+(1.0f/4.0f)*mass*r*r; m[7]=0;
  m[2]=0;  m[5]=0;        m[8]=0.5f*mass*r*r;
}

void glmGenerateInertiaMatrixForCylinder3d(double m[9], double r, double h, double mass)
{
  m[0]=(1.0/12.0)*mass*h*h+(1.0/4.0)*mass*r*r; m[3]=0;  m[6]=0;
  m[1]=0;  m[4]=(1.0/12.0)*mass*h*h+(1.0/4.0)*mass*r*r; m[7]=0;
  m[2]=0;  m[5]=0;        m[8]=0.5f*mass*r*r;
}


void glmGenerateInertiaMatrixForCylinder4f(float m[16], float r, float h, float mass)
{
  m[0]=(1.0f/12.0f)*mass*h*h+(1.0f/4.0f)*mass*r*r; m[4]=0;  m[8]=0;  m[12]=0;
  m[1]=0;  m[5]=(1.0f/12.0f)*mass*h*h+(1.0f/4.0f)*mass*r*r; m[9]=0;  m[13]=0;
  m[2]=0;  m[6]=0;        m[10]=0.5f*mass*r*r;              m[14]=0;
  m[3]=0;  m[7]=0;                                          m[11]=0; m[15]=1;
}

void glmGenerateInertiaMatrixForCylinder4d(double m[16], double r, double h, double mass)
{
  m[0]=(1.0/12.0)*mass*h*h+(1.0/4.0)*mass*r*r; m[4]=0;  m[8]=0;  m[12]=0;
  m[1]=0;  m[5]=(1.0/12.0)*mass*h*h+(1.0/4.0)*mass*r*r; m[9]=0;  m[13]=0;
  m[2]=0;  m[6]=0;        m[10]=0.5f*mass*r*r;              m[14]=0;
  m[3]=0;  m[7]=0;                                          m[11]=0; m[15]=1;
}

void glmGenerateInertiaMatrixForTorus3f(float m[9], float a, float c, float mass)
{
  m[0]=0.125f*(5*a*a+4*c*c)*mass; m[3]=0;        m[6]=0;
  m[1]=0;        m[4]=0.125f*(5*a*a+4*c*c)*mass; m[7]=0;
  m[2]=0;        m[5]=0;        m[8]=(0.75f*a*a+c*c)*mass;
}

void glmGenerateInertiaMatrixForTorus3d(double m[9], double a, double c, double mass)
{
  m[0]=0.125*(5*a*a+4*c*c)*mass; m[3]=0;        m[6]=0;
  m[1]=0;        m[4]=0.125*(5*a*a+4*c*c)*mass; m[7]=0;
  m[2]=0;        m[5]=0;        m[8]=(0.75*a*a+c*c)*mass;
}


void glmGenerateInertiaMatrixForTorus4f(float m[16], float a, float c, float mass)
{
  m[0]=0.125f*(5*a*a+4*c*c)*mass; m[4]=0;        m[8]=0;  m[12]=0;
  m[1]=0;        m[5]=0.125f*(5*a*a+4*c*c)*mass; m[9]=0;  m[13]=0;
  m[2]=0;        m[6]=0;        m[10]=(0.75f*a*a+c*c)*mass;           m[14]=0;
  m[3]=0;        m[7]=0;                         m[11]=0; m[15]=1;
}

void glmGenerateInertiaMatrixForTorus4d(double m[16], double a, double c, double mass)
{
  m[0]=0.125*(5*a*a+4*c*c)*mass; m[4]=0;        m[8]=0;  m[12]=0;
  m[1]=0;        m[5]=0.125*(5*a*a+4*c*c)*mass; m[9]=0;  m[13]=0;
  m[2]=0;        m[6]=0;        m[10]=(0.75*a*a+c*c)*mass;           m[14]=0;
  m[3]=0;        m[7]=0;                         m[11]=0; m[15]=1;
}

void glmGenerateScaleMatrix2f(float m[4], const float s[2])
{
  m[0]=s[0]; m[2]=0;
  m[1]=0;    m[3]=s[1];
}

void glmGenerateScaleMatrix2d(double m[4], const double s[2])
{
  m[0]=s[0]; m[2]=0;
  m[1]=0;    m[3]=s[1];
}

void glmGenerateScaleMatrix3f(float m[9], const float s[3])
{
  m[0]=s[0]; m[3]=0;    m[6]=0;
  m[1]=0;    m[4]=s[1]; m[7]=0;
  m[2]=0;    m[5]=0;    m[8]=s[2];
}

void glmGenerateScaleMatrix3d(double m[9], const double s[3])
{
  m[0]=s[0]; m[3]=0;    m[6]=0;
  m[1]=0;    m[4]=s[1]; m[7]=0;
  m[2]=0;    m[5]=0;    m[8]=s[2];
}

void glmGenerateScaleMatrix4f(float m[16], const float s[3])
{
  m[0]=s[0]; m[4]=0;    m[8]=0;     m[12]=0;
  m[1]=0;    m[5]=s[1]; m[9]=0;     m[13]=0;
  m[2]=0;    m[6]=0;    m[10]=s[2]; m[14]=0;
  m[3]=0;    m[7]=0;    m[11]=0;    m[15]=1;
}

void glmGenerateScaleMatrix4d(double m[16], const double s[3])
{
  m[0]=s[0]; m[4]=0;    m[8]=0;     m[12]=0;
  m[1]=0;    m[5]=s[1]; m[9]=0;     m[13]=0;
  m[2]=0;    m[6]=0;    m[10]=s[2]; m[14]=0;
  m[3]=0;    m[7]=0;    m[11]=0;    m[15]=1;
}

void glmGenerateAxisAngleRotationMatrix3fRad(float m[9], const float axis[3], float angle)
{
  float rcos, rsin;
  rcos = cos(angle);
  rsin = sin(angle);

  m[0] =          rcos+ axis[0]*axis[0]*(1.0f-rcos);
  m[3] =  axis[2]*rsin+ axis[1]*axis[0]*(1.0f-rcos);
  m[6] = -axis[1]*rsin+ axis[2]*axis[0]*(1.0f-rcos);
  m[1] = -axis[2]*rsin+ axis[0]*axis[1]*(1.0f-rcos);
  m[4] =          rcos+ axis[1]*axis[1]*(1.0f-rcos);
  m[7] =  axis[0]*rsin+ axis[2]*axis[1]*(1.0f-rcos);
  m[2] =  axis[1]*rsin+ axis[0]*axis[2]*(1.0f-rcos);
  m[5] = -axis[0]*rsin+ axis[1]*axis[2]*(1.0f-rcos);
  m[8]=           rcos+ axis[2]*axis[2]*(1.0f-rcos);
}

void glmGenerateAxisAngleRotationMatrix3fDeg(float m[9], const float axis[3], float angle)
{
  float rcos, rsin;
  angle=glmRadf(angle);
  rcos = cos(angle);
  rsin = sin(angle);

  m[0] =          rcos+ axis[0]*axis[0]*(1.0f-rcos);
  m[3] =  axis[2]*rsin+ axis[1]*axis[0]*(1.0f-rcos);
  m[6] = -axis[1]*rsin+ axis[2]*axis[0]*(1.0f-rcos);
  m[1] = -axis[2]*rsin+ axis[0]*axis[1]*(1.0f-rcos);
  m[4] =          rcos+ axis[1]*axis[1]*(1.0f-rcos);
  m[7] =  axis[0]*rsin+ axis[2]*axis[1]*(1.0f-rcos);
  m[2] =  axis[1]*rsin+ axis[0]*axis[2]*(1.0f-rcos);
  m[5] = -axis[0]*rsin+ axis[1]*axis[2]*(1.0f-rcos);
  m[8]=           rcos+ axis[2]*axis[2]*(1.0f-rcos);
}

void glmGenerateAxisAngleRotationMatrix3dRad(double m[9], const double axis[3], double angle)
{
  double rcos, rsin;
  rcos = cos(angle);
  rsin = sin(angle);

  m[0] =          rcos+ axis[0]*axis[0]*(1.0f-rcos);
  m[3] =  axis[2]*rsin+ axis[1]*axis[0]*(1.0f-rcos);
  m[6] = -axis[1]*rsin+ axis[2]*axis[0]*(1.0f-rcos);
  m[1] = -axis[2]*rsin+ axis[0]*axis[1]*(1.0f-rcos);
  m[4] =          rcos+ axis[1]*axis[1]*(1.0f-rcos);
  m[7] =  axis[0]*rsin+ axis[2]*axis[1]*(1.0f-rcos);
  m[2] =  axis[1]*rsin+ axis[0]*axis[2]*(1.0f-rcos);
  m[5] = -axis[0]*rsin+ axis[1]*axis[2]*(1.0f-rcos);
  m[8]=           rcos+ axis[2]*axis[2]*(1.0f-rcos);
}

void glmGenerateAxisAngleRotationMatrix3dDeg(double m[9], const double axis[3], double angle)
{
  double rcos, rsin;
  angle=glmRadd(angle);
  rcos = cos(angle);
  rsin = sin(angle);

  m[0] =          rcos+ axis[0]*axis[0]*(1.0f-rcos);
  m[3] =  axis[2]*rsin+ axis[1]*axis[0]*(1.0f-rcos);
  m[6] = -axis[1]*rsin+ axis[2]*axis[0]*(1.0f-rcos);
  m[1] = -axis[2]*rsin+ axis[0]*axis[1]*(1.0f-rcos);
  m[4] =          rcos+ axis[1]*axis[1]*(1.0f-rcos);
  m[7] =  axis[0]*rsin+ axis[2]*axis[1]*(1.0f-rcos);
  m[2] =  axis[1]*rsin+ axis[0]*axis[2]*(1.0f-rcos);
  m[5] = -axis[0]*rsin+ axis[1]*axis[2]*(1.0f-rcos);
  m[8]=           rcos+ axis[2]*axis[2]*(1.0f-rcos);
}

void glmGenerateAxisAngleRotationMatrix4fRad(float m[16], const float axis[3], float angle)
{
  float rcos = cos(angle);
  float rsin = sin(angle);

  m[0] =          rcos+ axis[0]*axis[0]*(1.0f-rcos);
  m[4] =  axis[2]*rsin+ axis[1]*axis[0]*(1.0f-rcos);
  m[8] = -axis[1]*rsin+ axis[2]*axis[0]*(1.0f-rcos);
  m[12]= 0;
  m[1] = -axis[2]*rsin+ axis[0]*axis[1]*(1.0f-rcos);
  m[5] =          rcos+ axis[1]*axis[1]*(1.0f-rcos);
  m[9] =  axis[0]*rsin+ axis[2]*axis[1]*(1.0f-rcos);
  m[13]= 0;
  m[2] =  axis[1]*rsin+ axis[0]*axis[2]*(1.0f-rcos);
  m[6] = -axis[0]*rsin+ axis[1]*axis[2]*(1.0f-rcos);
  m[10]=          rcos+ axis[2]*axis[2]*(1.0f-rcos);
  m[14]= 0;
  m[3] = 0;
  m[7] = 0;
  m[11]= 0;
  m[15]= 1;
}

void glmGenerateAxisAngleRotationMatrix4fDeg(float m[16], const float axis[3], float angle)
{
  float rcos, rsin;
  angle=glmRadf(angle);
  rcos = cos(angle);
  rsin = sin(angle);

  m[0] =          rcos+ axis[0]*axis[0]*(1.0f-rcos);
  m[4] =  axis[2]*rsin+ axis[1]*axis[0]*(1.0f-rcos);
  m[8] = -axis[1]*rsin+ axis[2]*axis[0]*(1.0f-rcos);
  m[12]= 0;
  m[1] = -axis[2]*rsin+ axis[0]*axis[1]*(1.0f-rcos);
  m[5] =          rcos+ axis[1]*axis[1]*(1.0f-rcos);
  m[9] =  axis[0]*rsin+ axis[2]*axis[1]*(1.0f-rcos);
  m[13]= 0;
  m[2] =  axis[1]*rsin+ axis[0]*axis[2]*(1.0f-rcos);
  m[6] = -axis[0]*rsin+ axis[1]*axis[2]*(1.0f-rcos);
  m[10]=          rcos+ axis[2]*axis[2]*(1.0f-rcos);
  m[14]= 0;
  m[3] = 0;
  m[7] = 0;
  m[11]= 0;
  m[15]= 1;
}

void glmGenerateAxisAngleRotationMatrix4dRad(double m[16], const double axis[3], double angle)
{
  double rcos, rsin;
  rcos = cos(angle);
  rsin = sin(angle);

  m[0] =          rcos+ axis[0]*axis[0]*(1.0-rcos);
  m[4] =  axis[2]*rsin+ axis[1]*axis[0]*(1.0-rcos);
  m[8] = -axis[1]*rsin+ axis[2]*axis[0]*(1.0-rcos);
  m[12]= 0;
  m[1] = -axis[2]*rsin+ axis[0]*axis[1]*(1.0-rcos);
  m[5] =          rcos+ axis[1]*axis[1]*(1.0-rcos);
  m[9] =  axis[0]*rsin+ axis[2]*axis[1]*(1.0-rcos);
  m[13]= 0;
  m[2] =  axis[1]*rsin+ axis[0]*axis[2]*(1.0-rcos);
  m[6] = -axis[0]*rsin+ axis[1]*axis[2]*(1.0-rcos);
  m[10]=          rcos+ axis[2]*axis[2]*(1.0-rcos);
  m[14]= 0;
  m[3] = 0;
  m[7] = 0;
  m[11]= 0;
  m[15]= 1;
}

void glmGenerateAxisAngleRotationMatrix4dDeg(double m[16], const double axis[3], double angle)
{
  double rcos, rsin;
  angle=glmRadd(angle);
  rcos = cos(angle);
  rsin = sin(angle);

  m[0] =          rcos+ axis[0]*axis[0]*(1.0-rcos);
  m[4] =  axis[2]*rsin+ axis[1]*axis[0]*(1.0-rcos);
  m[8] = -axis[1]*rsin+ axis[2]*axis[0]*(1.0-rcos);
  m[12]= 0;
  m[1] = -axis[2]*rsin+ axis[0]*axis[1]*(1.0-rcos);
  m[5] =          rcos+ axis[1]*axis[1]*(1.0-rcos);
  m[9] =  axis[0]*rsin+ axis[2]*axis[1]*(1.0-rcos);
  m[13]= 0;
  m[2] =  axis[1]*rsin+ axis[0]*axis[2]*(1.0-rcos);
  m[6] = -axis[0]*rsin+ axis[1]*axis[2]*(1.0-rcos);
  m[10]=          rcos+ axis[2]*axis[2]*(1.0-rcos);
  m[14]= 0;
  m[3] = 0;
  m[7] = 0;
  m[11]= 0;
  m[15]= 1;
}

void glmOrthonormalizeMatrix3f(float m[9])
{
  glmNormalizeVector3f(m);
  glmCrossProduct3f(&m[6], m, &m[3]);
  glmNormalizeVector3f(&m[6]);
  glmCrossProduct3f(&m[3], &m[6], m);
  glmNormalizeVector3f(&m[3]);
}

void glmOrthonormalizeMatrix3d(double m[9])
{
  glmNormalizeVector3d(m);
  glmCrossProduct3d(&m[6], m, &m[3]);
  glmNormalizeVector3d(&m[6]);
  glmCrossProduct3d(&m[3], &m[6], m);
  glmNormalizeVector3d(&m[3]);
}

void glmPrintMatrix2f(const float m[4])
{
  printf("[ %f %f ]\n",   m[0], m[2]);
  printf("[ %f %f ]\n\n", m[1], m[3]);
}

void glmPrintMatrix2d(const double m[4])
{
  printf("[ %f %f ]\n",   m[0], m[2]);
  printf("[ %f %f ]\n\n", m[1], m[3]);
}

void glmPrintMatrix3f(const float m[9])
{
  printf("[ %f %f %f ]\n",   m[0], m[3], m[6]);
  printf("[ %f %f %f ]\n",   m[1], m[4], m[7]);
  printf("[ %f %f %f ]\n\n", m[2], m[5], m[8]);
}

void glmPrintMatrix3d(const double m[9])
{
  printf("[ %f %f %f ]\n",   m[0], m[3], m[6]);
  printf("[ %f %f %f ]\n",   m[1], m[4], m[7]);
  printf("[ %f %f %f ]\n\n", m[2], m[5], m[8]);
}

void glmPrintMatrix4f(const float m[16])
{
  printf("[ %f %f %f %f ]\n", m[0], m[4], m[8],  m[12]);
  printf("[ %f %f %f %f ]\n", m[1], m[5], m[9],  m[13]);
  printf("[ %f %f %f %f ]\n", m[2], m[6], m[10], m[14]);
  printf("[ %f %f %f %f ]\n\n", m[3], m[7], m[11], m[15]);
}

void glmPrintMatrix4d(const double m[16])
{
  printf("[ %f %f %f %f ]\n", m[0], m[4], m[8],  m[12]);
  printf("[ %f %f %f %f ]\n", m[1], m[5], m[9],  m[13]);
  printf("[ %f %f %f %f ]\n", m[2], m[6], m[10], m[14]);
  printf("[ %f %f %f %f ]\n\n", m[3], m[7], m[11], m[15]);
}

void glmExpandMatrix24f(float result[16], const float m[4])
{
  result[0]=m[0];
  result[1]=m[1];
  result[2]=0;
  result[3]=0;

  result[4]=m[2];
  result[5]=m[3];
  result[6]=0;
  result[7]=0;

  result[8]=0;
  result[9]=0;
  result[10]=0;
  result[11]=0;

  result[12]=0;
  result[13]=0;
  result[14]=0;
  result[15]=1;
}

void glmExpandMatrix24d(double result[16], const double m[4])
{
  result[0]=m[0];
  result[1]=m[1];
  result[2]=0;
  result[3]=0;

  result[4]=m[2];
  result[5]=m[3];
  result[6]=0;
  result[7]=0;

  result[8]=0;
  result[9]=0;
  result[10]=0;
  result[11]=0;

  result[12]=0;
  result[13]=0;
  result[14]=0;
  result[15]=1;
}

void glmExpandMatrix23f(float result[9], const float m[4])
{
  result[0]=m[0];
  result[1]=m[1];
  result[2]=0;

  result[3]=m[2];
  result[4]=m[3];
  result[5]=0;

  result[6]=0;
  result[7]=0;
  result[8]=1;
}

void glmExpandMatrix23d(double result[9], const double m[4])
{
  result[0]=m[0];
  result[1]=m[1];
  result[2]=0;

  result[3]=m[2];
  result[4]=m[3];
  result[5]=0;

  result[6]=0;
  result[7]=0;
  result[8]=1;
}

void glmExpandMatrix34f(float result[16], const float m[9])
{
  result[0]=m[0];
  result[1]=m[1];
  result[2]=m[2];
  result[3]=0;

  result[4]=m[3];
  result[5]=m[4];
  result[6]=m[5];
  result[7]=0;

  result[8]=m[6];
  result[9]=m[7];
  result[10]=m[8];
  result[11]=0;

  result[12]=0;
  result[13]=0;
  result[14]=0;
  result[15]=1;
}

void glmExpandMatrix34d(double result[16], const double m[9])
{
  result[0]=m[0];
  result[1]=m[1];
  result[2]=m[2];
  result[3]=0;

  result[4]=m[3];
  result[5]=m[4];
  result[6]=m[5];
  result[7]=0;

  result[8]=m[6];
  result[9]=m[7];
  result[10]=m[8];
  result[11]=0;

  result[12]=0;
  result[13]=0;
  result[14]=0;
  result[15]=1;
}


void glmTransformMeshf(float *vertex, DWORD vStride,
                       float *normal, DWORD nStride,
                       DWORD vCount, const float m[16])
{
  float temp[3];
  BYTE *v, *n;
  int i;
  float *tm;
  if (vStride==0) vStride=3*sizeof(float);
  if (nStride==0) nStride=3*sizeof(float);

  if (vertex)
  {
    v=(BYTE *)vertex;
    for (i=0; i<vCount; i++, v+=vStride)
    {
      glmCopyVector3f((float *)v, temp);
	  glmMultVectorByMatrix43f((float *)v, temp, m);
    }
  }

  if (normal)
  {
    n=(BYTE *)normal;
    tm=(float *)malloc(16*sizeof(float));
	glmGenerateInverseMatrix4f(tm, m);
    glmTransposeMatrix4f(tm);
	for (i=0; i<vCount; i++, n+=nStride)
    {
      glmCopyVector3f((float *)n, temp);
	  glmMultVectorByMatrix43f((float *)n, temp, tm);
    }
    free(tm);
  }
}

void glmTransformMeshd(double *vertex, DWORD vStride,
                       double *normal, DWORD nStride,
                       DWORD vCount, const double m[16])
{
  double temp[3];
  BYTE *v, *n;
  int i;
  double *tm;
  if (vStride==0) vStride=3*sizeof(double);
  if (nStride==0) nStride=3*sizeof(double);


  if (vertex)
  {
    v=(BYTE *)vertex;
    for (i=0; i<vCount; i++, v+=vStride)
    {
      glmCopyVector3d((double *)v, temp);
	  glmMultVectorByMatrix43d((double *)v, temp, m);
    }
  }

  if (normal)
  {
    n=(BYTE *)normal;
    tm=(double *)malloc(16*sizeof(double));
	glmGenerateInverseMatrix4d(tm, m);
    glmTransposeMatrix4d(tm);
	for (i=0; i<vCount; i++, n+=nStride)
    {
      glmCopyVector3d((double *)n, temp);
	  glmMultVectorByMatrix43d((double *)n, temp, tm);
    }
    free(tm);
  }
}

void glmGenerateTransformedMeshf(float *resultVertex, DWORD rvStride,
                       const float *vertex, DWORD vStride,
                       float *resultNormal, DWORD rnStride,
					   const float *normal, DWORD nStride,
                       DWORD vCount, const float m[16])
{
  BYTE *v, *n, *rv, *rn;
  int i;
  float *tm;
  if (vStride==0) vStride=3*sizeof(float);
  if (nStride==0) nStride=3*sizeof(float);
  if (rvStride==0) rvStride=3*sizeof(float);
  if (rnStride==0) rnStride=3*sizeof(float);

  if (vertex&&resultVertex)
  {
    v=(BYTE *)vertex;
    rv=(BYTE *)resultVertex;
    for (i=0; i<vCount; i++, v+=vStride, rv+=rvStride)
    {
      glmMultVectorByMatrix43f((float *)rv, (float *)v, m);
    }
  }

  if (normal&&resultNormal)
  {
    n=(BYTE *)normal;
    rn=(BYTE *)resultNormal;
    tm=(float *)malloc(16*sizeof(float));
	glmGenerateInverseMatrix4f(tm, m);
    glmTransposeMatrix4f(tm);
    for (i=0; i<vCount; i++, n+=nStride, rn+=rnStride)
    {
      glmMultVectorByMatrix43f((float *)rn, (float *)n, tm);
    }
    free(tm);
  }
}

void glmGenerateTransformedMeshd(double *resultVertex, DWORD rvStride,
                       const double *vertex, DWORD vStride,
                       double *resultNormal, DWORD rnStride,
					   const double *normal, DWORD nStride,
                       DWORD vCount, const double m[16])
{
  BYTE *v, *n, *rv, *rn;
  int i;
  double *tm;
  if (vStride==0) vStride=3*sizeof(double);
  if (nStride==0) nStride=3*sizeof(double);
  if (rvStride==0) rvStride=3*sizeof(double);
  if (rnStride==0) rnStride=3*sizeof(double);

  if (vertex&&resultVertex)
  {
    v=(BYTE *)vertex;
    rv=(BYTE *)resultVertex;
    for (i=0; i<vCount; i++, v+=vStride, rv+=rvStride)
    {
      glmMultVectorByMatrix43d((double *)rv, (double *)v, m);
    }
  }

  if (normal&&resultNormal)
  {
    n=(BYTE *)normal;
    rn=(BYTE *)resultNormal;
    tm=(double *)malloc(16*sizeof(double));
    glmGenerateInverseMatrix4d(tm, m);
    glmTransposeMatrix4d(tm);
    for (i=0; i<vCount; i++, n+=nStride, rn+=rnStride)
    {
      glmMultVectorByMatrix43d((double *)rn, (double *)n, tm);
    }
    free(tm);
  }
}

void glmGenerateShadowMatrix4f(float m[16], const float planeOrigin[3],
                               const float planeNormal[3], const float light[3])
{
  float d=-glmDotProduct3f(planeOrigin, planeNormal);
  float p=light[0]*planeNormal[0]+light[1]*planeNormal[1]+light[2]*planeNormal[2]+d;

  m[0]=p-light[0]*planeNormal[0];
  m[4]= -light[0]*planeNormal[1];
  m[8]= -light[0]*planeNormal[2];
  m[12]= -light[0]*d;

  m[1]= -light[1]*planeNormal[0];
  m[5]=p-light[1]*planeNormal[1];
  m[9]= -light[1]*planeNormal[2];
  m[13]= -light[1]*d;

  m[2]=  -light[2]*planeNormal[0];
  m[6]=  -light[2]*planeNormal[1];
  m[10]=p-light[2]*planeNormal[2];
  m[14]= -light[2]*d;

  m[3]= -planeNormal[0];
  m[7]= -planeNormal[1];
  m[11]= -planeNormal[2];
  m[15]=p-d;
}

void glmGenerateShadowMatrix4d(double m[16], const double planeOrigin[3],
                               const double planeNormal[3], const double light[3])
{
  double d=-glmDotProduct3d(planeOrigin, planeNormal);
  double p=light[0]*planeNormal[0]+light[1]*planeNormal[1]+light[2]*planeNormal[2]+d;

  m[0]=p-light[0]*planeNormal[0];
  m[4]= -light[0]*planeNormal[1];
  m[8]= -light[0]*planeNormal[2];
  m[12]= -light[0]*d;

  m[1]= -light[1]*planeNormal[0];
  m[5]=p-light[1]*planeNormal[1];
  m[9]= -light[1]*planeNormal[2];
  m[13]= -light[1]*d;

  m[2]=  -light[2]*planeNormal[0];
  m[6]=  -light[2]*planeNormal[1];
  m[10]=p-light[2]*planeNormal[2];
  m[14]= -light[2]*d;

  m[3]= -planeNormal[0];
  m[7]= -planeNormal[1];
  m[11]= -planeNormal[2];
  m[15]=p-d;
}

void glmGeneratePerspectiveProjectionMatrixf(float m[16], float left, float right,
                                                          float bottom, float top,
                                                          float near, float far)
{
  m[0]=(2*near)/(right-left);
  m[1]=0;
  m[2]=0;
  m[3]=0;
  m[4]=0;
  m[5]=(2*near)/(top-bottom);
  m[6]=0;
  m[7]=0;
  m[8]=(right+left)/(right-left);
  m[9]=(top+bottom)/(top-bottom);
  m[10]=-(far+near)/(far-near);
  m[11]=-1;
  m[12]=0;
  m[13]=0;
  m[14]=(-2*far*near)/(far-near);
  m[15]=0;
}

void glmGeneratePerspectiveProjectionMatrixd(double m[16], double left, double right,
                                                           double bottom, double top,
                                                           double near, double far)
{
  m[0]=(2*near)/(right-left);
  m[1]=0;
  m[2]=0;
  m[3]=0;
  m[4]=0;
  m[5]=(2*near)/(top-bottom);
  m[6]=0;
  m[7]=0;
  m[8]=(right+left)/(right-left);
  m[9]=(top+bottom)/(top-bottom);
  m[10]=-(far+near)/(far-near);
  m[11]=-1;
  m[12]=0;
  m[13]=0;
  m[14]=(-2*far*near)/(far-near);
  m[15]=0;
}

