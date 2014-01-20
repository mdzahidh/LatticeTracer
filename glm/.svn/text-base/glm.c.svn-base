// GLM - a math library for OpenGL
// Programmed by Markus Ilmola
// markus.ilmola@pp.inet.fi
// http://markus_ilmola.tripod.com
// http://personal.inet.fi/koti/markus.ilmola

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#define GLM_VERSION 1.3
typedef int BOOL;
typedef unsigned char BYTE;
typedef unsigned long int DWORD;
#define TRUE 1
#define FALSE 0
#define GLM_2PI 6.283185307179586476925286766559
#define GLM_PI 3.1415926535897932384626433832795


float glmRandomNumberf(float min, float max)
{
  return min+(max-min)*rand()/(float)RAND_MAX;
}

double glmRandomNumberd(double min, double max)
{
  return min+(max-min)*rand()/(double)RAND_MAX;
}

float glmRadf(float d)
{
  return d*0.017453292519943295769236907684886f;
}

double glmRadd(double d)
{
  return d*0.017453292519943295769236907684886;
}

float glmDegf(float r)
{
  return r*57.295779513082320876798154814105f;
}

double glmDegd(double r)
{
  return r*57.295779513082320876798154814105;
}

float glmMax2f(float a, float b)
{
  if (a>b) return a;
  return b;
}

double glmMax2d(double a, double b)
{
  if (a>b) return a;
  return b;
}

float glmMax3f(float a, float b, float c)
{
  if (a>b)
  {
    if (a>c) return a;
    else return c;
  }
  else
  {
    if (b>c) return b;
    else return c;
  }
}

double glmMax3d(double a, double b, double c)
{
  if (a>b)
  {
    if (a>c) return a;
    else return c;
  }
  else
  {
    if (b>c) return b;
    else return c;
  }
}

float glmMin2f(float a, float b)
{
  if (a<b) return a;
  return b;
}

double glmMin2d(double a, double b)
{
  if (a<b) return a;
  return b;
}

float glmMin3f(float a, float b, float c)
{
  if (a<b)
  {
    if (a<c) return a;
    else return c;
  }
  else
  {
    if (b<c) return b;
    else return c;
  }
}

double glmMin3d(double a, double b, double c)
{
  if (a<b)
  {
    if (a<c) return a;
    else return c;
  }
  else
  {
    if (b<c) return b;
    else return c;
  }
}

int glmMaxIndex2f(float a, float b, int i1, int i2)
{
  if (a>b) return i1;
  return i2;
}

int glmMaxIndex2d(double a, double b, int i1, int i2)
{
  if (a>b) return i1;
  return i2;
}

int glmMaxIndex3f(float a, float b, float c, int i1, int i2, int i3)
{
  if (a>b)
  {
    if (a>c) return i1;
    else return i3;
  }
  else
  {
    if (b>c) return i2;
    else return i3;
  }
}

int glmMaxIndex3d(double a, double b, double c, int i1, int i2, int i3)
{
  if (a>b)
  {
    if (a>c) return i1;
    else return i3;
  }
  else
  {
    if (b>c) return i2;
    else return i3;
  }
}

int glmMinIndex2f(float a, float b, int i1, int i2)
{
  if (a<b) return i1;
  return i2;
}

int glmMinIndex2d(double a, double b, int i1, int i2)
{
  if (a<b) return i1;
  return i2;
}

int glmMinIndex3f(float a, float b, float c, int i1, int i2, int i3)
{
  if (a<b)
  {
    if (a<c) return i1;
    else return i3;
  }
  else
  {
    if (b<c) return i2;
    else return i3;
  }
}

int glmMinIndex3d(double a, double b, double c, int i1, int i2, int i3)
{
  if (a<b)
  {
    if (a<c) return i1;
    else return i3;
  }
  else
  {
    if (b<c) return i2;
    else return i3;
  }
}

BOOL glmIsInListi(DWORD a, DWORD *list, DWORD iCount)
{
  DWORD i;
  for (i=0; i<iCount; i++)
    if (list[i]==a) return TRUE;
  return FALSE;
}

BOOL glmIsInListf(float a, float *list, DWORD iCount)
{
  DWORD i;
  for (i=0; i<iCount; i++)
    if (list[i]==a) return TRUE;
  return FALSE;
}

BOOL glmIsInListd(double a, double *list, DWORD iCount)
{
  DWORD i;
  for (i=0; i<iCount; i++)
    if (list[i]==a) return TRUE;
  return FALSE;
}
