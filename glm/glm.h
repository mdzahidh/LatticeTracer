/* lcc generated prototype file. Do not edit */
#ifndef __glm_h__
#define __glm_h__

#define GLM_VERSION 1.3
typedef int BOOL;
typedef unsigned char BYTE;
typedef unsigned long int DWORD;
#define TRUE 1
#define FALSE 0
#define GLM_2PI 6.283185307179586476925286766559
#define GLM_PI 3.1415926535897932384626433832795

float glmRandomNumberf(float, float);
double glmRandomNumberd(double, double);
float glmRadf(float);
double glmRadd(double);
float glmDegf(float);
double glmDegd(double);
float glmMax2f(float, float);
double glmMax2d(double, double);
float glmMax3f(float, float, float);
double glmMax3d(double, double, double);
float glmMin2f(float, float);
double glmMin2d(double, double);
float glmMin3f(float, float, float);
double glmMin3d(double, double, double);
int glmMaxIndex2f(float, float, int, int);
int glmMaxIndex2d(double, double, int, int);
int glmMaxIndex3f(float, float, float, int, int, int);
int glmMaxIndex3d(double, double, double, int, int, int);
int glmMinIndex2f(float, float, int, int);
int glmMinIndex2d(double, double, int, int);
int glmMinIndex3f(float, float, float, int, int, int);
int glmMinIndex3d(double, double, double, int, int, int);
int glmIsInListi(unsigned long, unsigned long *, unsigned long);
int glmIsInListf(float, float *, unsigned long);
int glmIsInListd(double, double *, unsigned long);
#endif
