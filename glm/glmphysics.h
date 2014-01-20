/* lcc generated prototype file. Do not edit */
#ifndef __glmphysics_h__
#define __glmphysics_h__
void glmSolveOdeEulerf(float *, const float *, unsigned long, float, void (*)(float *, const float *, float), float *);
void glmSolveOdeEulerd(double *, const double *, unsigned long, double, void (*)(double *, const double *, double), double *);
void glmSolveOdeMidpointf(float *, const float *, unsigned long, float, void (*)(float *, const float *, float), float *);
void glmSolveOdeMidpointd(double *, const double *, unsigned long, double, void (*)(double *, const double *, double), double *);
void glmSolveOdeRungeKuttaf(float *, const float *, unsigned long, float, void (*)(float *, const float *, float), float *);
float glmComputeImpulse2f(const float *, const float *, float, float, const float *, float);
double glmComputeImpulse2d(const double *, const double *, float, float, const double *, double);
float glmComputeImpulse3f(const float *, const float *, float, float, const float *, float);
float glmComputeImpulse3d(const double *, const double *, double, double, const double *, double);
void glmGenerateSkewSymmetricMatrix3f(float *, const float *);
void glmGenerateSkewSymmetricMatrix3d(double *, const double *);
#endif
