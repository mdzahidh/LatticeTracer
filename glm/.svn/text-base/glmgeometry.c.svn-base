// GLM - a math library for OpenGL
// Programmed by Markus Ilmola
// markus.ilmola@pp.inet.fi
// http://markus_ilmola.tripod.com
// http://personal.inet.fi/koti/markus.ilmola

// Geometry stuff

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "glm.h"
#include "glmvector.h"
#include "glmmatrix.h"

// Tests if a ray intersects a line.
// Returns distance to the intersection point or -1 if no intersection.
// Intersection point is stored to "intersection" if it is not NULL.
// "lineNormal" and "rayDirection" must be normalized. (lenght=1)
float glmRayLineIntersectionf(const float rayOrigin[2], const float rayDirection[2],
                              const float lineOrigin[2], const float lineNormal[2],
                              float intersection[2])
{
  // Locals.
  float numer, denom, d;

  // Calculate the distance to the intersection point.
  d=glmDotProduct2f(lineNormal, lineOrigin);
  numer=-glmDotProduct2f(lineNormal, rayOrigin)+d;
  denom=glmDotProduct2f(lineNormal, rayDirection);
  if (denom==0) return -1;
  d=numer/denom;
  if (d<0) return -1;

  // Store the intersection point.
  if (intersection)
  {
    intersection[0]=rayOrigin[0]+d*rayDirection[0];
    intersection[1]=rayOrigin[1]+d*rayDirection[1];
  }

  // Return the distance.
  return d;
}

double glmRayLineIntersectiond(const double rayOrigin[2], const double rayDirection[2],
                            const double lineOrigin[2], const double lineNormal[2],
                            double intersection[2])
{
  // Locals.
  double numer, denom, d;

  // Calculate the distance to the intersection point.
  d=glmDotProduct2d(lineNormal, lineOrigin);
  numer=-glmDotProduct2d(lineNormal, rayOrigin)+d;
  denom=glmDotProduct2d(lineNormal, rayDirection);
  if (denom==0) return -1;
  d=numer/denom;
  if (d<0) return -1;

  // Store the intersection point.
  if (intersection)
  {
    intersection[0]=rayOrigin[0]+d*rayDirection[0];
    intersection[1]=rayOrigin[1]+d*rayDirection[1];
  }

  // Return the distance.
  return d;
}

// Tests if a ray intersects a plane.
// Returns distance to the intersection point or -1 if no intersection.
// Intersection point is stored to "intersection" if it is not NULL.
// "planeNormal" and "rayDirection" must be normalized. (lenght=1)
float glmRayPlaneIntersectionf(const float rayOrigin[3], const float rayDirection[3],
      const float planeOrigin[3], const float planeNormal[3], float intersection[3])
{
  // Locals.
  float numer, denom, d;

  // Calculate the distance to the intersection point.
  d=glmDotProduct3f(planeNormal, planeOrigin);
  numer=-glmDotProduct3f(planeNormal, rayOrigin)+d;
  denom=glmDotProduct3f(planeNormal, rayDirection);
  if (denom==0) return -1;
  d=numer/denom;
  if (d<0) return -1;

  // Store the intersection point.
  if (intersection)
  {
    intersection[0]=rayOrigin[0]+d*rayDirection[0];
    intersection[1]=rayOrigin[1]+d*rayDirection[1];
    intersection[2]=rayOrigin[2]+d*rayDirection[2];
  }

  // Return the distance.
  return d;
}

double glmRayPlaneIntersectiond(const double rayOrigin[3], const double rayDirection[3],
       const double planeOrigin[3], const double planeNormal[3], double intersection[3])
{
  // Locals.
  double numer, denom, d;

  // Calculate the distance to the intersection point.
  d=glmDotProduct3d(planeNormal, planeOrigin);
  numer=-glmDotProduct3d(planeNormal, rayOrigin)+d;
  denom=glmDotProduct3d(planeNormal, rayDirection);
  if (denom==0) return -1;
  d=numer/denom;
  if (d<0) return -1;

  // Store the intersection point.
  if (intersection)
  {
    intersection[0]=rayOrigin[0]+d*rayDirection[0];
    intersection[1]=rayOrigin[1]+d*rayDirection[1];
    intersection[2]=rayOrigin[2]+d*rayDirection[2];
  }

  // Return the distance.
  return d;
}

float glmRaySpaceIntersectionf(const float rayOrigin[4], const float rayDirection[4],
      const float spaceOrigin[4], const float spaceNormal[4], float intersection[4])
{
  // Locals.
  float numer, denom, d;

  // Calculate the distance to the intersection point.
  d=glmDotProduct4f(spaceNormal, spaceOrigin);
  numer=-glmDotProduct4f(spaceNormal, rayOrigin)+d;
  denom=glmDotProduct4f(spaceNormal, rayDirection);
  if (denom==0) return -1;
  d=numer/denom;
  if (d<0) return -1;

  // Store the intersection point.
  if (intersection)
  {
    intersection[0]=rayOrigin[0]+d*rayDirection[0];
    intersection[1]=rayOrigin[1]+d*rayDirection[1];
    intersection[2]=rayOrigin[2]+d*rayDirection[2];
    intersection[3]=rayOrigin[3]+d*rayDirection[3];
  }

  // Return the distance.
  return d;
}

double glmRaySpaceIntersectiond(const double rayOrigin[4], const double rayDirection[4],
      const double spaceOrigin[4], const double spaceNormal[4], double intersection[4])
{
  // Locals.
  double numer, denom, d;

  // Calculate the distance to the intersection point.
  d=glmDotProduct4d(spaceNormal, spaceOrigin);
  numer=-glmDotProduct4d(spaceNormal, rayOrigin)+d;
  denom=glmDotProduct4d(spaceNormal, rayDirection);
  if (denom==0) return -1;
  d=numer/denom;
  if (d<0) return -1;

  // Store the intersection point.
  if (intersection)
  {
    intersection[0]=rayOrigin[0]+d*rayDirection[0];
    intersection[1]=rayOrigin[1]+d*rayDirection[1];
    intersection[2]=rayOrigin[2]+d*rayDirection[2];
    intersection[3]=rayOrigin[3]+d*rayDirection[3];
  }

  // Return the distance.
  return d;
}

float glmRayHyperSpaceIntersectionf(const float *rayOrigin, const float *rayDirection,
      const float *spaceOrigin, const float *spaceNormal, float *intersection, DWORD size)
{
  // Locals.
  float numer, denom, d;
  DWORD i;

  // Calculate the distance to the intersection point.
  d=glmDotProductf(spaceNormal, spaceOrigin, size);
  numer=-glmDotProductf(spaceNormal, rayOrigin, size)+d;
  denom=glmDotProductf(spaceNormal, rayDirection, size);
  if (denom==0) return -1;
  d=numer/denom;
  if (d<0) return -1;

  // Store the intersection point.
  if (intersection)
  {
    for (i=0; i<size; i++) intersection[i]=rayOrigin[i]+d*rayDirection[i];
  }

  // Return the distance.
  return d;
}

double glmRayHyperSpaceIntersectiond(const double *rayOrigin, const double *rayDirection,
      const double *spaceOrigin, const double *spaceNormal, double *intersection, DWORD size)
{
  // Locals.
  double numer, denom, d;
  DWORD i;

  // Calculate the distance to the intersection point.
  d=glmDotProductd(spaceNormal, spaceOrigin, size);
  numer=-glmDotProductd(spaceNormal, rayOrigin, size)+d;
  denom=glmDotProductd(spaceNormal, rayDirection, size);
  if (denom==0) return -1;
  d=numer/denom;
  if (d<0) return -1;

  // Store the intersection point.
  if (intersection)
  {
    for (i=0; i<size; i++) intersection[i]=rayOrigin[i]+d*rayDirection[i];
  }

  // Return the distance.
  return d;
}


BOOL glmLineSegmentPlaneIntersectionf(const float p1[3], const float p2[3],
      const float planeOrigin[3], const float planeNormal[3], float intersection[3], float *a)
{
  // Locals.
  float numer, denom, d;
  float v[3];
  float l;

  // Make the lines direction vector
  glmMakeVector3f(v, p1, p2);
  l=glmNormalizeVector3f(v);

  // Calculate the distance to the intersection point.
  d=glmDotProduct3f(planeNormal, planeOrigin);
  numer=-glmDotProduct3f(planeNormal, p1)+d;
  denom=glmDotProduct3f(planeNormal, v);
  if (denom==0) return FALSE;
  d=numer/denom;
  if (d<0) return FALSE;
  if (d>l) return FALSE;

  // Store the intersection point.
  if (intersection)
  {
    intersection[0]=p1[0]+d*v[0];
    intersection[1]=p1[1]+d*v[1];
    intersection[2]=p1[2]+d*v[2];
  }

  if (a) *a=d/l;

  return TRUE;
}

BOOL glmLineSegmentPlaneIntersectiond(const double p1[3], const double p2[3],
      const double planeOrigin[3], const double planeNormal[3], double intersection[3], double *a)
{
  // Locals.
  double numer, denom, d;
  double v[3];
  double l;

  // Make the lines direction vector
  glmMakeVector3d(v, p1, p2);
  l=glmNormalizeVector3d(v);

  // Calculate the distance to the intersection point.
  d=glmDotProduct3d(planeNormal, planeOrigin);
  numer=-glmDotProduct3d(planeNormal, p1)+d;
  denom=glmDotProduct3d(planeNormal, v);
  if (denom==0) return FALSE;
  d=numer/denom;
  if (d<0) return FALSE;
  if (d>l) return FALSE;

  // Store the intersection point.
  if (intersection)
  {
    intersection[0]=p1[0]+d*v[0];
    intersection[1]=p1[1]+d*v[1];
    intersection[2]=p1[2]+d*v[2];
  }

  if (a) *a=d/l;

  return TRUE;
}

BOOL glmLineSegmentSpaceIntersectionf(const float p1[4], const float p2[4],
      const float spaceOrigin[4], const float spaceNormal[4], float intersection[4], float *a)
{
  // Locals.
  float numer, denom, d;
  float v[4];
  float l;

  // Make the lines direction vector
  glmMakeVector4f(v, p1, p2);
  l=glmNormalizeVector4f(v);

  // Calculate the distance to the intersection point.
  d=glmDotProduct4f(spaceNormal, spaceOrigin);
  numer=-glmDotProduct4f(spaceNormal, p1)+d;
  denom=glmDotProduct4f(spaceNormal, v);
  if (denom==0) return FALSE;
  d=numer/denom;
  if (d<0) return FALSE;
  if (d>l) return FALSE;

  // Store the intersection point.
  if (intersection)
  {
    intersection[0]=p1[0]+d*v[0];
    intersection[1]=p1[1]+d*v[1];
    intersection[2]=p1[2]+d*v[2];
    intersection[3]=p1[3]+d*v[3];
  }

  if (a) *a=d/l;

  return TRUE;
}

BOOL glmLineSegmentSpaceIntersectiond(const double p1[4], const double p2[4],
      const double spaceOrigin[4], const double spaceNormal[4], double intersection[4], double *a)
{
  // Locals.
  double numer, denom, d;
  double v[4];
  double l;

  // Make the lines direction vector
  glmMakeVector4d(v, p1, p2);
  l=glmNormalizeVector4d(v);

  // Calculate the distance to the intersection point.
  d=glmDotProduct4d(spaceNormal, spaceOrigin);
  numer=-glmDotProduct4d(spaceNormal, p1)+d;
  denom=glmDotProduct4d(spaceNormal, v);
  if (denom==0) return FALSE;
  d=numer/denom;
  if (d<0) return FALSE;
  if (d>l) return FALSE;

  // Store the intersection point.
  if (intersection)
  {
    intersection[0]=p1[0]+d*v[0];
    intersection[1]=p1[1]+d*v[1];
    intersection[2]=p1[2]+d*v[2];
    intersection[3]=p1[3]+d*v[3];
  }

  if (a) *a=d/l;

  return TRUE;
}


// Tests if a point lies inside a trianle.
// If it does TRUE is returned.
// Otherwise FALSE is returned.
// Point "p" must lie on the same plane as the trianle.
BOOL glmIsPointInTriangle2f(const float p[2], const float c1[2], const float c2[2], const float c3[2])
{
  // Locals
  float v1[3], v2[3], v3[3];
  float angle;

  // Create tree vectors from the point to all of the corners.
  glmMakeVector2f(v1, p, c1);
  glmMakeVector2f(v2, p, c2);
  glmMakeVector2f(v3, p, c3);
  glmNormalizeVector2f(v1);
  glmNormalizeVector2f(v2);
  glmNormalizeVector2f(v3);

  // Sum of the angles should be 2*PI.
  angle=0;
  angle+=glmComputeNormalizedVectorAngle2fRad(v1, v2);
  angle+=glmComputeNormalizedVectorAngle2fRad(v2, v3);
  angle+=glmComputeNormalizedVectorAngle2fRad(v3, v1);
  if (angle>6.282f && angle<6.284f) return TRUE;
  else return FALSE;
}

BOOL glmIsPointInTriangle2d(const double p[2], const double c1[2], const double c2[2], const double c3[2])
{
  // Locals
  double v1[3], v2[3], v3[3];
  double angle;

  // Create tree vectors from the point to all of the corners.
  glmMakeVector2d(v1, p, c1);
  glmMakeVector2d(v2, p, c2);
  glmMakeVector2d(v3, p, c3);
  glmNormalizeVector2d(v1);
  glmNormalizeVector2d(v2);
  glmNormalizeVector2d(v3);

  // Sum of the angles should be 2*PI.
  angle=0;
  angle+=glmComputeNormalizedVectorAngle2dRad(v1, v2);
  angle+=glmComputeNormalizedVectorAngle2dRad(v2, v3);
  angle+=glmComputeNormalizedVectorAngle2dRad(v3, v1);
  if (angle>6.282 && angle<6.284) return TRUE;
  else return FALSE;
}

BOOL glmIsPointInTriangle3f(const float p[3], const float c1[3], const float c2[3], const float c3[3])
{
  // Locals
  float v1[3], v2[3], v3[3];
  float angle;

  // Create tree vectors from the point to all of the corners.
  glmMakeVector3f(v1, p, c1);
  glmMakeVector3f(v2, p, c2);
  glmMakeVector3f(v3, p, c3);
  glmNormalizeVector3f(v1);
  glmNormalizeVector3f(v2);
  glmNormalizeVector3f(v3);

  // Sum of the angles should be 2*PI.
  angle=0;
  angle+=glmComputeNormalizedVectorAngle3fRad(v1, v2);
  angle+=glmComputeNormalizedVectorAngle3fRad(v2, v3);
  angle+=glmComputeNormalizedVectorAngle3fRad(v3, v1);
  if (angle>6.282f && angle<6.284f) return TRUE;
  else return FALSE;
}

BOOL glmIsPointInTriangle3d(const double p[3], const double c1[3], const double c2[3], const double c3[3])
{
  // Locals
  double v1[3], v2[3], v3[3];
  double angle;

  // Create tree vectors from the point to all of the corners.
  glmMakeVector3d(v1, p, c1);
  glmMakeVector3d(v2, p, c2);
  glmMakeVector3d(v3, p, c3);
  glmNormalizeVector3d(v1);
  glmNormalizeVector3d(v2);
  glmNormalizeVector3d(v3);

  // Sum of the angles should be 2*PI.
  angle=0;
  angle+=glmComputeNormalizedVectorAngle3dRad(v1, v2);
  angle+=glmComputeNormalizedVectorAngle3dRad(v2, v3);
  angle+=glmComputeNormalizedVectorAngle3dRad(v3, v1);
  if (angle>6.282 && angle<6.284) return TRUE;
  else return FALSE;
}


// Tests if a ray intersects a trianle.
// Returns the distence to the intersection point or -1 if no intersection.
// The intersection point is stored to "intersection".
float glmRayTriangleIntersectionf(const float rayOrigin[3], const float rayVector[3],
              const float c1[3], const float c2[3], const float c3[3], const float normal[3],
							  float intersection[3])
{
  // Locals
  float distance;

  // Calculete the intersection
  distance=glmRayPlaneIntersectionf(rayOrigin, rayVector, c1, normal, intersection);
  if (distance<0) return -1;
  if (!glmIsPointInTriangle3f(intersection, c1, c2, c3)) return -1;

  return distance;
}

double glmRayTriangleIntersectiond(const double rayOrigin[3], const double rayVector[3],
              const double c1[3], const double c2[3], const double c3[3], const double normal[3],
							  double intersection[3])
{
  // Locals
  double distance;

  // Calculete the intersection
  distance=glmRayPlaneIntersectiond(rayOrigin, rayVector, c1, normal, intersection);
  if (distance<0) return -1;
  if (!glmIsPointInTriangle3d(intersection, c1, c2, c3)) return -1;

  return distance;
}

void glmClosestPointOnLineSegment2f(float result[2], const float a[2], const float b[2], const float p[2])
{
  // Locals
  float c[2], v[2];
  float d, t;

  // Determine the length of the vector from a to b
  glmMakeVector2f(c, a, p);
  glmMakeVector2f(v, a, b);
  d=glmNormalizeVector2f(v);
  t=glmDotProduct2f(v, c);

  // Check to see if ‘t’ is beyond the extents of the line segment
  if (t < 0)
  {
    glmCopyVector2f(a, result);
    return;
  }
  if (t > d)
  {
    glmCopyVector2f(b, result);
    return;
  }

  // Return the point between ‘a’ and ‘b’
  glmMultVectorByScalar2f(v, t);
  glmAddVector2f(result, a, v);
}

void glmClosestPointOnLineSegment2d(double result[2], const double a[2], const double b[2], const double p[2])
{
  // Locals
  double c[2], v[2];
  double d, t;

  // Determine the length of the vector from a to b
  glmMakeVector2d(c, a, p);
  glmMakeVector2d(v, a, b);
  d=glmNormalizeVector2d(v);
  t=glmDotProduct2d(v, c);

  // Check to see if ‘t’ is beyond the extents of the line segment
  if (t < 0)
  {
    glmCopyVector2d(a, result);
    return;
  }
  if (t > d)
  {
    glmCopyVector2d(b, result);
    return;
  }

  // Return the point between ‘a’ and ‘b’
  glmMultVectorByScalar2d(v, t);
  glmAddVector2d(result, a, v);
}



void glmClosestPointOnLineSegment3f(float result[3], const float a[3], const float b[3], const float p[3])
{
  // Locals
  float c[3], v[3];
  float d, t;

  // Determine the length of the vector from a to b
  glmMakeVector3f(c, a, p);
  glmMakeVector3f(v, a, b);
  d=glmNormalizeVector3f(v);
  t=glmDotProduct3f(v, c);

  // Check to see if ‘t’ is beyond the extents of the line segment
  if (t < 0)
  {
    glmCopyVector3f(a, result);
    return;
  }
  if (t > d)
  {
    glmCopyVector3f(b, result);
    return;
  }

  // Return the point between ‘a’ and ‘b’
  glmMultVectorByScalar3f(v, t);
  glmAddVector3f(result, a, v);
}

void glmClosestPointOnLineSegment3d(double result[3], const double a[3], const double b[3], const double p[3])
{
  // Locals
  double c[3], v[3];
  double d, t;

  // Determine the length of the vector from a to b
  glmMakeVector3d(c, a, p);
  glmMakeVector3d(v, a, b);
  d=glmNormalizeVector3d(v);
  t=glmDotProduct3d(v, c);

  // Check to see if ‘t’ is beyond the extents of the line segment
  if (t < 0)
  {
    glmCopyVector3d(a, result);
    return;
  }
  if (t > d)
  {
    glmCopyVector3d(b, result);
    return;
  }

  // Return the point between ‘a’ and ‘b’
  glmMultVectorByScalar3d(v, t);
  glmAddVector3d(result, a, v);
}


float glmRayCircleIntersectionf(const float rO[2], const float rV[2],
                                const float sO[2], float sR, float intersection[2])
{
  float c, v, d;
  float q[2];

  glmMakeVector2f(q, rO, sO);
  c=glmComputeVectorLength2f(q);
  if (c<sR) return 0;
  v=glmDotProduct2f(q, rV);
  d=sR*sR-(c*c-v*v);

  // If there was no intersection, return -1
  if (d<0) return -1;
  d=v-sqrt(d);

  if (intersection)
  {
    intersection[0]=rO[0]+d*rV[0];
    intersection[1]=rO[1]+d*rV[1];
  }

  return d;
}

double glmRayCircleIntersectiond(const double rO[2], const double rV[2],
                               const double sO[2], double sR, double intersection[2])
{
  double c, v, d;
  double q[2];

  glmMakeVector2d(q, rO, sO);
  c=glmComputeVectorLength2d(q);
  if (c<sR) return 0;
  v=glmDotProduct2d(q, rV);
  d=sR*sR-(c*c-v*v);

  // If there was no intersection, return -1
  if (d<0) return -1;
  d=v-sqrt(d);

  if (intersection)
  {
    intersection[0]=rO[0]+d*rV[0];
    intersection[1]=rO[1]+d*rV[1];
  }

  return d;
}

float glmRaySphereIntersectionf(const float rO[3], const float rV[3],
                              const float sO[3], float sR, float intersection[3])
{
  float c, v, d;
  float q[3];

  glmMakeVector3f(q, rO, sO);
  c=glmComputeVectorLength3f(q);
  if (c<sR) return 0;
  v=glmDotProduct3f(q, rV);
  d=sR*sR-(c*c-v*v);

  // If there was no intersection, return -1
  if (d<0) return -1;
  d=v-sqrt(d);

  if (intersection)
  {
    intersection[0]=rO[0]+d*rV[0];
    intersection[1]=rO[1]+d*rV[1];
    intersection[2]=rO[2]+d*rV[2];
  }

  return d;
}

double glmRaySphereIntersectiond(const double rO[3], const double rV[3],
                               const double sO[3], double sR, double intersection[3])
{
  double c, v, d;
  double q[3];

  glmMakeVector3d(q, rO, sO);
  c=glmComputeVectorLength3d(q);
  if (c<sR) return 0;
  v=glmDotProduct3d(q, rV);
  d=sR*sR-(c*c-v*v);

  // If there was no intersection, return -1
  if (d<0) return -1;
  d=v-sqrt(d);

  if (intersection)
  {
    intersection[0]=rO[0]+d*rV[0];
    intersection[1]=rO[1]+d*rV[1];
    intersection[2]=rO[2]+d*rV[2];
  }

  return d;
}


float glmRayHyperSphereIntersectionf(const float *rO, const float *rV,
                              const float *sO, float sR, float *intersection, DWORD size)
{
  float c, v, d;
  float *q;
  DWORD i;

  q=(float *)malloc(size*sizeof(float));

  glmMakeVectorf(q, rO, sO, size);
  c=glmComputeVectorLengthf(q, size);
  if (c<sR) return 0;
  v=glmDotProductf(q, rV, size);
  d=sR*sR-(c*c-v*v);

  // If there was no intersection, return -1
  if (d<0) return -1;
  d=v-sqrt(d);

  if (intersection)
  {
    for (i=0; i<size; i++) intersection[i]=rO[i]+d*rV[i];
  }

  free(q);
  return d;
}

double glmRayHyperSphereIntersectiond(const double *rO, const double *rV,
                              const double *sO, double sR, double *intersection, DWORD size)
{
  double c, v, d;
  double *q;
  DWORD i;

  q=(double *)malloc(size*sizeof(double));

  glmMakeVectord(q, rO, sO, size);
  c=glmComputeVectorLengthd(q, size);
  if (c<sR) return 0;
  v=glmDotProductd(q, rV, size);
  d=sR*sR-(c*c-v*v);

  // If there was no intersection, return -1
  if (d<0) return -1;
  d=v-sqrt(d);

  if (intersection)
  {
    for (i=0; i<size; i++) intersection[i]=rO[i]+d*rV[i];
  }

  free(q);
  return d;
}

float glmRayMeshIntersectionf(const float rO[3], const float rV[3],
                           const float *vertex, DWORD stride,
                           const DWORD *indicates, DWORD iCount, float intersection[3], float normal[3])
{
  float distance=-1, tempDistance;
  BYTE *v;
  float tempNormal[3], tempIntersection[3];
  float *vert1, *vert2, *vert3;
  int i;

  if (stride==0) stride=3*sizeof(float);

  v=(BYTE *)vertex;
  for (i=0; i<iCount; i+=3)
  {
    vert1=(float *)(v+indicates[i]*stride);
    vert2=(float *)(v+indicates[i+1]*stride);
    vert3=(float *)(v+indicates[i+2]*stride);
    glmComputeTriangleNormalf(tempNormal, vert1, vert2, vert3);
    tempDistance=glmRayTriangleIntersectionf(rO, rV, vert1, vert2, vert3, tempNormal, tempIntersection);
    if (tempDistance>=0&&(tempDistance<distance||distance<0))
    {
      distance=tempDistance;
      if (intersection)
      {
        glmCopyVector3f(tempIntersection, intersection);
      }
      if (normal)
      {
        glmCopyVector3f(tempNormal, normal);
      }
    }
  }
  return distance;
}

double glmRayMeshIntersectiond(const double rO[3], const double rV[3],
                            const double *vertex, DWORD stride,
                            const DWORD *indicates, DWORD iCount, double intersection[3], double normal[3])
{
  double distance=-1, tempDistance;
  BYTE *v;
  double tempNormal[3], tempIntersection[3];
  double *vert1, *vert2, *vert3;
  int i;

  if (stride==0) stride=3*sizeof(double);

  v=(BYTE *)vertex;
  for (i=0; i<iCount; i+=3)
  {
    vert1=(double *)(v+indicates[i]*stride);
    vert2=(double *)(v+indicates[i+1]*stride);
	vert3=(double *)(v+indicates[i+2]*stride);
	glmComputeTriangleNormald(tempNormal, vert1, vert2, vert3);
	tempDistance=glmRayTriangleIntersectiond(rO, rV, vert1, vert2, vert3, tempNormal, tempIntersection);
    if (tempDistance>=0&&(tempDistance<distance||distance<0))
    {
      distance=tempDistance;
      if (intersection)
      {
        glmCopyVector3d(tempIntersection, intersection);
      }
      if (normal)
      {
        glmCopyVector3d(tempNormal, normal);
      }
    }
  }
  return distance;
}


BOOL glmCircleLineIntersectionf(const float cO[2], const float cR,
                                const float lO[2], const float lN[2],
                                float intersection[2], float *r)
{
  // Locals.
  float d;

  // Calculate the distance to the line.
  d=-glmDotProduct2f(lN, cO)+glmDotProduct2f(lN, lO);
  if (d<-cR||d>cR) return FALSE;

  // Store the intersection point.
  if (intersection)
  {
    intersection[0]=cO[0]+d*lN[0];
    intersection[1]=cO[1]+d*lN[1];
  }

  if (r)
  {
    *r=cR*sin(acos(d/cR));
  }

  return TRUE;
}

BOOL glmCircleLineIntersectiond(const double cO[2], const double cR,
                                const double lO[2], const double lN[2],
                                double intersection[2], double *r)
{
  // Locals.
  double d;

  // Calculate the distance to the line.
  d=-glmDotProduct2d(lN, cO)+glmDotProduct2d(lN, lO);
  if (d<-cR||d>cR) return FALSE;

  // Store the intersection point.
  if (intersection)
  {
    intersection[0]=cO[0]+d*lN[0];
    intersection[1]=cO[1]+d*lN[1];
  }

  if (r)
  {
    *r=cR*sin(acos(d/cR));
  }

  return TRUE;
}

BOOL glmSpherePlaneIntersectionf(const float sO[3], const float sR,
                                 const float pO[3], const float pN[3],
                                 float intersection[3], float *r)
{
  // Locals.
  float d;

  // Calculate the distance to the line.
  d=-glmDotProduct3f(pN, sO)+glmDotProduct3f(pN, pO);
  if (d<-sR||d>sR) return FALSE;

  // Store the intersection point.
  if (intersection)
  {
    intersection[0]=sO[0]+d*pN[0];
    intersection[1]=sO[1]+d*pN[1];
    intersection[2]=sO[2]+d*pN[2];
  }

  if (r)
  {
    *r=sR*sin(acos(d/sR));
  }

  return TRUE;
}


BOOL glmSpherePlaneIntersectiond(const double sO[3], const double sR,
                                 const double pO[3], const double pN[3],
                                 double intersection[3], double *r)
{
  // Locals.
  double d;

  // Calculate the distance to the line.
  d=-glmDotProduct3d(pN, sO)+glmDotProduct3d(pN, pO);
  if (d<-sR||d>sR) return FALSE;

  // Store the intersection point.
  if (intersection)
  {
    intersection[0]=sO[0]+d*pN[0];
    intersection[1]=sO[1]+d*pN[1];
    intersection[2]=sO[2]+d*pN[2];
  }

  if (r)
  {
    *r=sR*sin(acos(d/sR));
  }

  return TRUE;
}

BOOL glmHyperSphereSpaceIntersectionf(const float sO[4], const float sR,
                                 const float pO[4], const float pN[4],
                                 float intersection[4], float *r)
{
  // Locals.
  float d;

  // Calculate the distance to the line.
  d=-glmDotProduct4f(pN, sO)+glmDotProduct4f(pN, pO);
  if (d<-sR||d>sR) return FALSE;

  // Store the intersection point.
  if (intersection)
  {
    intersection[0]=sO[0]+d*pN[0];
    intersection[1]=sO[1]+d*pN[1];
    intersection[2]=sO[2]+d*pN[2];
    intersection[3]=sO[3]+d*pN[3];
  }

  if (r)
  {
    *r=sR*sin(acos(d/sR));
  }

  return TRUE;
}

BOOL glmHyperSphereSpaceIntersectiond(const double sO[4], const double sR,
                                 const double pO[4], const double pN[4],
                                 double intersection[4], double *r)
{
  // Locals.
  float d;

  // Calculate the distance to the line.
  d=-glmDotProduct4d(pN, sO)+glmDotProduct4d(pN, pO);
  if (d<-sR||d>sR) return FALSE;

  // Store the intersection point.
  if (intersection)
  {
    intersection[0]=sO[0]+d*pN[0];
    intersection[1]=sO[1]+d*pN[1];
    intersection[2]=sO[2]+d*pN[2];
    intersection[3]=sO[3]+d*pN[3];
  }

  if (r)
  {
    *r=sR*sin(acos(d/sR));
  }

  return TRUE;
}

void glmExtractFrustumf(float frustum[6][4], const float view[16], const float proj[16])
{
  float clip[16];
  float t;

  glmMultMatrix4f(clip, proj, view);

  // The RIGHT plane
  frustum[0][0] = clip[ 3] - clip[ 0];
  frustum[0][1] = clip[ 7] - clip[ 4];
  frustum[0][2] = clip[11] - clip[ 8];
  frustum[0][3] = clip[15] - clip[12];
  t=glmNormalizeVector3f(frustum[0]);
  frustum[0][3]/=t;

  // The LEFT plane
  frustum[1][0] = clip[ 3] + clip[ 0];
  frustum[1][1] = clip[ 7] + clip[ 4];
  frustum[1][2] = clip[11] + clip[ 8];
  frustum[1][3] = clip[15] + clip[12];
  t=glmNormalizeVector3f(frustum[1]);
  frustum[1][3]/= t;

  // The BOTTOM plane
  frustum[2][0] = clip[ 3] + clip[ 1];
  frustum[2][1] = clip[ 7] + clip[ 5];
  frustum[2][2] = clip[11] + clip[ 9];
  frustum[2][3] = clip[15] + clip[13];
  t=glmNormalizeVector3f(frustum[2]);
  frustum[2][3] /= t;

  // The TOP plane
  frustum[3][0] = clip[ 3] - clip[ 1];
  frustum[3][1] = clip[ 7] - clip[ 5];
  frustum[3][2] = clip[11] - clip[ 9];
  frustum[3][3] = clip[15] - clip[13];
  t=glmNormalizeVector3f(frustum[3]);
  frustum[3][3] /= t;

  // The FAR plane
  frustum[4][0] = clip[ 3] - clip[ 2];
  frustum[4][1] = clip[ 7] - clip[ 6];
  frustum[4][2] = clip[11] - clip[10];
  frustum[4][3] = clip[15] - clip[14];
  t=glmNormalizeVector3f(frustum[4]);
  frustum[4][3] /= t;

  // The NEAR plane
  frustum[5][0] = clip[ 3] + clip[ 2];
  frustum[5][1] = clip[ 7] + clip[ 6];
  frustum[5][2] = clip[11] + clip[10];
  frustum[5][3] = clip[15] + clip[14];
  t=glmNormalizeVector3f(frustum[5]);
  frustum[5][3] /= t;
}

void glmExtractFrustumd(double frustum[6][4], const double view[16], const double proj[16])
{
  double clip[16];
  double t;

  glmMultMatrix4d(clip,  proj, view);

  // The RIGHT plane
  frustum[0][0] = clip[ 3] - clip[ 0];
  frustum[0][1] = clip[ 7] - clip[ 4];
  frustum[0][2] = clip[11] - clip[ 8];
  frustum[0][3] = clip[15] - clip[12];
  t=glmNormalizeVector3d(frustum[0]);
  frustum[0][3]/=t;

  // The LEFT plane
  frustum[1][0] = clip[ 3] + clip[ 0];
  frustum[1][1] = clip[ 7] + clip[ 4];
  frustum[1][2] = clip[11] + clip[ 8];
  frustum[1][3] = clip[15] + clip[12];
  t=glmNormalizeVector3d(frustum[1]);
  frustum[1][3]/= t;

  // The BOTTOM plane
  frustum[2][0] = clip[ 3] + clip[ 1];
  frustum[2][1] = clip[ 7] + clip[ 5];
  frustum[2][2] = clip[11] + clip[ 9];
  frustum[2][3] = clip[15] + clip[13];
  t=glmNormalizeVector3d(frustum[2]);
  frustum[2][3] /= t;

  // The TOP plane
  frustum[3][0] = clip[ 3] - clip[ 1];
  frustum[3][1] = clip[ 7] - clip[ 5];
  frustum[3][2] = clip[11] - clip[ 9];
  frustum[3][3] = clip[15] - clip[13];
  t=glmNormalizeVector3d(frustum[3]);
  frustum[3][3] /= t;

  // The FAR plane
  frustum[4][0] = clip[ 3] - clip[ 2];
  frustum[4][1] = clip[ 7] - clip[ 6];
  frustum[4][2] = clip[11] - clip[10];
  frustum[4][3] = clip[15] - clip[14];
  t=glmNormalizeVector3d(frustum[4]);
  frustum[4][3] /= t;

  // The NEAR plane
  frustum[5][0] = clip[ 3] + clip[ 2];
  frustum[5][1] = clip[ 7] + clip[ 6];
  frustum[5][2] = clip[11] + clip[10];
  frustum[5][3] = clip[15] + clip[14];
  t=glmNormalizeVector3d(frustum[5]);
  frustum[5][3] /= t;
}

BOOL glmIsPointInFrustumf(const float frustum[6][4], const float p[3])
{
  int i;
  for (i=0; i<6; i++)
  {
    if(glmDotProduct3f(frustum[i], p)+frustum[i][3]<=0) return FALSE;
  }
  return TRUE;
}

BOOL glmIsPointInFrustumd(const double frustum[6][4], const double p[3])
{
  int i;
  for (i=0; i<6; i++)
  {
    if(glmDotProduct3d(frustum[i], p)+frustum[i][3]<=0) return FALSE;
  }
  return TRUE;
}

BOOL glmIsSphereInFrustumf(const float frustum[6][4], const float sphereOrigin[3], float radius)
{
  int i;
  for(i=0; i<6; i++)
  {
    if(glmDotProduct3f(frustum[i], sphereOrigin)+frustum[i][3]<=-radius) return FALSE;
  }
  return TRUE;
}

BOOL glmIsSphereInFrustumd(const double frustum[6][4], const double sphereOrigin[3], double radius)
{
  int i;
  for(i=0; i<6; i++)
  {
    if(glmDotProduct3d(frustum[i], sphereOrigin)+frustum[i][3]<=-radius) return FALSE;
  }
  return TRUE;
}

BOOL glmIsCubeInFrustumf(const float frustum[6][4], float x, float y, float z, float size)
{
  int p, c;

  for(p=0; p<6; p++)
  {
    c = 0;
    if(frustum[p][0]*(x-size)+frustum[p][1]*(y-size)+frustum[p][2]*(z-size)+frustum[p][3]>0) c++;
    if(frustum[p][0]*(x+size)+frustum[p][1]*(y-size)+frustum[p][2]*(z-size)+frustum[p][3]>0) c++;
    if(frustum[p][0]*(x-size)+frustum[p][1]*(y+size)+frustum[p][2]*(z-size)+frustum[p][3]>0) c++;
    if(frustum[p][0]*(x+size)+frustum[p][1]*(y+size)+frustum[p][2]*(z-size)+frustum[p][3]>0) c++;
    if(frustum[p][0]*(x-size)+frustum[p][1]*(y-size)+frustum[p][2]*(z+size)+frustum[p][3]>0) c++;
    if(frustum[p][0]*(x+size)+frustum[p][1]*(y-size)+frustum[p][2]*(z+size)+frustum[p][3]>0) c++;
    if(frustum[p][0]*(x-size)+frustum[p][1]*(y+size)+frustum[p][2]*(z+size)+frustum[p][3]>0) c++;
    if(frustum[p][0]*(x+size)+frustum[p][1]*(y+size)+frustum[p][2]*(z+size)+frustum[p][3]>0) c++;
    if(c==0) return FALSE;
  }

  return TRUE;
}

BOOL glmIsCubeInFrustumd(const double frustum[6][4], double x, double y, double z, double size)
{
  int p, c;

  for(p=0; p<6; p++)
  {
    c = 0;
    if(frustum[p][0]*(x-size)+frustum[p][1]*(y-size)+frustum[p][2]*(z-size)+frustum[p][3]>0) c++;
    if(frustum[p][0]*(x+size)+frustum[p][1]*(y-size)+frustum[p][2]*(z-size)+frustum[p][3]>0) c++;
    if(frustum[p][0]*(x-size)+frustum[p][1]*(y+size)+frustum[p][2]*(z-size)+frustum[p][3]>0) c++;
    if(frustum[p][0]*(x+size)+frustum[p][1]*(y+size)+frustum[p][2]*(z-size)+frustum[p][3]>0) c++;
    if(frustum[p][0]*(x-size)+frustum[p][1]*(y-size)+frustum[p][2]*(z+size)+frustum[p][3]>0) c++;
    if(frustum[p][0]*(x+size)+frustum[p][1]*(y-size)+frustum[p][2]*(z+size)+frustum[p][3]>0) c++;
    if(frustum[p][0]*(x-size)+frustum[p][1]*(y+size)+frustum[p][2]*(z+size)+frustum[p][3]>0) c++;
    if(frustum[p][0]*(x+size)+frustum[p][1]*(y+size)+frustum[p][2]*(z+size)+frustum[p][3]>0) c++;
    if(c==0) return FALSE;
  }

  return TRUE;
}

/*BOOL glmIsTriangleInFrustum(const float frustum[6][4], const float v1[3], const float v2[3], const float v3[3])
{
   int f, p;

   for( f = 0; f < 6; f++ )
   {
      for( p = 0; p < numpoints; p++ )
      {
         if( frustum[f][0] * pointlist[p].x + frustum[f][1] * pointlist[p].y + frustum[f][2] * pointlist[p].z + frustum[f][3] > 0 )
            break;
      }
      if( p == numpoints )
         return false;
   }
   return true;
}
*/



float glmComputeDistanceToPlane3f(const float planeOrigin[3], const float planeNormal[3], const float p[3])
{
  return glmDotProduct3f(planeNormal, p)-glmDotProduct3f(planeNormal, planeOrigin);
}

double glmComputeDistanceToPlane3d(const double planeOrigin[3], const double planeNormal[3], const double p[3])
{
  return glmDotProduct3d(planeNormal, p)-glmDotProduct3d(planeNormal, planeOrigin);
}

void glmClosestPointOnPlane3f(float result[3], const float planeOrigin[3], const float planeNormal[3], const float p[3])
{
  float temp[3];
  float d;

  glmSubtractVector3f(temp, p, planeOrigin);
  d=glmDotProduct3f(temp, planeNormal);

  result[0]=p[0]-d*planeNormal[0];
  result[1]=p[1]-d*planeNormal[1];
  result[2]=p[2]-d*planeNormal[2];
}

void glmClosestPointOnPlane3d(double result[3], const double planeOrigin[3], const double planeNormal[3], const double p[3])
{
  double temp[3];
  double d;

  glmSubtractVector3d(temp, p, planeOrigin);
  d=glmDotProduct3d(temp, planeNormal);

  result[0]=p[0]-d*planeNormal[0];
  result[1]=p[1]-d*planeNormal[1];
  result[2]=p[2]-d*planeNormal[2];
}

float glmComputeDistanceToLine3f(const float lineOrigin[3], const float lineDirection[3], const float p[3])
{
  float c[3], temp[3];
  float t;

  glmMakeVector3f(c, lineOrigin, p);
  t=glmDotProduct3f(lineDirection, c);

  temp[0]=lineOrigin[0]+t*lineDirection[0];
  temp[1]=lineOrigin[1]+t*lineDirection[1];
  temp[2]=lineOrigin[2]+t*lineDirection[2];

  return glmComputeDistance3f(p, temp);
}

double glmComputeDistanceToLine3d(const double lineOrigin[3], const double lineDirection[3], const double p[3])
{
  double c[3], temp[3];
  double t;

  glmMakeVector3d(c, lineOrigin, p);
  t=glmDotProduct3d(lineDirection, c);

  temp[0]=lineOrigin[0]+t*lineDirection[0];
  temp[1]=lineOrigin[1]+t*lineDirection[1];
  temp[2]=lineOrigin[2]+t*lineDirection[2];

  return glmComputeDistance3d(p, temp);
}

void glmClosestPointOnTriangle3f(float result[3],
                                 const float v1[3], const float v2[3], const float v3[3],
                                 const float p[3])
{
  float temp[3];
  float normal[3];
  float d1, d2, d3;
  float t1[3], t2[3], t3[3];

  glmComputeTriangleNormalf(normal, v1, v2, v3);
  glmClosestPointOnPlane3f(temp, v1, normal, p);

  if (glmIsPointInTriangle3f(p, v1, v2, v3))
  {
    glmCopyVector3f(p, result);
    return;
  }

  glmClosestPointOnLineSegment3f(t1, v1, v2, temp);
  glmClosestPointOnLineSegment3f(t2, v2, v3, temp);
  glmClosestPointOnLineSegment3f(t3, v3, v1, temp);

  d1=glmComputeDistance3f(t1, p);
  d2=glmComputeDistance3f(t2, p);
  d3=glmComputeDistance3f(t3, p);

  switch (glmMinIndex3f(d1, d2, d3, 0, 1, 2))
  {
    case 0:
      glmCopyVector3f(t1, result);
      return;
    break;
    case 1:
      glmCopyVector3f(t2, result);
      return;
    break;
    case 2:
      glmCopyVector3f(t3, result);
      return;
    break;
  }
}

float glmComputeTriangleArea2f(const float v1[2], const float v2[2], const float v3[2])
{
  float a, b, c, p;

  a=glmComputeDistance2f(v2, v3);
  b=glmComputeDistance2f(v1, v3);
  c=glmComputeDistance2f(v1, v2);

  p=(a+b+c)/2.0;

  return sqrt(p*(p-a)*(p-b)*(p-c));
}

double glmComputeTriangleArea2d(const double v1[2], const double v2[2], const double v3[2])
{
  float a, b, c, p;

  a=glmComputeDistance2d(v2, v3);
  b=glmComputeDistance2d(v1, v3);
  c=glmComputeDistance2d(v1, v2);

  p=(a+b+c)/2.0;

  return sqrt(p*(p-a)*(p-b)*(p-c));
}

float glmComputeTriangleArea3f(const float v1[3], const float v2[3], const float v3[3])
{
  float a, b, c, p;

  a=glmComputeDistance3f(v2, v3);
  b=glmComputeDistance3f(v1, v3);
  c=glmComputeDistance3f(v1, v2);

  p=(a+b+c)/2.0;

  return sqrt(p*(p-a)*(p-b)*(p-c));
}

double glmComputeTriangleArea3d(const double v1[3], const double v2[3], const double v3[3])
{
  float a, b, c, p;

  a=glmComputeDistance3d(v2, v3);
  b=glmComputeDistance3d(v1, v3);
  c=glmComputeDistance3d(v1, v2);

  p=(a+b+c)/2.0;

  return sqrt(p*(p-a)*(p-b)*(p-c));
}

void glmComputeBarycentricCoordinates2f(float result[3],
                                        const float v1[2], const float v2[2], const float v3[2],
                                        const float p[2])
{
  float a, b, c, u, v, D;
  float tv1[2], tv2[2], tv3[2];

  glmMakeVector2f(tv1, v3, v1);
  glmMakeVector2f(tv2, v3, v2);
  a=glmDotProduct2f(tv1, tv1);
  b=glmDotProduct2f(tv2, tv1);
  c=glmDotProduct2f(tv2, tv2);

  glmMakeVector2f(tv3, v3, p);
  u=glmDotProduct2f(tv3, tv1);
  v=glmDotProduct2f(tv3, tv2);

  D = a*c-b*b;

  if (D==0)
  {
    result[0]=0;
	result[1]=0;
    result[2]=0;
    return;
  }

  D = 1.0/D;

  result[0]=(c*u-b*v)*D;
  result[1]=(-b*u+a*v)*D;
  result[2]=1-result[0]-result[1];
}

void glmComputeBarycentricCoordinates2d(double result[3],
                                        const double v1[2], const double v2[2], const double v3[2],
                                        const double p[2])
{
  double a, b, c, u, v, D;
  double tv1[2], tv2[2], tv3[2];

  glmMakeVector2d(tv1, v3, v1);
  glmMakeVector2d(tv2, v3, v2);
  a=glmDotProduct2d(tv1, tv1);
  b=glmDotProduct2d(tv2, tv1);
  c=glmDotProduct2d(tv2, tv2);

  glmMakeVector2d(tv3, v3, p);
  u=glmDotProduct2d(tv3, tv1);
  v=glmDotProduct2d(tv3, tv2);

  D = a*c-b*b;

  if (D==0)
  {
    result[0]=0;
	result[1]=0;
    result[2]=0;
    return;
  }

  D = 1.0/D;

  result[0]=(c*u-b*v)*D;
  result[1]=(-b*u+a*v)*D;
  result[2]=1-result[0]-result[1];
}



void glmComputeBarycentricCoordinates3f(float result[3],
                                        const float v1[3], const float v2[3], const float v3[3],
                                        const float p[3])
{
  float a, b, c, u, v, D;
  float tv1[3], tv2[3], tv3[3];

  glmMakeVector3f(tv1, v3, v1);
  glmMakeVector3f(tv2, v3, v2);
  a=glmDotProduct3f(tv1, tv1);
  b=glmDotProduct3f(tv2, tv1);
  c=glmDotProduct3f(tv2, tv2);

  glmMakeVector3f(tv3, v3, p);
  u=glmDotProduct3f(tv3, tv1);
  v=glmDotProduct3f(tv3, tv2);

  D = a*c-b*b;

  if (D==0)
  {
    result[0]=0;
	result[1]=0;
    result[2]=0;
    return;
  }

  D = 1.0/D;

  result[0]=(c*u-b*v)*D;
  result[1]=(-b*u+a*v)*D;
  result[2]=1-result[0]-result[1];
}

void glmComputeBarycentricCoordinates3d(double result[3],
                                        const double v1[3], const double v2[3], const double v3[3],
                                        const double p[3])
{
  double a, b, c, u, v, D;
  double tv1[3], tv2[3], tv3[3];

  glmMakeVector3d(tv1, v3, v1);
  glmMakeVector3d(tv2, v3, v2);
  a=glmDotProduct3d(tv1, tv1);
  b=glmDotProduct3d(tv2, tv1);
  c=glmDotProduct3d(tv2, tv2);

  glmMakeVector3d(tv3, v3, p);
  u=glmDotProduct3d(tv3, tv1);
  v=glmDotProduct3d(tv3, tv2);

  D = a*c-b*b;

  if (D==0)
  {
    result[0]=0;
	result[1]=0;
    result[2]=0;
    return;
  }

  D = 1.0/D;

  result[0]=(c*u-b*v)*D;
  result[1]=(-b*u+a*v)*D;
  result[2]=1-result[0]-result[1];
}

/* This function tests if a moving sphere hits a triangle.
   The sphere's center point is "p" and radius "radius". Sphere's movement vector is "movement"
   and it will tested against triangle that has vertices "v1", "v2", "v3", and normal "normal".*/
BOOL glmMovingSphereTriangleIntersectionf(
  const float p[3], float radius, const float movement[3],
  const float v1[3], const float v2[3], const float v3[3], const float normal[3],
  float *distance, float result[3], float slideNormal[3], float slide[3])
{
  // Locals.
  float direction[3];       // Unit lenght direction vector of the movement.
  float invertDirection[3]; // Points to opposite direction than "direction".
  float inverseNormal[3];   // Points to opposite direction than the triangles normal.
  float sphereIntersectionPoint[3]; // The point on the sphere that will collide with the triangle.
  float planeIntersectionPoint[3];  // The point on the triangle's plane closest to the sphere.
  float triangleIntersectionPoint[3]; // the point on the triangle that collides with the triangle.
  float distanceToTravel=0;   // Lenght of the movement vector.
  float planeDistance=0;  // Sphere's center distance to the plane.

  // Calculate the movement's unit vector.
  distanceToTravel=glmComputeVectorLength3f(movement);
  if (distanceToTravel<=0) return FALSE;    // No movement-> no collision.
  direction[0]=movement[0]/distanceToTravel;
  direction[1]=movement[1]/distanceToTravel;
  direction[2]=movement[2]/distanceToTravel;

  // Ignore backfaces.
  if (glmDotProduct3f(normal, direction)>=0) return FALSE;
  glmGenerateReverseVector3f(inverseNormal, normal);

  // Determine the distance from the plane to the center of the sphere.
  planeDistance=glmRayPlaneIntersectionf(p, inverseNormal, v1, normal, planeIntersectionPoint);
  if (planeDistance<0) return FALSE;  // No intersection?

  // If the plane doesn't intersects the sphere.
  if (planeDistance>radius)
  {
    // Calculate the sphere intersection point.
    sphereIntersectionPoint[0]=p[0]-radius*normal[0];
    sphereIntersectionPoint[1]=p[1]-radius*normal[1];
    sphereIntersectionPoint[2]=p[2]-radius*normal[2];

    // Calculate the plane intersection point
    *distance=glmRayPlaneIntersectionf(sphereIntersectionPoint, direction, v1, normal,
                                       planeIntersectionPoint);
    if (*distance>distanceToTravel) return FALSE;
  }

  // Find a point on the triangle closest to the plane intersection point.
  if (glmIsPointInTriangle3f(planeIntersectionPoint, v1, v2, v3))
    glmCopyVector3f(planeIntersectionPoint, triangleIntersectionPoint);
  else
    glmClosestPointOnTriangle3f(triangleIntersectionPoint, v1, v2, v3, planeIntersectionPoint);


  if (glmComputeDistance3f(triangleIntersectionPoint, p)<=radius)
  {
    *distance=0;
    glmCopyVector3f(p, result);
  }
  else
  {
    // Invert the movement direction vector
    glmGenerateReverseVector3f(invertDirection, direction);

    // Using the polygonIntersectionPoint, we need to reverse-intersect with the sphere
    *distance=glmRaySphereIntersectionf(triangleIntersectionPoint, invertDirection, p, radius, NULL);

    // Was there an intersection with the sphere?
    if (*distance<0||*distance>distanceToTravel) return FALSE;

    // Where is the center of the sphere at time of the intersection?
    result[0]=p[0]+direction[0]**distance;
    result[1]=p[1]+direction[1]**distance;
    result[2]=p[2]+direction[2]**distance;
  }

  // Calculate the sliding normal
  slideNormal[0]=result[0]-triangleIntersectionPoint[0];
  slideNormal[1]=result[1]-triangleIntersectionPoint[1];
  slideNormal[2]=result[2]-triangleIntersectionPoint[2];
  if (glmNormalizeVector3f(slideNormal)==0)
  {
    glmCopyVector3f(normal, slideNormal);
  }

  // Last escape test.
  if (glmDotProduct3f(slideNormal, direction)>=0) return FALSE;

  // Calculate the sliding vector.
  direction[0]*=distanceToTravel-*distance;
  direction[1]*=distanceToTravel-*distance;
  direction[2]*=distanceToTravel-*distance;
  glmProjectVector3f(slide, direction, slideNormal);

  // Intersection found!
  return TRUE;
}

int glmMovingSphereMeshIntersectionf(
       const float p[3], float radius, const float movement[3],
       const float *vertex,    DWORD stride,
       const DWORD *indicates,  DWORD iCount,
       float *distance, float result[3], float intersectionNormal[3], float slide[3])
{
  float normal[3];
  float tempDistance;
  float shortestDistance=-1;
  float tempResult[3];
  float tempiNormal[3];
  float tempSlide[3];
  int num=-1;
  DWORD i;
  BYTE *v;
  float *vert1, *vert2, *vert3;

  if (stride==0) stride=3*sizeof(float);

  v=(BYTE *)vertex;
  for (i=0; i<iCount; i+=3)
  {
    vert1=(float *)(v+indicates[i]*stride);
    vert2=(float *)(v+indicates[i+1]*stride);
    vert3=(float *)(v+indicates[i+2]*stride);

    glmComputeTriangleNormalf(normal, vert1, vert2, vert3);

    if (glmMovingSphereTriangleIntersectionf(p, radius, movement,
          vert1, vert2, vert3, normal,
          &tempDistance, tempResult, tempiNormal, tempSlide))
    {
      if (shortestDistance<0||tempDistance<shortestDistance)
      {
        num=i/3;
		shortestDistance=tempDistance;
        if (distance) *distance=shortestDistance;
        if (result) glmCopyVector3f(tempResult, result);
        if (intersectionNormal) glmCopyVector3f(tempiNormal, intersectionNormal);
        if (slide) glmCopyVector3f(tempSlide, slide);
      }
    }
  }

  return num;
}

BOOL glmIsPointInFrontOfPlanef(const float p[3], const float planeOrigin[3], const float planeNormal[3])
{
  if (glmDotProduct3f(p, planeNormal)-glmDotProduct3f(planeOrigin, planeNormal)<0) return FALSE;
  return TRUE;
}

BOOL glmIsPointInFrontOfPlaned(const double p[3], const double planeOrigin[3], const double planeNormal[3])
{
  if (glmDotProduct3d(p, planeNormal)-glmDotProduct3d(planeOrigin, planeNormal)<0) return FALSE;
  return TRUE;
}

void glmCutTriangleByPlanef(const float p1[3], const float p2[3], const float p3[3],
                            const float *t1, const float *t2, const float *t3, DWORD tComponents,
                            const float n1[3], const float n2[3], const float n3[3],

                            const float planeOrigin[3], const float planeNormal[3],

                            float *frontVertex, DWORD fvStride,
							float *frontTexCoord, DWORD ftStride,
							float *frontNormal,  DWORD fnStride,
							DWORD *frontIndicates, DWORD *frontVCount, DWORD *frontICount,

                            float *backVertex, DWORD bvStride,
							float *backTexCoord, DWORD btStride,
							float *backNormal,  DWORD bnStride,
							DWORD *backIndicates, DWORD *backVCount, DWORD *backICount)
{
  const float *p[3];
  const float *t[3];
  const float *n[3];
  DWORD i;
  BOOL sideFlag[3]={FALSE, FALSE, FALSE};
  float ip[3];
  float d;
  DWORD fvc=0, bvc=0;

  p[0]=p1; p[1]=p2; p[2]=p3;
  t[0]=t1; t[1]=t2; t[2]=t3;
  n[0]=n1; n[1]=n2; n[2]=n3;

  if (fvStride<3*sizeof(float)) fvStride=3*sizeof(float);
  if (ftStride<tComponents*sizeof(float)) ftStride=tComponents*sizeof(float);
  if (fnStride<3*sizeof(float)) fnStride=3*sizeof(float);
  if (bvStride<3*sizeof(float)) bvStride=3*sizeof(float);
  if (btStride<tComponents*sizeof(float)) btStride=tComponents*sizeof(float);
  if (bnStride<3*sizeof(float)) bnStride=3*sizeof(float);

  // Determine vertex sides
  for (i=0; i<3; i++)
  {
    if (glmIsPointInFrontOfPlanef(p[i], planeOrigin, planeNormal)) sideFlag[i]=TRUE;
  }

  // Cut the triangle
  for (i=0; i<3; i++)
  {
    if (sideFlag[i])
    {
      if (frontVertex) glmCopyVector3f(p[i], (BYTE *)frontVertex+fvStride*fvc);
      if (frontTexCoord && t[i]) glmCopyVectorf(t[i], (BYTE *)frontTexCoord+ftStride*fvc, tComponents);
      if (frontNormal && n[i]) glmCopyVector3f(n[i], (BYTE *)frontNormal+fnStride*fvc);
      fvc++;
    }
    else
    {
      if (backVertex) glmCopyVector3f(p[i], (BYTE *)backVertex+bvStride*bvc);
      if (backTexCoord && t[i]) glmCopyVectorf(t[i], (BYTE *)backTexCoord+btStride*bvc, tComponents);
      if (backNormal && n[i]) glmCopyVector3f(n[i], (BYTE *)backNormal+bnStride*bvc);
      bvc++;
    }

    if (sideFlag[i]!=sideFlag[(i+1)%3])
    {
      glmLineSegmentPlaneIntersectionf(p[i], p[(i+1)%3], planeOrigin, planeNormal, ip, &d);
      if (frontVertex) glmCopyVector3f(ip, (BYTE *)frontVertex+fvStride*fvc);
      if (frontTexCoord && t[i] && t[(i+1)%3])
        glmInterpolateVectorf((BYTE *)frontTexCoord+ftStride*fvc, t[i], t[(i+1)%3], d, tComponents);
      if (frontNormal && n[i] && n[(i+1)%3])
      {
        glmInterpolateVector3f((BYTE *)frontNormal+fnStride*fvc, n[i], n[(i+1)%3], d);
        glmNormalizeVector3f((BYTE *)frontNormal+fnStride*fvc);
      }
      if (backVertex) glmCopyVector3f(ip, (BYTE *)backVertex+bvStride*bvc);
      if (backTexCoord && t[i] && t[(i+1)%3])
        glmInterpolateVectorf((BYTE *)backTexCoord+btStride*bvc, t[i], t[(i+1)%3], d, tComponents);
      if (backNormal && n[i] && n[(i+1)%3])
      {
        glmInterpolateVector3f((BYTE *)backNormal+bnStride*bvc, n[i], n[(i+1)%3], d);
        glmNormalizeVector3f((BYTE *)backNormal+bnStride*bvc);
      }
      fvc++;
      bvc++;
    }
  }

  // Set indicates
  if (frontICount) *frontICount=0;
  if (fvc>=3)
  {
    if (frontICount) *frontICount=3;
    if (frontIndicates)
    {
      frontIndicates[0]=0;
      frontIndicates[1]=1;
      frontIndicates[2]=2;
    }
    if (fvc>=4)
    {
      if (frontICount) *frontICount=6;
      if (frontIndicates)
      {
        frontIndicates[3]=0;
        frontIndicates[4]=2;
        frontIndicates[5]=3;
      }
    }
  }
  if (backICount) *backICount=0;
  if (bvc>=3)
  {
    if (backICount) *backICount=3;
    if (backIndicates)
    {
      backIndicates[0]=0;
      backIndicates[1]=1;
      backIndicates[2]=2;
    }
    if (bvc>=4)
    {
      if (backICount) *backICount=6;
      if (backIndicates)
      {
        backIndicates[3]=0;
        backIndicates[4]=2;
        backIndicates[5]=3;
      }
    }
  }

  if (frontVCount) *frontVCount=fvc;
  if (backVCount) *backVCount=bvc;
}


void glmCutTriangleByPlaned(const double p1[3], const double p2[3], const double p3[3],
                            const double *t1, const double *t2, const double *t3, DWORD tComponents,
                            const double n1[3], const double n2[3], const double n3[3],

                            const double planeOrigin[3], const double planeNormal[3],

                            double *frontVertex, DWORD fvStride,
							double *frontTexCoord, DWORD ftStride,
							double *frontNormal,  DWORD fnStride,
							DWORD *frontIndicates, DWORD *frontVCount, DWORD *frontICount,

                            double *backVertex, DWORD bvStride,
							double *backTexCoord, DWORD btStride,
							double *backNormal,  DWORD bnStride,
							DWORD *backIndicates, DWORD *backVCount, DWORD *backICount)
{
  const double *p[3];
  const double *t[3];
  const double *n[3];
  DWORD i;
  BOOL sideFlag[3]={FALSE, FALSE, FALSE};
  double ip[3];
  double d;
  DWORD fvc=0, bvc=0;

  p[0]=p1; p[1]=p2; p[2]=p3;
  t[0]=t1; t[1]=t2; t[2]=t3;
  n[0]=n1; n[1]=n2; n[2]=n3;

  if (fvStride<3*sizeof(double)) fvStride=3*sizeof(double);
  if (ftStride<tComponents*sizeof(double)) ftStride=tComponents*sizeof(double);
  if (fnStride<3*sizeof(double)) fnStride=3*sizeof(double);
  if (bvStride<3*sizeof(double)) bvStride=3*sizeof(double);
  if (btStride<tComponents*sizeof(double)) btStride=tComponents*sizeof(double);
  if (bnStride<3*sizeof(double)) bnStride=3*sizeof(double);

  // Determine vertex sides
  for (i=0; i<3; i++)
  {
    if (glmIsPointInFrontOfPlaned(p[i], planeOrigin, planeNormal)) sideFlag[i]=TRUE;
  }

  // Cut the triangle
  for (i=0; i<3; i++)
  {
    if (sideFlag[i])
    {
      if (frontVertex) glmCopyVector3d(p[i], (BYTE *)frontVertex+fvStride*fvc);
      if (frontTexCoord && t[i]) glmCopyVectord(t[i], (BYTE *)frontTexCoord+ftStride*fvc, tComponents);
      if (frontNormal && n[i]) glmCopyVector3d(n[i], (BYTE *)frontNormal+fnStride*fvc);
      fvc++;
    }
    else
    {
      if (backVertex) glmCopyVector3d(p[i], (BYTE *)backVertex+bvStride*bvc);
      if (backTexCoord && t[i]) glmCopyVectord(t[i], (BYTE *)backTexCoord+btStride*bvc, tComponents);
      if (backNormal && n[i]) glmCopyVector3d(n[i], (BYTE *)backNormal+bnStride*bvc);
      bvc++;
    }

    if (sideFlag[i]!=sideFlag[(i+1)%3])
    {
      glmLineSegmentPlaneIntersectiond(p[i], p[(i+1)%3], planeOrigin, planeNormal, ip, &d);
      if (frontVertex) glmCopyVector3d(ip, (BYTE *)frontVertex+fvStride*fvc);
      if (frontTexCoord && t[i] && t[(i+1)%3])
        glmInterpolateVectord((BYTE *)frontTexCoord+ftStride*fvc, t[i], t[(i+1)%3], d, tComponents);
      if (frontNormal && n[i] && n[(i+1)%3])
      {
        glmInterpolateVector3d((BYTE *)frontNormal+fnStride*fvc, n[i], n[(i+1)%3], d);
        glmNormalizeVector3d((BYTE *)frontNormal+fnStride*fvc);
      }
      if (backVertex) glmCopyVector3d(ip, (BYTE *)backVertex+bvStride*bvc);
      if (backTexCoord && t[i] && t[(i+1)%3])
        glmInterpolateVectord((BYTE *)backTexCoord+btStride*bvc, t[i], t[(i+1)%3], d, tComponents);
      if (backNormal && n[i] && n[(i+1)%3])
      {
        glmInterpolateVector3d((BYTE *)backNormal+bnStride*bvc, n[i], n[(i+1)%3], d);
        glmNormalizeVector3d((BYTE *)backNormal+bnStride*bvc);
      }
      fvc++;
      bvc++;
    }
  }

  // Set indicates
  if (frontICount) *frontICount=0;
  if (fvc>=3)
  {
    if (frontICount) *frontICount=3;
    if (frontIndicates)
    {
      frontIndicates[0]=0;
      frontIndicates[1]=1;
      frontIndicates[2]=2;
    }
    if (fvc>=4)
    {
      if (frontICount) *frontICount=6;
      if (frontIndicates)
      {
        frontIndicates[3]=0;
        frontIndicates[4]=2;
        frontIndicates[5]=3;
      }
    }
  }
  if (backICount) *backICount=0;
  if (bvc>=3)
  {
    if (backICount) *backICount=3;
    if (backIndicates)
    {
      backIndicates[0]=0;
      backIndicates[1]=1;
      backIndicates[2]=2;
    }
    if (bvc>=4)
    {
      if (backICount) *backICount=6;
      if (backIndicates)
      {
        backIndicates[3]=0;
        backIndicates[4]=2;
        backIndicates[5]=3;
      }
    }
  }

  if (frontVCount) *frontVCount=fvc;
  if (backVCount) *backVCount=bvc;
}

void glmCutPolygonByPlanef(const float *vertex, DWORD vStride,
                           const float *texCoord, DWORD tStride, DWORD tComponents,
                           const float *normal, DWORD nStride,
						   DWORD vCount,

                           const float planeOrigin[3], const float planeNormal[3],

                           float *frontVertex, DWORD fvStride,
                           float *frontTexCoord, DWORD ftStride,
                           float *frontNormal,  DWORD fnStride,
                           DWORD *frontVCount,

                           float *backVertex, DWORD bvStride,
                           float *backTexCoord, DWORD btStride,
                           float *backNormal,  DWORD bnStride,
                           DWORD *backVCount)
{
  DWORD i;
  BOOL *sideFlag;
  float ip[3];
  float d;
  DWORD fvc=0, bvc=0;

  sideFlag = (BOOL*) malloc( vCount * sizeof(BOOL) );

  memset( sideFlag, 0, sizeof( BOOL ) * vCount );

  if (vStride<3*sizeof(float)) vStride=3*sizeof(float);
  if (tStride<tComponents*sizeof(float)) tStride=tComponents*sizeof(float);
  if (nStride<3*sizeof(float)) nStride=3*sizeof(float);
  if (fvStride<3*sizeof(float)) fvStride=3*sizeof(float);
  if (ftStride<tComponents*sizeof(float)) ftStride=tComponents*sizeof(float);
  if (fnStride<3*sizeof(float)) fnStride=3*sizeof(float);
  if (bvStride<3*sizeof(float)) bvStride=3*sizeof(float);
  if (btStride<tComponents*sizeof(float)) btStride=tComponents*sizeof(float);
  if (bnStride<3*sizeof(float)) bnStride=3*sizeof(float);

  // Determine vertex sides
  for (i=0; i<vCount; i++)
  {
    if (glmIsPointInFrontOfPlanef((BYTE *)vertex+i*vStride,
                                  planeOrigin, planeNormal)) sideFlag[i]=TRUE;
  }

  // Cut the polygon
  for (i=0; i<vCount; i++)
  {
    if (sideFlag[i])
    {
      if (frontVertex)
        glmCopyVector3f((BYTE *)vertex+i*vStride,
                        (BYTE *)frontVertex+fvStride*fvc);
      if (frontTexCoord && texCoord)
        glmCopyVectorf((BYTE *)texCoord+i*tStride,
                       (BYTE *)frontTexCoord+ftStride*fvc, tComponents);
      if (frontNormal && normal)
        glmCopyVector3f((BYTE *)normal+i*nStride,
                        (BYTE *)frontNormal+fnStride*fvc);
      fvc++;
    }
    else
    {
      if (backVertex)
        glmCopyVector3f((BYTE *)vertex+i*vStride,
                        (BYTE *)backVertex+bvStride*bvc);
      if (backTexCoord && texCoord)
        glmCopyVectorf((BYTE *)texCoord+i*tStride,
                       (BYTE *)backTexCoord+btStride*bvc, tComponents);
      if (backNormal && normal)
        glmCopyVector3f((BYTE *)normal+i*nStride,
                        (BYTE *)backNormal+bnStride*bvc);
      bvc++;
    }

    if (sideFlag[i]!=sideFlag[(i+1)%vCount])
    {
      glmLineSegmentPlaneIntersectionf((BYTE *)vertex+i*vStride,
                                       (BYTE *)vertex+((i+1)%vCount)*vStride,
                                       planeOrigin, planeNormal, ip, &d);
      if (frontVertex) glmCopyVector3f(ip, (BYTE *)frontVertex+fvStride*fvc);
      if (frontTexCoord && texCoord)
        glmInterpolateVectorf((BYTE *)frontTexCoord+ftStride*fvc,
                              (BYTE *)texCoord+i*tStride,
                              (BYTE *)texCoord+((i+1)%vCount)*tStride,
                              d, tComponents);
      if (frontNormal && normal)
      {
        glmInterpolateVector3f((BYTE *)frontNormal+fnStride*fvc,
                               (BYTE *)normal+i*nStride,
                               (BYTE *)normal+((i+1)%vCount)*nStride, d);
        glmNormalizeVector3f((BYTE *)frontNormal+fnStride*fvc);
      }
      if (backVertex) glmCopyVector3f(ip, (BYTE *)backVertex+bvStride*bvc);
      if (backTexCoord && texCoord)
        glmInterpolateVectorf((BYTE *)backTexCoord+btStride*bvc,
                              (BYTE *)texCoord+i*tStride,
                              (BYTE *)texCoord+((i+1)%vCount)*tStride,
                              d, tComponents);
      if (backNormal && normal)
      {
        glmInterpolateVector3f((BYTE *)backNormal+bnStride*bvc,
                               (BYTE *)normal+i*nStride,
                               (BYTE *)normal+((i+1)%vCount)*nStride, d);
        glmNormalizeVector3f((BYTE *)backNormal+bnStride*bvc);
      }
      fvc++;
      bvc++;
    }
  }

  if (frontVCount) *frontVCount=fvc;
  if (backVCount) *backVCount=bvc;

  free( sideFlag );

}
