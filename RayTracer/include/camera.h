#ifndef __CAMERA_H__
#define __CAMERA_H__

#include "utils.h"

enum{
	PERSPECTIVE,
	ORTHOGONAL
};

typedef struct Camera_s{
	
		// Perspective parameters
	elem_t		perspectiveAngle;
	elem_t		nearPlaneDistance;
	elem_t		aspectRatio;

		// Orthogonal parameters
	elem_t		left;
	elem_t		right;
	elem_t		bottom;
	elem_t		top;
	
		// Screen dimension
	elem_t		screenWidth;
	elem_t		screenHeight;
	
		// All kinds of matrices
	matrix44_t  matCamSpaceToWorldSpace;
	matrix44_t  matScreenSpaceToCamSpace;
	matrix44_t	matScreenSpaceToWorldSpace;
	matrix44_t  matCamSpaceToScreenSpace;
	matrix44_t  matWorldSpaceToCamSpace;

	int			mode;
} Camera_t;

typedef struct Ray_s {
	vector4_t start;
	vector4_t dir;
} Ray_t;


#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "camera.h"
#include "glmvector.h"
#include "glmmatrix.h"


#ifdef __cplusplus
	extern "C" {
#endif

void camTransformRay( const Ray_t *inRay, const matrix44_t mat, Ray_t *out );
void camGetWorldSpaceRayFromScreenSpace( double x, double y, Ray_t *ray );
void camGetCamSpaceRayFromScreenSpace( double x, double y, Ray_t *ray );
void camInitPerspectiveView();
void camInitOrthogonalView();
void camSetImageParameters( const elem_t screenWidth, const elem_t screenHeight );
void camSetOrthogonalParameters( const elem_t Left, const elem_t Right, const elem_t Bottom, const elem_t Top );
void camSetPerspectiveParameters( const elem_t angle, const elem_t aspectRatio, const elem_t nearPlane );
void camSetOrientation( const vector3_t Pos, const vector3_t LookAt, const vector3_t Up );
void camGetPosition( vector3_t Pos );

void camGetScreenSpacePointFromWorldSpace(vector2_t screenSpace, const vector3_t worldSpace);
void camGetScreenSpacePointFromCamSpace( vector2_t screenSpace, const vector3_t camSpace);

#ifdef __cplusplus
	}
#endif


#endif

