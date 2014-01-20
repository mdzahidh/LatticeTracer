#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "camera.h"
#include "glmvector.h"
#include "glmmatrix.h"

Camera_t g_Camera = {0};

void camGetPosition( vector3_t Pos )
{	
	Pos[0] = g_Camera.matCamSpaceToWorldSpace[12];
	Pos[1] = g_Camera.matCamSpaceToWorldSpace[13];
	Pos[2] = g_Camera.matCamSpaceToWorldSpace[14];
}

void camTransformRay( const Ray_t *inRay, const matrix44_t mat, Ray_t *out )
{
	glmMultVectorByMatrix4d( out->start, mat, inRay->start );
	glmMultVectorByMatrix4d( out->dir, mat, inRay->dir );
}

void camGetCamSpaceRayFromScreenSpace( double x, double y, Ray_t *ray )
{
	vector4_t in = {x,y,0,1};
	glmMultVectorByMatrix4d( ray->start, in, g_Camera.matScreenSpaceToCamSpace );
	
	if( g_Camera.mode == PERSPECTIVE ){
		glmCopyVector3d( ray->start, ray->dir );
		glmNormalizeVector3d( ray->dir );
		ray->dir[W] = 0;
	}
	else{
		// Just copy the Look vector of the camera
		glmCopyVector3d( &g_Camera.matCamSpaceToWorldSpace[12], ray->dir ); 
		ray->dir[W] = 0;
	}

}

void camGetWorldSpaceRayFromScreenSpace( double x, double y, Ray_t *ray )
{	
	vector4_t in = {x,y,0,1};
	glmMultVectorByMatrix4d( ray->start, in, g_Camera.matScreenSpaceToWorldSpace );
	
	if( g_Camera.mode == PERSPECTIVE ){
		glmCopyVector3d( ray->start, ray->dir );
		ray->dir[X] -= g_Camera.matCamSpaceToWorldSpace[12];
		ray->dir[Y] -= g_Camera.matCamSpaceToWorldSpace[13];
		ray->dir[Z] -= g_Camera.matCamSpaceToWorldSpace[14];
		glmNormalizeVector3d( ray->dir );
		ray->dir[W] = 0;
	}
	else{
		// Just copy the Look vector of the camera
		glmCopyVector3d( &g_Camera.matCamSpaceToWorldSpace[12], ray->dir ); 
		ray->dir[W] = 0;
	}
}

void camGetScreenSpacePointFromCamSpace( vector2_t screenSpace, const vector3_t camSpace) 
{
	vector4_t CS;
	vector4_t ret;

	glmCopyVector3d( camSpace, CS );
	CS[W] = 1;

	if( g_Camera.mode == PERSPECTIVE){
		
		CS[X] = CS[X] * g_Camera.nearPlaneDistance / CS[Z];
		CS[Y] = CS[Y] * g_Camera.nearPlaneDistance / CS[Z];
		CS[Z] = g_Camera.nearPlaneDistance;
		

		glmMultVectorByMatrix4d( ret, CS, g_Camera.matCamSpaceToScreenSpace );
		
	}
	else{
		
		glmMultVectorByMatrix4d( ret, CS, g_Camera.matCamSpaceToScreenSpace );

	}

	screenSpace[X] = ret[X];
	screenSpace[Y] = ret[Y];
}

void camGetScreenSpacePointFromWorldSpace( vector3_t screenSpace, const vector3_t worldSpace)
{
	vector4_t WS;
	vector4_t camSpace;

	glmCopyVector3d( worldSpace, WS);
	WS[W] = 1;

	glmMultVectorByMatrix4d( camSpace,WS,g_Camera.matWorldSpaceToCamSpace );
	camGetScreenSpacePointFromCamSpace( screenSpace, camSpace );
}


void camInitPerspectiveView()
{
	elem_t tanHalfAngle = tan( g_Camera.perspectiveAngle / 2.0f );

	g_Camera.matScreenSpaceToCamSpace[ 0] = 2 * g_Camera.nearPlaneDistance * tanHalfAngle / g_Camera.screenWidth;
	g_Camera.matScreenSpaceToCamSpace[ 1] = 0;
	g_Camera.matScreenSpaceToCamSpace[ 2] = 0;
	g_Camera.matScreenSpaceToCamSpace[ 3] = 0;
	
	g_Camera.matScreenSpaceToCamSpace[ 4] = 0;
	g_Camera.matScreenSpaceToCamSpace[ 5] = -2 * g_Camera.nearPlaneDistance * tanHalfAngle  / (  g_Camera.screenHeight * g_Camera.aspectRatio );
	g_Camera.matScreenSpaceToCamSpace[ 6] = 0;
	g_Camera.matScreenSpaceToCamSpace[ 7] = 0;
	
	g_Camera.matScreenSpaceToCamSpace[ 8] = 0;
	g_Camera.matScreenSpaceToCamSpace[ 9] = 0;
	g_Camera.matScreenSpaceToCamSpace[10] = 1;
	g_Camera.matScreenSpaceToCamSpace[11] = 0;
	
	g_Camera.matScreenSpaceToCamSpace[12] = -g_Camera.nearPlaneDistance * tanHalfAngle;
	g_Camera.matScreenSpaceToCamSpace[13] =  g_Camera.nearPlaneDistance * tanHalfAngle / g_Camera.aspectRatio;
	g_Camera.matScreenSpaceToCamSpace[14] =  g_Camera.nearPlaneDistance;
	g_Camera.matScreenSpaceToCamSpace[15] =  1;

	
	g_Camera.matCamSpaceToScreenSpace[0] = g_Camera.screenWidth / (2 * g_Camera.nearPlaneDistance * tanHalfAngle);
	g_Camera.matCamSpaceToScreenSpace[1] = 0;
	g_Camera.matCamSpaceToScreenSpace[2] = 0;
	g_Camera.matCamSpaceToScreenSpace[3] = 0;

	g_Camera.matCamSpaceToScreenSpace[4] = 0;
	g_Camera.matCamSpaceToScreenSpace[5] = - (g_Camera.screenHeight * g_Camera.aspectRatio) / (2 * g_Camera.nearPlaneDistance * tanHalfAngle);
	g_Camera.matCamSpaceToScreenSpace[6] = 0;
	g_Camera.matCamSpaceToScreenSpace[7] = 0;

	g_Camera.matCamSpaceToScreenSpace[8]  = 0;
	g_Camera.matCamSpaceToScreenSpace[9]  = 0;
	g_Camera.matCamSpaceToScreenSpace[10] = 1;
	g_Camera.matCamSpaceToScreenSpace[11] = 0;

	g_Camera.matCamSpaceToScreenSpace[12] = g_Camera.screenWidth/2.0;
	g_Camera.matCamSpaceToScreenSpace[13] = g_Camera.screenHeight/2.0;
	g_Camera.matCamSpaceToScreenSpace[14] = -g_Camera.nearPlaneDistance;
	g_Camera.matCamSpaceToScreenSpace[15] = 1;

	glmMultMatrix4d( g_Camera.matScreenSpaceToWorldSpace, g_Camera.matCamSpaceToWorldSpace, g_Camera.matScreenSpaceToCamSpace );

	g_Camera.mode = PERSPECTIVE;

}

void camInitOrthogonalView()
{	
	elem_t viewWidth = g_Camera.right - g_Camera.left;
	elem_t viewHeight = g_Camera.top - g_Camera.bottom;

	g_Camera.matScreenSpaceToCamSpace[ 0] = viewWidth / g_Camera.screenWidth;
	g_Camera.matScreenSpaceToCamSpace[ 1] = 0;
	g_Camera.matScreenSpaceToCamSpace[ 2] = 0;
	g_Camera.matScreenSpaceToCamSpace[ 3] = 0;
	
	g_Camera.matScreenSpaceToCamSpace[ 4] = 0;
	g_Camera.matScreenSpaceToCamSpace[ 5] = -viewHeight / g_Camera.screenHeight;
	g_Camera.matScreenSpaceToCamSpace[ 6] = 0;
	g_Camera.matScreenSpaceToCamSpace[ 7] = 0;
	
	g_Camera.matScreenSpaceToCamSpace[ 8] = 0;
	g_Camera.matScreenSpaceToCamSpace[ 9] = 0;
	g_Camera.matScreenSpaceToCamSpace[10] = 1;
	g_Camera.matScreenSpaceToCamSpace[11] = 0;
	
	g_Camera.matScreenSpaceToCamSpace[12] = g_Camera.left;
	g_Camera.matScreenSpaceToCamSpace[13] = g_Camera.top;
	g_Camera.matScreenSpaceToCamSpace[14] = 0;
	g_Camera.matScreenSpaceToCamSpace[15] = 1;

	g_Camera.matCamSpaceToScreenSpace[0] = g_Camera.screenWidth / viewWidth;
	g_Camera.matCamSpaceToScreenSpace[1] = 0;
	g_Camera.matCamSpaceToScreenSpace[2] = 0;
	g_Camera.matCamSpaceToScreenSpace[3] = 0;

	g_Camera.matCamSpaceToScreenSpace[4] = 0;
	g_Camera.matCamSpaceToScreenSpace[5] = - (g_Camera.screenHeight) / viewHeight;
	g_Camera.matCamSpaceToScreenSpace[6] = 0;
	g_Camera.matCamSpaceToScreenSpace[7] = 0;

	g_Camera.matCamSpaceToScreenSpace[8]  = 0;
	g_Camera.matCamSpaceToScreenSpace[9]  = 0;
	g_Camera.matCamSpaceToScreenSpace[10] = 1;
	g_Camera.matCamSpaceToScreenSpace[11] = 0;

	g_Camera.matCamSpaceToScreenSpace[12] = - g_Camera.left * g_Camera.screenWidth / viewWidth;
	g_Camera.matCamSpaceToScreenSpace[13] = g_Camera.top * g_Camera.screenHeight / viewHeight;
	g_Camera.matCamSpaceToScreenSpace[14] = 0;
	g_Camera.matCamSpaceToScreenSpace[15] = 1;

	glmMultMatrix4d( g_Camera.matScreenSpaceToWorldSpace, g_Camera.matCamSpaceToWorldSpace, g_Camera.matScreenSpaceToCamSpace );

	g_Camera.mode = ORTHOGONAL;
}

void camSetImageParameters( const elem_t screenWidth, const elem_t screenHeight )
{
	g_Camera.screenWidth  = screenWidth;
	g_Camera.screenHeight = screenHeight;
}

void camSetOrthogonalParameters( const elem_t Left, const elem_t Right, const elem_t Bottom, const elem_t Top )
{
	g_Camera.left = Left;
	g_Camera.right = Right;
	g_Camera.bottom = Bottom;
	g_Camera.top = Top;
}

void camSetPerspectiveParameters( const elem_t angle, const elem_t aspectRatio, const elem_t nearPlane )
{
	g_Camera.aspectRatio = aspectRatio;
	g_Camera.perspectiveAngle = (elem_t) DEG2RAD(angle);
	g_Camera.nearPlaneDistance = nearPlane;
}

void camSetOrientation( const vector3_t Pos, const vector3_t LookAt, const vector3_t Up )
{
	vector3_t Right;
	vector3_t Look;
	vector3_t CorrectUp;

	glmSubtractVector3d( Look, LookAt, Pos );
	glmNormalizeVector3d( Look );

	
	glmCrossProduct3d( Right, Look, Up );
	glmNormalizeVector3d( Right );
	glmCrossProduct3d( CorrectUp, Right, Look );


	g_Camera.matCamSpaceToWorldSpace[ 0] = Right[X];	
	g_Camera.matCamSpaceToWorldSpace[ 1] = Right[Y];
	g_Camera.matCamSpaceToWorldSpace[ 2] = Right[Z];
	g_Camera.matCamSpaceToWorldSpace[ 3] = 0;
	
	g_Camera.matCamSpaceToWorldSpace[ 4] = CorrectUp[X];
	g_Camera.matCamSpaceToWorldSpace[ 5] = CorrectUp[Y];
	g_Camera.matCamSpaceToWorldSpace[ 6] = CorrectUp[Z];
	g_Camera.matCamSpaceToWorldSpace[ 7] = 0;
	
	g_Camera.matCamSpaceToWorldSpace[ 8] = Look[X];
	g_Camera.matCamSpaceToWorldSpace[ 9] = Look[Y];
	g_Camera.matCamSpaceToWorldSpace[10] = Look[Z];
	g_Camera.matCamSpaceToWorldSpace[11] = 0;
	
	g_Camera.matCamSpaceToWorldSpace[12] = Pos[X];
	g_Camera.matCamSpaceToWorldSpace[13] = Pos[Y];
	g_Camera.matCamSpaceToWorldSpace[14] = Pos[Z];
	g_Camera.matCamSpaceToWorldSpace[15] = 1;

	
	g_Camera.matWorldSpaceToCamSpace[ 0] = Right[X];
	g_Camera.matWorldSpaceToCamSpace[ 1] = CorrectUp[X];
	g_Camera.matWorldSpaceToCamSpace[ 2] = Look[X];
	g_Camera.matWorldSpaceToCamSpace[ 3] = 0;

	g_Camera.matWorldSpaceToCamSpace[ 4] = Right[Y];
	g_Camera.matWorldSpaceToCamSpace[ 5] = CorrectUp[Y];
	g_Camera.matWorldSpaceToCamSpace[ 6] = Look[Y];
	g_Camera.matWorldSpaceToCamSpace[ 7] = 0;

	g_Camera.matWorldSpaceToCamSpace[ 8] = Right[Z];
	g_Camera.matWorldSpaceToCamSpace[ 9] = CorrectUp[Z];
	g_Camera.matWorldSpaceToCamSpace[10] = Look[Z];
	g_Camera.matWorldSpaceToCamSpace[11] = 0;

	g_Camera.matWorldSpaceToCamSpace[12] = -glmDotProduct3d( Right, Pos );
	g_Camera.matWorldSpaceToCamSpace[13] = -glmDotProduct3d( CorrectUp,Pos);
	g_Camera.matWorldSpaceToCamSpace[14] = -glmDotProduct3d( Look,Pos);
	g_Camera.matWorldSpaceToCamSpace[15] = 1;

}


