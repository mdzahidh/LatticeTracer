#ifndef __CVOLUME_H__
#define __CVOLUME_H__

#include "utils.h"
#include "glmvector.h"
#include "glmmatrix.h"
#include "camera.h"


class CVolume {
public:
	
	// These are same as grid spacing, i.e the distance between two grid points
	// in each axis direction.
	double m_xScale;
	double m_yScale;
	double m_zScale;

	unsigned int m_nRaysCulled;
	
	vector3_t m_corners[8];
	double     m_XSize;
	double     m_YSize;
	double     m_ZSize;
	
	double	   m_HXSize; //half of m_XSize;
	double	   m_HYSize; //half of m_XSize;
	double	   m_HZSize; //half of m_XSize;

	double	  m_mappingScaleTrans[3];
	double	  m_invScale[3];	
	
		//There will be 6 faces and each one will have three indices of [ Vmin, Vmax, Silent Ordinate ]
	int		  m_faces[6][3];	
		
		//Normals of each faces;
	vector3_t	  m_faceNormals[6];

		//Each faces will have a Normal Dot V = a
	double		  m_faceDistances[6];

	double		  m_radiusSquared;
	
	inline void CalculateFaceNormal( int v0, int v1, int v2, vector3_t N )
	{
		vector3_t v0v1, v0v2;
		glmSubtractVector3d( v0v1, m_corners[v1], m_corners[v0] );
		glmSubtractVector3d( v0v2, m_corners[v2], m_corners[v0] );
		glmCrossProduct3d( N, v0v1, v0v2 );
	}

	CVolume( double lx, double ly, double lz, double xScale, double yScale, double zScale )
	{
		m_xScale = xScale;
		m_yScale = yScale;
		m_zScale = zScale;
		m_nRaysCulled = 0;
		
		double halfX = lx / 2.0f;
		double halfY = ly / 2.0f;
		double halfZ = lz / 2.0f;
		
		m_HXSize = halfX;
		m_HYSize = halfY;
		m_HZSize = halfZ;

		m_XSize = halfX * 2;
		m_YSize = halfY * 2;
		m_ZSize = halfZ * 2;				
		
		glmInitVector3d( m_corners[0], -halfX, -halfY, halfZ );
		glmInitVector3d( m_corners[1],  halfX, -halfY, halfZ );
		glmInitVector3d( m_corners[2],  halfX,  halfY, halfZ );
		glmInitVector3d( m_corners[3], -halfX,  halfY, halfZ );				

		glmInitVector3d( m_corners[4], -halfX, -halfY, -halfZ );
		glmInitVector3d( m_corners[5],  halfX, -halfY, -halfZ );
		glmInitVector3d( m_corners[6],  halfX,  halfY, -halfZ );
		glmInitVector3d( m_corners[7], -halfX,  halfY, -halfZ );
		
		printf("Calculated Bounding Box:\n\tMin: [%f %f %f]\n\tMax: [%f %f %f]\n", -halfX, -halfY, -halfZ, halfX, halfY, halfZ ); 

		vector3_t diagonal;		
		glmSubtractVector3d( diagonal, m_corners[6], m_corners[0] );
		m_radiusSquared = glmDotProduct3d( diagonal, diagonal ) / 4;
/*
		m_mappingScaleTrans[X] = (dataXSize-1) * 0.5f;
		m_mappingScaleTrans[Y] = (dataYSize-1) * 0.5f;
		m_mappingScaleTrans[Z] = (dataZSize-1) * 0.5f;
*/
		m_mappingScaleTrans[X] = (lx/xScale) * 0.5f;
		m_mappingScaleTrans[Y] = (ly/yScale) * 0.5f;
		m_mappingScaleTrans[Z] = (lz/zScale) * 0.5f;

		m_invScale[X] = 1.0f / xScale;
		m_invScale[Y] = 1.0f / yScale;
		m_invScale[Z] = 1.0f / zScale;

		m_faces[0][0] = 0; m_faces[0][1] = 2; m_faces[0][2] = Z;
		CalculateFaceNormal( 0,1,3, m_faceNormals[0] );
		m_faceDistances[0] = glmDotProduct3d( m_faceNormals[0], m_corners[0] );


		m_faces[1][0] = 5; m_faces[1][1] = 2; m_faces[1][2] = X;
		CalculateFaceNormal( 1,5,2, m_faceNormals[1] );
		m_faceDistances[1] = glmDotProduct3d( m_faceNormals[1], m_corners[1] );

		m_faces[2][0] = 4; m_faces[2][1] = 6; m_faces[2][2] = Z;
		CalculateFaceNormal( 5,4,6, m_faceNormals[2] );
		m_faceDistances[2] = glmDotProduct3d( m_faceNormals[2], m_corners[5] );

		
		m_faces[3][0] = 4; m_faces[3][1] = 3; m_faces[3][2] = X;
		CalculateFaceNormal( 4,0,7, m_faceNormals[3] );
		m_faceDistances[3] = glmDotProduct3d( m_faceNormals[3], m_corners[4] );

		m_faces[4][0] = 4; m_faces[4][1] = 1; m_faces[4][2] = Y;
		CalculateFaceNormal( 4,5,0, m_faceNormals[4] );
		m_faceDistances[4] = glmDotProduct3d( m_faceNormals[4], m_corners[4] );

		m_faces[5][0] = 7; m_faces[5][1] = 2; m_faces[5][2] = Y;
		CalculateFaceNormal( 7,3,6, m_faceNormals[5] );
		m_faceDistances[5] = glmDotProduct3d( m_faceNormals[5], m_corners[7] );
	}

	CVolume( int dataXSize, int dataYSize, int dataZSize, double xScale, double yScale, double zScale )
	{
		m_xScale = xScale;
		m_yScale = yScale;
		m_zScale = zScale;
		m_nRaysCulled = 0;
		
		double halfX = xScale * (dataXSize-1) / 2.0f;
		double halfY = yScale * (dataYSize-1) / 2.0f;
		double halfZ = zScale * (dataZSize-1) / 2.0f;
		
		m_HXSize = halfX;
		m_HYSize = halfY;
		m_HZSize = halfZ;

		m_XSize = halfX * 2;
		m_YSize = halfY * 2;
		m_ZSize = halfZ * 2;				
		
		glmInitVector3d( m_corners[0], -halfX, -halfY, halfZ );
		glmInitVector3d( m_corners[1],  halfX, -halfY, halfZ );
		glmInitVector3d( m_corners[2],  halfX,  halfY, halfZ );
		glmInitVector3d( m_corners[3], -halfX,  halfY, halfZ );				

		glmInitVector3d( m_corners[4], -halfX, -halfY, -halfZ );
		glmInitVector3d( m_corners[5],  halfX, -halfY, -halfZ );
		glmInitVector3d( m_corners[6],  halfX,  halfY, -halfZ );
		glmInitVector3d( m_corners[7], -halfX,  halfY, -halfZ );
		
		printf("Calculated Bounding Box:\n\tMin: [%f %f %f]\n\tMax: [%f %f %f]\n", -halfX, -halfY, -halfZ, halfX, halfY, halfZ );

		vector3_t diagonal;		
		glmSubtractVector3d( diagonal, m_corners[6], m_corners[0] );
		m_radiusSquared = glmDotProduct3d( diagonal, diagonal ) / 4;

		m_mappingScaleTrans[X] = (dataXSize-1) * 0.5f;
		m_mappingScaleTrans[Y] = (dataYSize-1) * 0.5f;
		m_mappingScaleTrans[Z] = (dataZSize-1) * 0.5f;

		m_invScale[X] = 1.0f / xScale;
		m_invScale[Y] = 1.0f / yScale;
		m_invScale[Z] = 1.0f / zScale;

		m_faces[0][0] = 0; m_faces[0][1] = 2; m_faces[0][2] = Z;
		CalculateFaceNormal( 0,1,3, m_faceNormals[0] );
		m_faceDistances[0] = glmDotProduct3d( m_faceNormals[0], m_corners[0] );


		m_faces[1][0] = 5; m_faces[1][1] = 2; m_faces[1][2] = X;
		CalculateFaceNormal( 1,5,2, m_faceNormals[1] );
		m_faceDistances[1] = glmDotProduct3d( m_faceNormals[1], m_corners[1] );

		m_faces[2][0] = 4; m_faces[2][1] = 6; m_faces[2][2] = Z;
		CalculateFaceNormal( 5,4,6, m_faceNormals[2] );
		m_faceDistances[2] = glmDotProduct3d( m_faceNormals[2], m_corners[5] );

		
		m_faces[3][0] = 4; m_faces[3][1] = 3; m_faces[3][2] = X;
		CalculateFaceNormal( 4,0,7, m_faceNormals[3] );
		m_faceDistances[3] = glmDotProduct3d( m_faceNormals[3], m_corners[4] );

		m_faces[4][0] = 4; m_faces[4][1] = 1; m_faces[4][2] = Y;
		CalculateFaceNormal( 4,5,0, m_faceNormals[4] );
		m_faceDistances[4] = glmDotProduct3d( m_faceNormals[4], m_corners[4] );

		m_faces[5][0] = 7; m_faces[5][1] = 2; m_faces[5][2] = Y;
		CalculateFaceNormal( 7,3,6, m_faceNormals[5] );
		m_faceDistances[5] = glmDotProduct3d( m_faceNormals[5], m_corners[7] );
	}	
	
	double GetXSize() { return m_XSize; };
	double GetYSize() { return m_YSize; };
	double GetZSize() { return m_ZSize; };

	double GetHXSize() { return m_HXSize; };
	double GetHYSize() { return m_HYSize; };
	double GetHZSize() { return m_HZSize; };

	void FindExtentPoints( const double *mat, vector3_t min, vector3_t max )
	{
		vector4_t transV;		
		
		glmMultVectorByMatrix4d( transV, m_corners[0], mat );
		glmMultVectorByScalar3d( transV, 1.0f / transV[W] );
		
		glmInitVector3d( min, transV[X], transV[Y], transV[Z] );
		glmInitVector3d( max , transV[X], transV[Y], transV[Z] );
						
		for( int i=1;i<8;++i){
			
			glmMultVectorByMatrix4d( transV, m_corners[i], mat );
			glmMultVectorByScalar3d( transV, 1.0f / transV[W] );
			
			if( min[X] > transV[X] ) min[X] = transV[X];
			if( min[Y] > transV[Y] ) min[Y] = transV[Y];
			if( min[Z] > transV[Z] ) min[Z] = transV[Z];
		
			if( max[X] < transV[X] ) max[X] = transV[X];
			if( max[Y] < transV[Y] ) max[Y] = transV[Y];
			if( max[Z] < transV[Z] ) max[Z] = transV[Z];
		}		
	}
	
	inline void FindUVW( const vector4_t objectSpaceP, vector3_t uvw )
	{
					
		uvw[X] = objectSpaceP[X] / m_XSize + 0.5f;
		uvw[Y] = objectSpaceP[Y] / m_YSize + 0.5f;
		uvw[Z] = objectSpaceP[Z] / m_ZSize + 0.5f;
		
	}

	void FindUVW( const vector4_t objectSpaceP[], vector3_t uvw[], int num )
	{
		for(int i=0;i<num;++i){
			
			uvw[i][X] = objectSpaceP[i][X] / m_XSize + 0.5f;
			uvw[i][Y] = objectSpaceP[i][Y] / m_YSize + 0.5f;
			uvw[i][Z] = objectSpaceP[i][Z] / m_ZSize + 0.5f;
		}
	}

	inline void FindXYZScaleOffset( vector3_t offset )
	{
		offset[X] = m_mappingScaleTrans[X];
		offset[Y] = m_mappingScaleTrans[Y];
		offset[Z] = m_mappingScaleTrans[Z];
	}
	
	inline void FindXYZScaleGradient( vector3_t grad )
	{
		grad[X] *= m_invScale[X];
		grad[Y] *= m_invScale[Y];
		grad[Z] *= m_invScale[Z];
	}
	
	
	inline void FindObjectSpacePoint( const vector3_t gridXYZ, vector3_t objectSpaceP)
	{
		objectSpaceP[X] = (gridXYZ[X] - m_mappingScaleTrans[X]) * m_xScale;
		objectSpaceP[Y] = (gridXYZ[Y] - m_mappingScaleTrans[Y]) * m_yScale;
		objectSpaceP[Z] = (gridXYZ[Z] - m_mappingScaleTrans[Z]) * m_zScale;
	}

	inline void FindXYZ( const vector3_t objectSpaceP, vector3_t xyz )
	{
		
		xyz[X] = objectSpaceP[X] * m_invScale[X] + m_mappingScaleTrans[X];
		xyz[Y] = objectSpaceP[Y] * m_invScale[Y] + m_mappingScaleTrans[Y];
		xyz[Z] = objectSpaceP[Z] * m_invScale[Z] + m_mappingScaleTrans[Z];			
	}

	bool FindEntryExitPoints( const Ray_t *pRay, vector3_t entry, vector3_t exit );

	inline int DoesPassSphereTest( const Ray_t *pRay )
	{
		double abDot = glmDotProduct3d( pRay->start, pRay->dir );
		double aaDot = glmDotProduct3d( pRay->start, pRay->start );
		double bbDot = glmDotProduct3d( pRay->dir, pRay->dir );

		double e = abDot * abDot - (( aaDot - m_radiusSquared) * bbDot);

		if( e < EPSILON && e > -EPSILON ) return 1;
		else if( e > EPSILON ) return 2;
		
		return 0;
	}	

	bool CheckBound( int face, vector3_t r );
	void GetScales( vector3_t scales)
	{
		scales[0] = m_xScale;
		scales[1] = m_yScale;
		scales[2] = m_zScale;
	}
};
		
#endif

