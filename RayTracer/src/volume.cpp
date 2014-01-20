#include "volume.h"

bool CVolume::CheckBound( int face, vector3_t r )
{
	vector3_t &vMin = m_corners[ m_faces[face][0] ];
	vector3_t &vMax = m_corners[ m_faces[face][1] ];

	switch( m_faces[face][2] ){
		case X:
			if( ((r[Y] + EPSILON) < vMin[Y] ) || ( (r[Y] - EPSILON) > vMax[Y]) || ( (r[Z] + EPSILON) < vMin[Z]) || ( (r[Z] - EPSILON) > vMax[Z]) ) return false;
			break;
		
		case Y:
			if( ((r[X] + EPSILON) < vMin[X] )|| ( (r[X] - EPSILON) > vMax[X]) || ( (r[Z] + EPSILON) < vMin[Z]) || ( (r[Z] - EPSILON) > vMax[Z]) ) return false;
			break;
		
		case Z:
			if( ((r[Y] + EPSILON)< vMin[Y] )|| ( (r[Y] - EPSILON) > vMax[Y]) || ( (r[X] + EPSILON) < vMin[X]) || ( (r[X] - EPSILON) > vMax[X]) ) return false;			
			break;
	};

	return true;
}

bool CVolume::FindEntryExitPoints( const Ray_t *pRay, vector3_t entry, vector3_t exit )
{
/*
	int nPoints = DoesPassSphereTest( pRay );
	if( !nPoints ) {
		m_nRaysCulled++;
		return false;
	}
*/	
	
	double maxLambda = 0;
	double minLambda = 0;

	bool firstFound = false;

	for( int i =0; i < 6; ++i ){
		
		double bDotN = glmDotProduct3d( pRay->dir, m_faceNormals[i] );
		double  aDotN = glmDotProduct3d( pRay->start, m_faceNormals[i] );	
	
		if( (bDotN > EPSILON) || (bDotN < - EPSILON) ){
			double lambda = (m_faceDistances[i] - aDotN) / bDotN;			
			
			vector3_t r;
			vector3_t scaledDir;

			glmCopyVector3d(pRay->dir, scaledDir);
			glmMultVectorByScalar3d( scaledDir, lambda );
			glmAddVector3d( r, pRay->start, scaledDir );

			if( CheckBound( i, r ) ){

				if( !firstFound ){
					firstFound = true;
					minLambda = maxLambda = lambda;
				}
				else{
					if( lambda < minLambda ) minLambda = lambda;
					if( lambda > maxLambda ) maxLambda = lambda;
				}
			}
		}		
		
	}

	if( !firstFound ) return false;

	if( maxLambda < 0) return false;

	if( minLambda <0 ) minLambda = 0;

	vector3_t scaledDir;

	glmCopyVector3d(pRay->dir, scaledDir);
	glmMultVectorByScalar3d( scaledDir, minLambda );
	glmAddVector3d( entry, pRay->start, scaledDir );

	glmCopyVector3d( pRay->dir, scaledDir);
	glmMultVectorByScalar3d( scaledDir, maxLambda );
	glmAddVector3d( exit, pRay->start, scaledDir );

	return true;
}
