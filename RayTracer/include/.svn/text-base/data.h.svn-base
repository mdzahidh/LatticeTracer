#ifndef __CDATA_H__
#define __CDATA_H__

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string>

#include "utils.h"
#include "glmvector.h"
#include "volume.h"

#define one_over_6 0.16666666666666666
#define two_over_24 0.083333333333333329
#define six_over_120 0.05

#define MAP_ENTRY_BY_FLAG( map, flag ) map[#flag] = flag

#define GUARD_BAND 32


elem_t *** LoadVUDFileByte( const char *fname, int &xsize, int &ysize, int &zsize );
elem_t *** LoadVUDFileFloat( const char *fname, int &xsize, int &ysize, int &zsize );

typedef unsigned char atomic_t;

class CBaseGradientComponentEstimator{
public:
	virtual ~CBaseGradientComponentEstimator()
	{
	}
	virtual elem_t GetGradientComponent( int x, int y, int z ) = 0;
};

class CBaseData{
public:
	int m_xSize;  // along x-axis , number of samples
	int m_ySize;  // along y-axis, """""
	int m_zSize;  // along z-axis, """"

	elem_t    ***m_pppData;
	vector3_t ***m_pppNormals;
	elem_t    ***m_sepNormals[3]; // To store normals separately

	//estimator for each X,Y and Z
	CBaseGradientComponentEstimator *m_gradEstimator[3];

	CVolume *m_pVolume;

	bool	m_bEpsilonNormalMode;
	double  m_Epsilon;


	std::string m_name;

	enum DATA_FORMAT{
		DATA_FORMAT_BYTE,
		DATA_FORMAT_FLOAT,
	};

	enum BOUNDARY_CONDITION{
		BOUNDARY_CONDITION_ZERO,
		BOUNDARY_CONDITION_PERIODIC,
		BOUNDARY_CONDITION_MIRROR
	};

	BOUNDARY_CONDITION  m_boundaryCondition;
	bool                m_bOnTheFlyNormals;
public:

	CBaseData( std::string name ) :  m_xSize(0), m_ySize(0), m_zSize(0),m_pppData(0), m_pppNormals(0),m_pVolume(0), m_bEpsilonNormalMode(0), m_Epsilon(DELTA), m_name(name), m_boundaryCondition( BOUNDARY_CONDITION_ZERO ),
	m_bOnTheFlyNormals(false)
	{
		m_gradEstimator[0] = NULL;
		m_gradEstimator[1] = NULL;
		m_gradEstimator[2] = NULL;
	}
	virtual ~CBaseData()
	{
		DestroyData();
		printf("Volume data destroyed\n");
	}

	virtual void Evaluate( const vector3_t min, const vector3_t max, const vector3_t step,
		elem_t (*fn)(const vector3_t p), void (*fnNormal)(const vector3_t p, vector3_t N) = 0 ) {};

	void CreateDataFromFile( const char* fname, DATA_FORMAT fmt, double xScale, double yScale, double zScale)
	{
		if( m_pppData ) {
			DestroyData();
		}

		switch( fmt ){
			case DATA_FORMAT_BYTE:
				m_pppData = LoadVUDFileByte( fname, m_xSize, m_ySize, m_zSize );
				break;
			case DATA_FORMAT_FLOAT:
				m_pppData = LoadVUDFileFloat( fname, m_xSize, m_ySize, m_zSize );
				break;
		}

		m_pVolume = new CVolume( m_xSize, m_ySize, m_zSize, xScale, yScale, zScale );

		// Creating normals memory

		m_pppNormals = new vector3_t** [ m_zSize + 2* GUARD_BAND];

		for(int z = 0; z< (m_zSize+2*GUARD_BAND);++z){

			m_pppNormals[z] = new vector3_t* [ m_ySize + 2* GUARD_BAND];

			for(int y=0;y<(m_ySize+2*GUARD_BAND);++y){
				m_pppNormals[z][y] = new vector3_t [ m_xSize + 2* GUARD_BAND ];
				memset(m_pppNormals[z][y], 0, sizeof(vector3_t) * (m_xSize+2*GUARD_BAND) );
			}
		}

		CalculateNormals();
	}

	void CreateData( int xSize, int ySize, int zSize, double xScale, double yScale, double zScale, atomic_t *pData )
	{

		if( m_pppData ) {
			DestroyData();
		}

		m_pppData = new elem_t ** [ zSize + 2*GUARD_BAND];
		for( int z = 0; z < (zSize+2*GUARD_BAND); z++ ){
			m_pppData[z] = new  elem_t *  [ ySize + 2*GUARD_BAND ];
			for( int y = 0; y< (ySize+2*GUARD_BAND); y++ ){
				//m_pppData[z][ySize - y - 1] = new double [ xSize ];

				m_pppData[z][y] = new elem_t [ xSize + 2*GUARD_BAND ];

				//memcpy( m_pppData[z][ ySize - y - 1 ], pOffset, sizeof( atomic_t ) * xSize );

                memset( m_pppData[z][y], 0, sizeof(elem_t) * (xSize + 2*GUARD_BAND) );

			}
		}


		atomic_t *pOffset = pData;

		for( int z = 0; z < zSize; z++ ){

			for( int y = 0; y< ySize; y++ ){

				/// TODO: Optimization possible in this loop
				for(int n = 0; n < xSize; ++n ){
					m_pppData[z + GUARD_BAND][y+GUARD_BAND][n+GUARD_BAND] = pOffset[n];
				}

				pOffset += xSize;
			}
		}


		m_xSize = xSize;
		m_ySize = ySize;
		m_zSize = zSize;


		m_pVolume = new CVolume( xSize, ySize, zSize, xScale, yScale, zScale );

		// Creating normals memory

		m_pppNormals = new vector3_t** [ zSize + 2*GUARD_BAND ];

		for(int z = 0; z< (zSize+2*GUARD_BAND);++z){

			m_pppNormals[z] = new vector3_t* [ ySize + 2*GUARD_BAND ];

			for(int y=0;y<(ySize+2*GUARD_BAND);++y){
				m_pppNormals[z][y] = new vector3_t [ xSize + 2*GUARD_BAND ];
				memset(m_pppNormals[z][y], 0, sizeof(vector3_t) * (xSize + 2*GUARD_BAND) );
			}
		}

		CalculateNormals();
	}

	void DestroyData()
	{
		if( m_pppData ) {
			for( int z = 0; z < (m_zSize+2*GUARD_BAND); z++ ){
				if( m_pppData[z] ){
					for( int y = 0; y < (m_ySize+2*GUARD_BAND); y++ ){
						if( m_pppData[z][y] ){
							DELETE_SAFE_ARRAY( m_pppData[z][y] );
						}
					}

					DELETE_SAFE_ARRAY( m_pppData[z] );
				}
			}

			DELETE_SAFE_ARRAY( m_pppData );
		};

		if( m_pppNormals ){

			for(int z=0;z< (m_zSize+2*GUARD_BAND);++z){

				if( m_pppNormals[z] ){

					for( int y = 0;y<(m_ySize+2*GUARD_BAND);++y){
						DELETE_SAFE_ARRAY( m_pppNormals[z][y] );
					}
				}

				DELETE_SAFE_ARRAY( m_pppNormals[z] );
			}

			DELETE_SAFE_ARRAY( m_pppNormals );
		}

		m_xSize = m_ySize = m_zSize = 0;

		DELETE_SAFE( m_pVolume );
	}

	inline CVolume* GetVolume(){ return m_pVolume; }

//	inline double FetchMirror( int x, int y, int z) const
//	{
//		int absXYZ[3] = { abs(x), abs(y), abs(z) };
//
//		if( x < 0 ) absXYZ[X] -= 1;
//		if( y < 0 ) absXYZ[Y] -= 1;
//		if( z < 0 ) absXYZ[Z] -= 1;
//
//		int odd[3] = { (absXYZ[X] / m_xSize) & 0x1, (absXYZ[Y] / m_ySize) & 0x1, (absXYZ[Z] / m_zSize ) & 0x1 };
//		int xyz[3];
//
//		if( odd[X] ){
//			xyz[X] = (m_xSize-1) - absXYZ[X] % m_xSize;
//		}
//		else{
//			xyz[X] = absXYZ[X] % m_xSize;
//		}
//
//		if( odd[Y] ){
//			xyz[Y] = (m_ySize-1) - absXYZ[Y] % m_ySize;
//		}
//		else{
//			xyz[Y] = absXYZ[Y] % m_ySize;
//		}
//
//		if( odd[Z] ){
//			xyz[Z] = (m_zSize-1) - absXYZ[Z] % m_zSize;
//		}
//		else{
//			xyz[Z] = absXYZ[Z] % m_zSize;
//		}
//
//		return m_pppData[xyz[Z]][xyz[Y]][xyz[X]];
//	}
//
//	inline elem_t FetchPeroidic( int x, int y, int z ) const
//	{
//		int xyz[3];
//		if( x > 0 ) xyz[X] = x % m_xSize;
//		else        xyz[X] = (m_xSize + x % m_xSize) % m_xSize;
//
//		if( y > 0 ) xyz[Y] = y % m_ySize;
//		else        xyz[Y] = (m_ySize + y % m_ySize) % m_ySize;
//
//		if( z > 0 ) xyz[Z] = z % m_zSize;
//		else        xyz[Z] = (m_zSize + z % m_zSize) % m_zSize;
//
//		return m_pppData[xyz[Z]][xyz[Y]][xyz[X]];
//	}

	inline elem_t Fetch( int X, int Y, int Z ) const
	{
//		if( (X >= m_xSize) || (X < 0) || (Y >= m_ySize) || Y < 0 || Z >= m_zSize  || Z < 0 )
//		{
//			//We have gone past the boundary in some axis
//			//I know this is a very lousy way to do things.. but hey
//			//who cares :)
//			switch( m_boundaryCondition ){
//				case BOUNDARY_CONDITION_ZERO:
//					return 0;
//					break;
//				case BOUNDARY_CONDITION_PERIODIC:
//					return FetchPeroidic( X,Y,Z );
//					break;
//				case BOUNDARY_CONDITION_MIRROR:
//					return FetchMirror( X, Y, Z );
//					break;
//			}
//		}

		return m_pppData[Z+GUARD_BAND][Y+GUARD_BAND][X+GUARD_BAND];
	}

	inline void FetchNormal( int X, int Y, int Z, vector3_t N ) const
	{
		vector3_t &v = m_pppNormals[Z+GUARD_BAND][Y+GUARD_BAND][X+GUARD_BAND];
		N[0] =  v[0];
		N[1] =  v[1];
		N[2] =  v[2];
	}

	inline void Set( int X, int Y, int Z, double v )
	{
		//assert( m_pppData );
		m_pppData[Z][Y][X] = v;
	}

	virtual void  CalculateNormals( void (*fnNormal)(const vector3_t p, vector3_t N) = 0 ){}
	virtual double Get( const vector3_t p ) const { return 0; }

	// The default GetNormal calculates the normal using numerical methods
	virtual void  GetNormal( const vector3_t p, vector3_t N ) const
	{
		vector3_t q[6];
		double delta2 = 2 * m_Epsilon;

		glmCopyVector3d( p, q[0]);
		glmCopyVector3d( p, q[1]);
		glmCopyVector3d( p, q[2]);
		glmCopyVector3d( p, q[3]);
		glmCopyVector3d( p, q[4]);
		glmCopyVector3d( p, q[5]);

		q[1][X] += m_Epsilon;
		q[0][X] -= m_Epsilon;

		q[3][Y] += m_Epsilon;
		q[2][Y] -= m_Epsilon;

		q[5][Z] += m_Epsilon;
		q[4][Z] -= m_Epsilon;

		N[X] = ( Get(q[1]) - Get(q[0]) ) / (delta2 );
		N[Y] = ( Get(q[3]) - Get(q[2]) ) / (delta2 );
		N[Z] = ( Get(q[5]) - Get(q[4]) ) / (delta2 );
	}

	virtual bool SetNormalEstimatorFilterByName( const char* name ) {return false; }
	virtual bool SetReconstructionFilterByName( const char* name ) { return false; }

	virtual bool SetNormalReconstructionFilterByName( const char *filter )
	{
		if( strcmp( filter, "EPSILON") == 0 ) {
			printf(" Gradient Estimation Set to EPSILON MODE\n");
			SetEpsilonNormalMode(true );
			return true;
		}
		return false;
	}

	virtual void HookGradientComponentEstimators()
	{ printf("CBaseData: Estimator not hooked!\n"); }; // must override this

	void SetBoundaryCondition( BOUNDARY_CONDITION cond ) { m_boundaryCondition = cond; }
	BOUNDARY_CONDITION GetBoundaryCondition( ) { return m_boundaryCondition; }
	void SetBoundaryCondition( const char *cond )
	{
		if( strcmp( cond, "BOUNDARY_CONDITION_ZERO") == 0 ) SetBoundaryCondition( BOUNDARY_CONDITION_ZERO );
		else if( strcmp( cond, "BOUNDARY_CONDITION_PERIODIC") == 0 ) SetBoundaryCondition( BOUNDARY_CONDITION_PERIODIC );
		else if( strcmp( cond, "BOUNDARY_CONDITION_MIRROR") == 0 ) SetBoundaryCondition( BOUNDARY_CONDITION_MIRROR );
		else{
			printf("CBaseData\n");
			printf("\tUnknown boundary condition %s\n",cond );
			printf("\tBoundary condtion is set to BOUNDARY_CONDITION_ZERO\n");
			SetBoundaryCondition( BOUNDARY_CONDITION_ZERO );
		}
	}
	void SetEpsilonNormalMode( bool mode ) { m_bEpsilonNormalMode = mode; };
	void SetEpsilon( double epsilon ) { m_Epsilon = epsilon; }
	void SetName( const char* name ) { m_name = std::string(name) ; }
	const char* GetName(){ return m_name.c_str(); }

	friend class CBaseGradientComponentEstimator;

};


#endif

