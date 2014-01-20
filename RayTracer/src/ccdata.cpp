#include <stdio.h>

#include "ccdata.h"
#include "camera.h"

double CubicBSpline( double t )
{
	#define ONE_OVER_6 0.166666666667f

	if( t <= -2 ) return 0;
	if( t >= 2  ) return 0;

	if( t > -2 && t <= -1 ){
		double tplus2 = t+2;
		return ONE_OVER_6 * tplus2*tplus2*tplus2;
	}

	if( t > -1 && t <= 0 ){
		double tplus2 = t + 2;
		double tplus1 = t + 1;
		return ONE_OVER_6 * ( tplus2 * tplus2 * tplus2 - 4 * tplus1 * tplus1 * tplus1 );
	}

	if( t > 0 && t <= 1 ){
		double tplus2inv = 2 - t;
		double tplus1inv = 1 - t;
		return ONE_OVER_6 * ( tplus2inv * tplus2inv * tplus2inv - 4 * tplus1inv * tplus1inv * tplus1inv );
	}

	if( t > 1 && t < 2 ){
		double tplus2inv = 2-t;
		return ONE_OVER_6 * tplus2inv*tplus2inv*tplus2inv;
	}

	return 0;
}
float weights[4][4] =
{
	{ 1.0f/6.0f, -3.0f/6.0f,  3.0f/6.0f, -1.0f/6.0f },
	{ 4.0f/6.0f,  0.0f/6.0f, -6.0f/6.0f,  3.0f/6.0f },
	{ 1.0f/6.0f,  3.0f/6.0f,  3.0f/6.0f, -3.0f/6.0f },
	{ 0.0f/6.0f,  0.0f/6.0f,  0.0f/6.0f,  1.0f/6.0f }
};


inline void CubicBSplineComputeWeights( double t, double w[])
{
#define one_over_six 0.166666667
#define four_over_six 0.666666667

	double t2 = t*t;	
 	double t3 = t2*t;
 	w[0] = 0.5 * (t2 -  t) + one_over_six*(1-t3);
	w[1] = 0.5* t3 - t2 + four_over_six;
	w[2] = 0.5*(t+t2-t3) + one_over_six;
	w[3] = t3 * one_over_six;
}

inline elem_t CubicBSplineInterpolateWeights( double w[], double p[] )
{
	return (p[0] *w[0]) + (p[1] *w[1]) + (p[2] *w[2]) + (p[3] *w[3]);
}


// inline elem_t CubicBSplineInterpolate(double t, double p0,double p1, double p2, double p3) 
// {
// 	return (1.0f * (weights[0][0] * p0 + weights[1][0] * p1 + weights[2][0]
// 			* p2 )) +
// 			(t * (weights[0][1] * p0 + weights[2][1] *
// 			p2 )) +
// 			(t * t * (weights[0][2] * p0 + weights[1][2] * p1 +
// 			weights[2][2] * p2 )) +
// 			(t * t * t * (weights[0][3] * p0 + weights[1][3] * p1 +
// 			weights[2][3] * p2 + weights[3][3] * p3));
// }

inline elem_t CubicBSplineInterpolate(double t[3], double p0,double p1, double p2, double p3) 
{
	return (1.0f * (weights[0][0] * p0 + weights[1][0] * p1 + weights[2][0]
			* p2 )) +
			(t[0] * (weights[0][1] * p0 + weights[2][1] *
			p2 )) +
			(t[1] * (weights[0][2] * p0 + weights[1][2] * p1 +
			weights[2][2] * p2 )) +
			(t[2] * (weights[0][3] * p0 + weights[1][3] * p1 +
			weights[2][3] * p2 + weights[3][3] * p3));
}

					     
std::map< std::string, CCCData::RECONSTRUCTION_FILTER >  CCCData::m_reconstructionFilterNameMap;
std::map< std::string, CCCData::NORMAL_DISCRETE_FILTER > CCCData::m_normalFilterNameMap;

bool CCCData::m_bMapInitialized = false;

void CCCData::InitializeMaps()
{
	if( m_bMapInitialized ) return;

	MAP_ENTRY_BY_FLAG( m_reconstructionFilterNameMap, CC_NEAREST );
	MAP_ENTRY_BY_FLAG( m_reconstructionFilterNameMap, CC_TRILINEAR );
	MAP_ENTRY_BY_FLAG( m_reconstructionFilterNameMap, CC_TRICUBICBSPLINE );

	MAP_ENTRY_BY_FLAG( m_normalFilterNameMap, CC_NORMAL_CD_10TH);
	MAP_ENTRY_BY_FLAG( m_normalFilterNameMap, CC_NORMAL_CD_6TH );
	MAP_ENTRY_BY_FLAG( m_normalFilterNameMap, CC_NORMAL_CD_4TH );
	MAP_ENTRY_BY_FLAG( m_normalFilterNameMap, CC_NORMAL_CD_2ND );

	m_bMapInitialized = true;
}

CCCData::CCCData( const char *name, const char *normalEstName, bool bOnTheFlyNormals, const char *gridPointDumpFile) : CBaseData( name )
{
	m_normalEst = CC_NORMAL_CD_4TH;
	m_reconFilter = m_normalReconFilter = CC_TRICUBICBSPLINE;
	m_bOnTheFlyNormals = bOnTheFlyNormals;

	m_xSize = m_ySize = m_zSize = 0;
	m_pppData = 0;
	m_gridPointDumpFile = gridPointDumpFile;

	SetNormalEstimatorFilterByName( normalEstName );

	m_invScale[X] = 1.0;
	m_invScale[Y] = 1.0;
	m_invScale[Z] = 1.0;

}

CCCData::CCCData( const char* name, const char* fname, DATA_FORMAT fmt, double xScale, double yScale, double zScale, const char* normalEstName, bool bOnTheFlyNormals, const char *gridPointDumpFile ) : CBaseData( name )
{
	m_normalEst = CC_NORMAL_CD_4TH;
	m_reconFilter = m_normalReconFilter = CC_TRICUBICBSPLINE;
	m_bOnTheFlyNormals = bOnTheFlyNormals;

	m_xSize = m_ySize = m_zSize = 0;
	m_pppData = 0;

	m_gridPointDumpFile = gridPointDumpFile;

	SetNormalEstimatorFilterByName( normalEstName );

	switch( fmt ){
		case DATA_FORMAT_BYTE:
			m_pppData  = LoadVUDFileByte( fname, m_xSize, m_ySize, m_zSize );
			break;
		case DATA_FORMAT_FLOAT:
			m_pppData  = LoadVUDFileFloat( fname, m_xSize, m_ySize, m_zSize );
			break;
	}

	m_pVolume = new CVolume( m_xSize, m_ySize, m_zSize, xScale, yScale, zScale );

	// Creating normals memory

	m_pppNormals = new vector3_t** [ m_zSize + 2*GUARD_BAND];

	for(int z = 0; z< (m_zSize+2*GUARD_BAND);++z){

		m_pppNormals[z] = new vector3_t* [ m_ySize + 2*GUARD_BAND];

		for(int y=0;y<(m_ySize+2*GUARD_BAND);++y){
			m_pppNormals[z][y] = new vector3_t [ m_xSize + 2*GUARD_BAND];
			memset(m_pppNormals[z][y], 0, sizeof(vector3_t) * (m_xSize+2*GUARD_BAND) );
		}
	}

	m_invScale[X] = 1.0/GetVolume()->m_xScale;
	m_invScale[Y] = 1.0/GetVolume()->m_yScale;
	m_invScale[Z] = 1.0/GetVolume()->m_zScale;

	CalculateNormals();
}

CCCData::CCCData( const char* name, const char* fname, DATA_FORMAT fmt, const char* xname, const char *yname, const char *zname, double xScale, double yScale, double zScale,  const char *gridPointDumpFile ) : CBaseData( name )
{
	m_normalEst = CC_NORMAL_CD_4TH;
	m_reconFilter = m_normalReconFilter = CC_TRICUBICBSPLINE;
	m_bOnTheFlyNormals = false;

	m_xSize = m_ySize = m_zSize = 0;
	m_pppData = 0;
	m_gridPointDumpFile = gridPointDumpFile;

	int dummyX, dummyY, dummyZ;
	elem_t ***pppXNormal = NULL, ***pppYNormal = NULL, ***pppZNormal = NULL;

	switch( fmt ){
		case DATA_FORMAT_BYTE:
			m_pppData  = LoadVUDFileByte( fname, m_xSize, m_ySize, m_zSize );

			pppXNormal = LoadVUDFileByte( xname, dummyX, dummyY, dummyZ );
			pppYNormal = LoadVUDFileByte( yname, dummyX, dummyY, dummyZ );
			pppZNormal = LoadVUDFileByte( zname, dummyX, dummyY, dummyZ );

			break;
		case DATA_FORMAT_FLOAT:
			m_pppData  = LoadVUDFileFloat( fname, m_xSize, m_ySize, m_zSize );

			pppXNormal = LoadVUDFileFloat( xname, dummyX, dummyY, dummyZ );
			pppYNormal = LoadVUDFileFloat( yname, dummyX, dummyY, dummyZ );
			pppZNormal = LoadVUDFileFloat( zname, dummyX, dummyY, dummyZ );

			break;
	}

	m_pVolume = new CVolume( m_xSize, m_ySize, m_zSize, xScale, yScale, zScale );

	// Creating normals memory

	m_pppNormals = new vector3_t** [ m_zSize + 2*GUARD_BAND ];

	for(int z = 0; z< (m_zSize+2*GUARD_BAND);++z){

		m_pppNormals[z] = new vector3_t* [ m_ySize +  2*GUARD_BAND];

		for(int y=0;y<(m_ySize+2*GUARD_BAND);++y){
			m_pppNormals[z][y] = new vector3_t [ m_xSize + 2*GUARD_BAND ];

            memset(m_pppNormals[z][y], 0, sizeof( vector3_t) * (m_xSize +2*GUARD_BAND));

		}
	}

	m_sepNormals[X] = new elem_t ** [m_zSize+2*GUARD_BAND];
	m_sepNormals[Y] = new elem_t ** [m_zSize+2*GUARD_BAND];
	m_sepNormals[Z] = new elem_t ** [m_zSize+2*GUARD_BAND];

	for(int z=0;z<(m_zSize+2*GUARD_BAND);++z){

		m_sepNormals[X][z] = new elem_t * [ m_ySize +2*GUARD_BAND];
		m_sepNormals[Y][z] = new elem_t * [ m_ySize +2*GUARD_BAND];
		m_sepNormals[Z][z] = new elem_t * [ m_ySize +2*GUARD_BAND];

		for(int y=0;y<(m_ySize+2*GUARD_BAND);++y){

			m_sepNormals[X][z][y] = new elem_t [ m_xSize +2*GUARD_BAND];
			m_sepNormals[Y][z][y] = new elem_t [ m_xSize +2*GUARD_BAND];
			m_sepNormals[Z][z][y] = new elem_t [ m_xSize +2*GUARD_BAND];

			memset( m_sepNormals[X][z][y], 0, sizeof( elem_t) * (m_xSize +2*GUARD_BAND) );
			memset( m_sepNormals[Y][z][y], 0, sizeof( elem_t) * (m_xSize +2*GUARD_BAND) );
			memset( m_sepNormals[Z][z][y], 0, sizeof( elem_t) * (m_xSize +2*GUARD_BAND) );
		}
	}
	
    for(int z = 0; z< m_zSize;++z){
		for(int y=0;y<m_ySize;++y){
			for(int x=0;x<m_xSize;++x){
				m_pppNormals[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND][X] = pppXNormal[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND];								
				m_pppNormals[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND][Y] = pppYNormal[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND];
				m_pppNormals[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND][Z] = pppZNormal[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND];
			}
		}
	}

	for(int z = 0; z< m_zSize;++z){
		for(int y=0;y<m_ySize;++y){
			for(int x=0;x<m_xSize;++x){
				m_sepNormals[X][z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND] = m_pppNormals[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND][X];
				
				m_sepNormals[Y][z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND] = m_pppNormals[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND][Y];
				
				m_sepNormals[Z][z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND] = m_pppNormals[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND][Z];
			}
		}
	}
	
	// Destroying all the normals
	if( pppXNormal ) {
		for( int z = 0; z < m_zSize; z++ ){
			if( pppXNormal[z] ){
				for( int y = 0; y < m_ySize; y++ ){
					if( pppXNormal[z][y] ){
						DELETE_SAFE_ARRAY( pppXNormal[z][y] );
					}
				}

				DELETE_SAFE_ARRAY( pppXNormal[z] );
			}
		}

		DELETE_SAFE_ARRAY( pppXNormal );
	};

	if( pppYNormal ) {
		for( int z = 0; z < m_zSize; z++ ){
			if( pppYNormal[z] ){
				for( int y = 0; y < m_ySize; y++ ){
					if( pppYNormal[z][y] ){
						DELETE_SAFE_ARRAY( pppYNormal[z][y] );
					}
				}

				DELETE_SAFE_ARRAY( pppYNormal[z] );
			}
		}

		DELETE_SAFE_ARRAY( pppYNormal );
	};

	if( pppZNormal ) {
		for( int z = 0; z < m_zSize; z++ ){
			if( pppZNormal[z] ){
				for( int y = 0; y < m_ySize; y++ ){
					if( pppZNormal[z][y] ){
						DELETE_SAFE_ARRAY( pppZNormal[z][y] );
					}
				}

				DELETE_SAFE_ARRAY( pppZNormal[z] );
			}
		}

		DELETE_SAFE_ARRAY( pppZNormal );
	};

	m_invScale[X] = 1.0/GetVolume()->m_xScale;
	m_invScale[Y] = 1.0/GetVolume()->m_yScale;
	m_invScale[Z] = 1.0/GetVolume()->m_zScale;	
}



CCCData::CCCData( const char* name, const char *fname, double xScale, double yScale, double zScale ) : CBaseData( std::string(name) )
{
	m_normalEst = CC_NORMAL_CD_4TH;
	m_reconFilter = m_normalReconFilter = CC_TRICUBICBSPLINE;

	m_xSize = m_ySize = m_zSize = 0;
	m_pppData = 0;

	int xSize,ySize,zSize;
	FILE *fp = fopen( fname, "rb" );

	if( !fp ) {
		printf("File \"%s\" cannot be opened!\n",fname );
		return;
	}

	char line[256];
	char *s;
	size_t readsize;

	for(int i=0;i<6;++i){
		s= fgets(line, 256, fp ); // reads the first commented line
	}

	char dummy[20];
	//line should have dimension information now
	sscanf(line,"%s %d %d %d", dummy, &xSize, &ySize, &zSize );

	//now just read off the next 7 lines
	for(int i=0;i<7;++i){
		s = fgets(line, 256, fp ); // reads the first commented line
	}


//	fscanf(fp,"%d %d %d",&xSize,&ySize,&zSize);

	printf("Reading %s (%dx%dx%d).... ", fname, xSize, ySize, zSize );

	int totalAtomic = xSize * ySize * zSize;
	atomic_t *pMem = new atomic_t [ totalAtomic ];

	readsize = fread(pMem, totalAtomic * sizeof(atomic_t),1,fp);

	CreateData( xSize, ySize, zSize, xScale, yScale, zScale, pMem );

	DELETE_SAFE_ARRAY( pMem );

	printf("Done!\n");

	m_invScale[X] = 1.0/GetVolume()->m_xScale;
	m_invScale[Y] = 1.0/GetVolume()->m_yScale;
	m_invScale[Z] = 1.0/GetVolume()->m_zScale;

}

void CCCData::Evaluate( const vector3_t min, const vector3_t max, const vector3_t step,
		elem_t (*fn)(const vector3_t p), void (*fnNormal)(const vector3_t p, vector3_t N))
{

	DestroyData();

	m_xSize = (int)((max[X] - min[X]) / step[X])+1;
	m_ySize = (int)((max[Y] - min[Y]) / step[Y])+1;
	m_zSize = (int)((max[Z] - min[Z]) / step[Z])+1;


	m_pppData = new elem_t ** [ m_zSize + 2*GUARD_BAND ];
	for( int z = 0; z < (m_zSize + 2*GUARD_BAND); z++ ){
		m_pppData[z] = new  elem_t *  [ m_ySize + 2*GUARD_BAND ];
		for( int y = 0; y< (m_ySize + 2*GUARD_BAND); y++ ){
			m_pppData[z][y] = new elem_t [ m_xSize + 2*GUARD_BAND ];
			memset( m_pppData[z][y], 0, sizeof(elem_t) * (m_xSize + 2*GUARD_BAND) );
		}
	}



	m_pVolume = new CVolume( m_xSize, m_ySize, m_zSize, step[X], step[Y], step[Z] );


	for(int z=0;z<m_zSize;++z){
		for(int y=0;y<m_ySize;++y){
			for(int x=0;x<m_xSize;++x){

				vector3_t p;

				p[X] = x*step[X] + min[X];
				p[Y] = y*step[Y] + min[Y];
				p[Z] = z*step[Z] + min[Z];

				m_pppData[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND] = fn( p );



			}
		}
	}

	// Creating normals memory

	m_pppNormals = new vector3_t** [ m_zSize + 2*GUARD_BAND ];

	for(int z = 0; z< (m_zSize+2*GUARD_BAND);++z){

		m_pppNormals[z] = new vector3_t* [ m_ySize + 2*GUARD_BAND ];

		for(int y=0;y<(m_ySize+2*GUARD_BAND);++y){
			m_pppNormals[z][y] = new vector3_t [ m_xSize + 2*GUARD_BAND];
			memset(m_pppNormals[z][y], 0, sizeof(vector3_t) * (m_xSize + 2*GUARD_BAND) );
		}
	}

	m_invScale[X] = 1.0/GetVolume()->m_xScale;
	m_invScale[Y] = 1.0/GetVolume()->m_yScale;
	m_invScale[Z] = 1.0/GetVolume()->m_zScale;


	CalculateNormals( fnNormal );


}

bool CCCData::SetNormalEstimatorFilterByName( const char* filter )
{
	if( CBaseData::SetNormalEstimatorFilterByName( filter ) ) return true;

	bool ret = true;

	InitializeMaps();

	if( m_normalFilterNameMap.find( std::string(filter) ) == m_normalFilterNameMap.end() ){
		printf(" CC Lattice ( %s )\n", GetName() );
		printf("\tThe normal discrete filter %s doesn't exists\n", filter);
		printf("\tNormal discrete filter will be set to default\n");

		ret = false;
	}

	m_normalEst = m_normalFilterNameMap[std::string(filter)];

	return ret;
}

bool CCCData::SetReconstructionFilterByName( const char* filter )
{
	if( CBaseData::SetReconstructionFilterByName( filter ) ){
		return true;
	}

	bool ret = true;

	InitializeMaps();

	if( m_reconstructionFilterNameMap.find( std::string(filter) ) == m_reconstructionFilterNameMap.end() ){
		printf(" CC Lattice ( %s )\n", GetName() );
		printf("\tThe reconstruction filter %s doesn't exists\n", filter);
		printf("\tReconstruction filter will be set to default\n");

		ret = false;
	}

	m_reconFilter = m_reconstructionFilterNameMap[std::string(filter)];

	return ret;
}

bool CCCData::SetNormalReconstructionFilterByName( const char *filter )
{

	if( CBaseData::SetNormalReconstructionFilterByName( filter )) {
		return true;
	}

	bool ret = true;

	InitializeMaps();

	if( m_reconstructionFilterNameMap.find( std::string(filter) ) == m_reconstructionFilterNameMap.end() ){
		printf(" CC Lattice ( %s )\n", GetName() );
		printf("\tThe normal reconstruction filter %s doesn't exists\n", filter);
		printf("\tNormal reconstruction filter will be set to default\n");

		ret = false;
	}

	m_normalReconFilter = m_reconstructionFilterNameMap[std::string(filter)];

	return ret;
}

elem_t CCCData::ComputeNormalXCentralDifference2nd( int x,int y, int z ) const
{
	return 0.5 * m_invScale[X] * (Fetch( x + 1, y, z ) - Fetch( x - 1, y, z ));
}

elem_t CCCData::ComputeNormalYCentralDifference2nd( int x,int y, int z ) const
{
	return 0.5 * m_invScale[Y] * (Fetch( x, y+1, z ) - Fetch( x, y-1, z ));
}

elem_t CCCData::ComputeNormalZCentralDifference2nd( int x,int y, int z ) const
{
	return 0.5 * m_invScale[Z] * (Fetch( x, y, z+1 ) - Fetch( x, y, z-1 ));
}

void   CCCData::PrecompNormalCentralDifference2nd( vector3_t N, int x,int y, int z )
{
	int offset[2]    = { 1, -1 };
	elem_t filter[2] = {1.0/2.0, -1/2.0};

	N[X] = (Fetch( x + offset[0], y, z ) * filter[0] + Fetch( x + offset[1], y, z ) * filter[1]) / GetVolume()->m_xScale;
	N[Y] = (Fetch( x, y + offset[0], z ) * filter[0] + Fetch( x, y + offset[1], z ) * filter[1]) / GetVolume()->m_yScale;
	N[Z] = (Fetch( x, y, z + offset[0] ) * filter[0] + Fetch( x, y, z + offset[1] ) * filter[1]) / GetVolume()->m_zScale;
}

elem_t CCCData::ComputeNormalXCentralDifference4th( int x, int y, int z ) const
{
	return m_invScale[X] * ((0.666666667) * (Fetch( x + 1, y, z ) - Fetch( x - 1, y, z ))
			+
			(0.083333333) * (Fetch( x - 2, y, z ) - Fetch( x + 2, y, z )));
}

elem_t CCCData::ComputeNormalYCentralDifference4th( int x, int y, int z ) const
{
	return m_invScale[Y] * ((0.666666667) * (Fetch( x, y+1, z ) - Fetch( x, y-1, z ))
			+
			(0.083333333) * (Fetch( x, y-2, z ) - Fetch( x, y+2, z )));
}
elem_t CCCData::ComputeNormalZCentralDifference4th( int x, int y, int z ) const
{
	return m_invScale[Z] * ((0.666666667) * (Fetch( x, y, z +1 ) - Fetch( x, y, z-1 ))
			+
			(0.083333333) * (Fetch( x, y, z-2 ) - Fetch( x, y, z+2 )));
}

void   CCCData::PrecompNormalCentralDifference4th( vector3_t N, int x,int y, int z )
{
	int offset[4]    = { 1, -1, 2, -2 };
	elem_t filter[4] = { 2.0/3.0, -2.0/3.0, -1.0/12.0, 1.0/12.0 };

	N[X] = N[Y] = N[Z] = 0;

	for(int i = 0; i < 4; ++i ){
		N[X] += Fetch( x + offset[i], y, z ) * filter[i];
		N[Y] += Fetch( x, y + offset[i], z ) * filter[i];
		N[Z] += Fetch( x, y, z + offset[i] ) * filter[i];
	}

	N[X] /= GetVolume()->m_xScale;
	N[Y] /= GetVolume()->m_yScale;
	N[Z] /= GetVolume()->m_zScale;
}
void   CCCData::PrecompNormalCentralDifference6th( vector3_t N, int x,int y, int z )
{
	int offset[6]    = { 1, -1, 2, -2, 3,-3 };
	elem_t filter[6] = { 3.0/4.0, -3.0/4.0, -3.0/20.0, 3.0/20.0, 1.0/60.0, -1.0/60.0 };

	N[X] = N[Y] = N[Z] = 0;

	for(int i = 0; i < 6; ++i ){
		N[X] += Fetch( x + offset[i], y, z ) * filter[i];
		N[Y] += Fetch( x, y + offset[i], z ) * filter[i];
		N[Z] += Fetch( x, y, z + offset[i] ) * filter[i];
	}

	N[X] /= GetVolume()->m_xScale;
	N[Y] /= GetVolume()->m_yScale;
	N[Z] /= GetVolume()->m_zScale;

}

void CCCData::PrecompNormalCentralDifference10th( vector3_t N, int x, int y, int z )
{
	int offset[10]    = { 1, -1, 2, -2, 3,-3,4,-4, 5,-5 };
	elem_t filter[10] = { 5.0/6.0, -5.0/6.0, -5.0/21.0, 5.0/21.0, 5.0/84.0, -5.0/84.0, -5.0/504.0, 5.0/504.0, 1.0/1260.0, -1.0/1260.0 };

	N[X] = N[Y] = N[Z] = 0;

	for(int i = 0; i < 10; ++i ){
		N[X] += Fetch( x + offset[i], y, z ) * filter[i];
		N[Y] += Fetch( x, y + offset[i], z ) * filter[i];
		N[Z] += Fetch( x, y, z + offset[i] ) * filter[i];
	}

	N[X] /= GetVolume()->m_xScale;
	N[Y] /= GetVolume()->m_yScale;
	N[Z] /= GetVolume()->m_zScale;
}


void CCCData::CalculateNormals( void (*fnNormal)(const vector3_t p, vector3_t N) )
{
	//double d[6] = {0};

	if( m_bOnTheFlyNormals || m_bEpsilonNormalMode){
		printf("CCCData: Skipping precomputation of gradients\n");
		m_sepNormals[X] = NULL;
		m_sepNormals[Y] = NULL;
		m_sepNormals[Z] = NULL;
		return;
	}

	FILE *fp = NULL;

	if( fnNormal ){
		fp = fopen( "CCData_Normal_Data.txt","w");
	}

	m_sepNormals[X] = new elem_t ** [m_zSize+2*GUARD_BAND];
	m_sepNormals[Y] = new elem_t ** [m_zSize+2*GUARD_BAND];
	m_sepNormals[Z] = new elem_t ** [m_zSize+2*GUARD_BAND];

	for(int z=0;z<(m_zSize+2*GUARD_BAND);++z){

		m_sepNormals[X][z] = new elem_t * [ m_ySize +2*GUARD_BAND];
		m_sepNormals[Y][z] = new elem_t * [ m_ySize +2*GUARD_BAND];
		m_sepNormals[Z][z] = new elem_t * [ m_ySize +2*GUARD_BAND];

		for(int y=0;y<(m_ySize+2*GUARD_BAND);++y){

			m_sepNormals[X][z][y] = new elem_t [ m_xSize +2*GUARD_BAND];
			m_sepNormals[Y][z][y] = new elem_t [ m_xSize +2*GUARD_BAND];
			m_sepNormals[Z][z][y] = new elem_t [ m_xSize +2*GUARD_BAND];

			memset( m_sepNormals[X][z][y], 0, sizeof( elem_t) * (m_xSize +2*GUARD_BAND) );
			memset( m_sepNormals[Y][z][y], 0, sizeof( elem_t) * (m_xSize +2*GUARD_BAND) );
			memset( m_sepNormals[Z][z][y], 0, sizeof( elem_t) * (m_xSize +2*GUARD_BAND) );
		}
	}

	for(int z=0;z<m_zSize;++z){
		for(int y=0;y<m_ySize;++y){
			for(int x=0;x<m_xSize;++x){

				/*
				if( (x+1) < m_xSize)	d[0] = (double)Fetch( x+1, y, z );
				if( (y+1) < m_ySize)	d[1] = (double)Fetch( x, y+1, z );
				if( (x-1) >= 0) 		d[2] = (double)Fetch( x-1,y,  z );
				if( (y-1) >= 0 )		d[3] = (double)Fetch( x, y-1, z );
				if( (z-1) >= 0 )		d[4] = (double)Fetch( x, y, z-1 );
				if( (z+1) < m_zSize)	d[5] = (double)Fetch( x, y, z+1 );

				vector3_t &v = m_pppNormals[z][y][x];



				v[X] = (d[0] - d[2]) / (2.0f * GetVolume()->m_xScale);
				v[Y] = (d[1] - d[3]) / (2.0f * GetVolume()->m_yScale);
				v[Z] = (d[5] - d[4]) / (2.0f * GetVolume()->m_zScale);

				//NOTE: Dont normalize
				//glmNormalizeVector3d( v );
				*/

				vector3_t &v = m_pppNormals[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND];

				switch( m_normalEst){
					case CC_NORMAL_CD_2ND:
						PrecompNormalCentralDifference2nd( v, x,y,z );
						break;
					case CC_NORMAL_CD_4TH:
						PrecompNormalCentralDifference4th( v, x,y,z );
						break;
					case CC_NORMAL_CD_6TH:
						PrecompNormalCentralDifference6th( v, x,y,z );
						break;
					case CC_NORMAL_CD_10TH:
						PrecompNormalCentralDifference10th( v, x, y, z );
						break;
				}

				if( fnNormal ){
					vector2_t screenSpace;
					vector3_t worldSpace;
					vector3_t gridSpace;

					gridSpace[X] = x;
					gridSpace[Y] = y;
					gridSpace[Z] = z;

					GetVolume()->FindObjectSpacePoint( gridSpace, worldSpace );
					camGetScreenSpacePointFromWorldSpace( screenSpace, worldSpace );

					vector3_t N;

					fnNormal( worldSpace, N );

					//glmNormalizeVector3d( N );

					fprintf(fp,"%d %d %d; %d %d; %f %f %f; %f %f %f;\n",x,y,z, (int)(screenSpace[X]+0.5),(int)(screenSpace[Y]+0.5), N[X], N[Y], N[Z], v[X],v[Y],v[Z]);

				}

				m_sepNormals[X][z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND] = v[X];
				m_sepNormals[Y][z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND] = v[Y];
				m_sepNormals[Z][z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND] = v[Z];
			}
		}
	}

	if( fnNormal ){
		fclose(fp);
	}
}

double CCCData::GetNearest( const vector3_t p ) const
{
	return m_pppData[(int)(p[2]+0.5f)+GUARD_BAND][(int)(p[1]+0.5f)+GUARD_BAND][(int)(p[0]+0.5f)+GUARD_BAND];
}

double CCCData::GetTriLinear( const vector3_t v ) const
{
	double p[8] = {0};

	int floorX = ((int)v[0]);
	int nextX  = floorX + 1;
	int floorY = ((int)v[1]);
	int nextY  = floorY+1;
	int floorZ = ((int)v[2]);
	int nextZ  = floorZ+1;


	bool nextXWithinBound = nextX < m_xSize;
	bool nextYWithinBound = nextY < m_ySize;
	bool nextZWithinBound = nextZ < m_zSize;

	if( nextZWithinBound ){
		p[0] = Fetch( floorX, floorY, nextZ );

		if( nextXWithinBound ) {

			p[1] = Fetch( nextX, floorY, nextZ );

			if( nextYWithinBound ){
				p[2] = Fetch(nextX, nextY,nextZ);
			}
		}

		if( nextYWithinBound ){
			p[3] = Fetch( floorX, nextY, nextZ );
		}

	}

	p[4] = Fetch( floorX, floorY, floorZ );

	if( nextXWithinBound ){
		p[5] = Fetch( nextX, floorY, floorZ );

		if( nextYWithinBound ){
			p[6] = Fetch( nextX, nextY, floorZ );
		}
	}

	if( nextYWithinBound ){
		p[7] = Fetch(floorX, nextY, floorZ );
	}

	double dx = v[0] - floorX;
	double dy = v[1] - floorY;
	double dz = v[2] - floorZ;

	double invDx = 1 - dx;
	double invDy = 1 - dy;
	double invDz = 1 - dz;

	return (( (p[0]*invDx + p[1]*dx) * dz ) + ( (p[4]*invDx + p[5]*dx) * invDz )) * invDy +

	       (( (p[3]*invDx + p[2]*dx) * dz ) + ( (p[7]*invDx + p[6]*dx) * invDz )) * dy;
}

// float CCSphereFlow::interpolate( const double pos[3], int gridId ) {
// 	float x = pos[0], y = pos[1], z = pos[2];
// 	int x2, y2, z2;
// 	float dx2, dy2, dz2;
// 
//    x1 = MYROUND(x + 1); dx1 = -(x - x1);
//    y1 = MYROUND(y + 1); dy1 = -(y - y1);
//    z1 = MYROUND(z + 1); dz1 = -(z - z1);
// 
// 
// 	x2 = (int)floor(x); dx2 = x - x2;
// 	y2 = (int)floor(y); dy2 = y - y2;
// 	z2 = (int)floor(z); dz2 = z - z2;
// 
//   assert(dx1>=0 && dy1 >=0 && dz1 >=0);
// 	assert(dx2>=0 && dy2 >=0 && dz2 >=0);
// 
// 	double pz[4];
// 	double py[4];
// 	for(int z = z2 - 1; z < z2 + 3; z++) {
// 		for(int y = y2 - 1; y < y2 + 3; y++) {
// 			double p0, p1, p2, p3;
// 			int index3D[3];
// 			ccToIndex( x2-1,y,z, index3D );
// 			p0 = getValue( index3D, gridId  );
// 			ccToIndex( x2, y, z, index3D );
// 			p1 = getValue( index3D, gridId );
// 			ccToIndex( x2+1,y,z, index3D );
// 			p2 = getValue( index3D, gridId );
// 			ccToIndex(x2+2, y, z, index3D );
// 			p3 = getValue( index3D, gridId );
// 
// 			py[y - y2 + 1] =
// 					cubicInterpolate(dx2, p0, p1, p2, p3);
// 		}
// 		pz[z - z2 + 1] = cubicInterpolate(dy2, py[0], py[1], py[2], py[3]);
// 	}
// 	return cubicInterpolate(dz2, pz[0], pz[1], pz[2], pz[3]);
// }

// inline elem_t CubicBSplineInterpolate(float t, double p0,double p1, double p2, double p3)
		
elem_t CCCData::GetTriCubicBSplineTurbo( const vector3_t p ) const
{
	int floorX = (int)p[0];
	int floorY = (int)p[1];
	int floorZ = (int)p[2];
	
	int xNode[4] = {
		floorX - 1,
		floorX,
		floorX + 1,
		floorX + 2,
	};

	int yNode[4] = {
		floorY - 1,
		floorY,
		floorY + 1,
		floorY + 2,
	};

	int zNode[4] = {
		floorZ - 1,
		floorZ,
		floorZ + 1,
		floorZ + 2,
	};
		
	elem_t dx = p[0] - floorX;
	elem_t dy = p[1] - floorY;
	elem_t dz = p[2] - floorZ;
	
	elem_t nx[4],ny[4],nz[4];
	//double t[3];
	double w[3][4];
	
	CubicBSplineComputeWeights(dx, w[0]);
	CubicBSplineComputeWeights(dy, w[1]);
	CubicBSplineComputeWeights(dz, w[2]);
	
	for( int z=0;z<4;z++){
		for(int y=0;y<4;y++){
			
			for( int x=0;x<4;x++){
				nx[x] = Fetch( xNode[x], yNode[y], zNode[z] );
			}
			
			ny[y] = CubicBSplineInterpolateWeights( w[0], nx );
		}
		
		nz[z] = CubicBSplineInterpolateWeights( w[1], ny );
	}		
	return CubicBSplineInterpolateWeights( w[2], nz );
}

void CCCData::GetNormalTriCubicBSplineTurbo( const vector3_t p, vector3_t N ) const
{
	int floorX = (int)p[0];
	int floorY = (int)p[1];
	int floorZ = (int)p[2];
	
	int xNode[4] = {
		floorX - 1,
		floorX,
		floorX + 1,
		floorX + 2,
	};

	int yNode[4] = {
		floorY - 1,
		floorY,
		floorY + 1,
		floorY + 2,
	};

	int zNode[4] = {
		floorZ - 1,
		floorZ,
		floorZ + 1,
		floorZ + 2,
	};
		
	elem_t dx = p[0] - floorX;
	elem_t dy = p[1] - floorY;
	elem_t dz = p[2] - floorZ;
	
	elem_t nx[3][4],ny[3][4],nz[3][4];
	//double t[3];
	double w[3][4];
	
	CubicBSplineComputeWeights(dx, w[0]);
	CubicBSplineComputeWeights(dy, w[1]);
	CubicBSplineComputeWeights(dz, w[2]);
	
	for( int z=0;z<4;z++){
		for(int y=0;y<4;y++){
			
			for( int x=0;x<4;x++){
				int index[3] = {xNode[x], yNode[y], zNode[z]};
				nx[0][x] = GetValueNorm( index, 0 );
				nx[1][x] = GetValueNorm( index, 1 );
				nx[2][x] = GetValueNorm( index, 2 );
			}
			
			ny[0][y] = CubicBSplineInterpolateWeights( w[0], nx[0] );
			ny[1][y] = CubicBSplineInterpolateWeights( w[0], nx[1] );
			ny[2][y] = CubicBSplineInterpolateWeights( w[0], nx[2] );
		}
		
		nz[0][z] = CubicBSplineInterpolateWeights( w[1], ny[0] );
		nz[1][z] = CubicBSplineInterpolateWeights( w[1], ny[1] );
		nz[2][z] = CubicBSplineInterpolateWeights( w[1], ny[2] );
	}		
	N[0] = CubicBSplineInterpolateWeights( w[2], nz[0] );
	N[1] = CubicBSplineInterpolateWeights( w[2], nz[1] );
	N[2] = CubicBSplineInterpolateWeights( w[2], nz[2] );
}

elem_t CCCData::GetTriCubicBSplineFast( const vector3_t p ) const
{
	int floorX = (int)p[0];
	int floorY = (int)p[1];
	int floorZ = (int)p[2];
	
	int xNode[4] = {
		floorX - 1,
		floorX,
		floorX + 1,
		floorX + 2,
	};

	int yNode[4] = {
		floorY - 1,
		floorY,
		floorY + 1,
		floorY + 2,
	};

	int zNode[4] = {
		floorZ - 1,
		floorZ,
		floorZ + 1,
		floorZ + 2,
	};
		
	elem_t dx = p[0] - floorX;
	elem_t dy = p[1] - floorY;
	elem_t dz = p[2] - floorZ;
	
	elem_t nx[4],ny[4],nz[4];
	double t[3];
	
	for( int z=0;z<4;z++){
		for(int y=0;y<4;y++){
			
			for( int x=0;x<4;x++){
				nx[x] = Fetch( xNode[x], yNode[y], zNode[z] );
			}
			t[0] = dx;
			t[1] = dx*dx;
			t[2] = t[1]*dx;
			ny[y] = CubicBSplineInterpolate( t, nx[0],nx[1],nx[2],nx[3] );
		}
		t[0] = dy;
		t[1] = dy*dy;
		t[2] = t[1]*dy;
		nz[z] = CubicBSplineInterpolate( t, ny[0],ny[1],ny[2],ny[3] );
	}
	
	t[0] = dz;
	t[1] = dz*dz;
	t[2] = t[1]*dz;
	return CubicBSplineInterpolate( t, nz[0],nz[1],nz[2],nz[3] );
}

void CCCData::GetNormalTricubicBSplineFast( const vector3_t p, vector3_t N) const
{
	int floorX = (int)p[0];
	int floorY = (int)p[1];
	int floorZ = (int)p[2];
	
	int xNode[4] = {
		floorX - 1,
		floorX,
		floorX + 1,
		floorX + 2,
	};

	int yNode[4] = {
		floorY - 1,
		floorY,
		floorY + 1,
		floorY + 2,
	};

	int zNode[4] = {
		floorZ - 1,
		floorZ,
		floorZ + 1,
		floorZ + 2,
	};
		
	elem_t dx = p[0] - floorX;
	elem_t dy = p[1] - floorY;
	elem_t dz = p[2] - floorZ;
	double t[3];
	//for(int ordinate=0;ordinate<3;ordinate++){
		elem_t nx[3][4],ny[3][4],nz[3][4];
		
		for( int z=0;z<4;z++){
			for(int y=0;y<4;y++){
				
				for( int x=0;x<4;x++){
					int index[3] = {xNode[x], yNode[y], zNode[z]};		
										
					nx[0][x] = GetValueNorm( index, 0 );
					nx[1][x] = GetValueNorm( index, 1 );
					nx[2][x] = GetValueNorm( index, 2 );
				}

				t[0] = dx;
				t[1] = dx*dx;	
				t[2] = t[1]*dx;

				ny[0][y] = CubicBSplineInterpolate( t, nx[0][0],nx[0][1],nx[0][2],nx[0][3] );
				ny[1][y] = CubicBSplineInterpolate( t, nx[1][0],nx[1][1],nx[1][2],nx[1][3] );
				ny[2][y] = CubicBSplineInterpolate( t, nx[2][0],nx[2][1],nx[2][2],nx[2][3] );
			}
			
			t[0] = dy;
			t[1] = dy*dy;	
			t[2] = t[1]*dy;
			
			nz[0][z] = CubicBSplineInterpolate( t, ny[0][0],ny[0][1],ny[0][2],ny[0][3] );
			nz[1][z] = CubicBSplineInterpolate( t, ny[1][0],ny[1][1],ny[1][2],ny[1][3] );
			nz[2][z] = CubicBSplineInterpolate( t, ny[2][0],ny[2][1],ny[2][2],ny[2][3] );
		}
		
		t[0] = dz;
		t[1] = dz*dz;	
		t[2] = t[1]*dz;
			
		N[0] = CubicBSplineInterpolate( t, nz[0][0],nz[0][1],nz[0][2],nz[0][3] );
		N[1] = CubicBSplineInterpolate( t, nz[1][0],nz[1][1],nz[1][2],nz[1][3] );
		N[2] = CubicBSplineInterpolate( t, nz[2][0],nz[2][1],nz[2][2],nz[2][3] );
	//}
}


// void CCCData::GetNormalTricubicBSplineFast( const vector3_t p, vector3_t N) const
// {
// 	int floorX = (int)p[0];
// 	int floorY = (int)p[1];
// 	int floorZ = (int)p[2];
// 	
// 	int xNode[4] = {
// 		floorX - 1,
// 		floorX,
// 		floorX + 1,
// 		floorX + 2,
// 	};
// 
// 	int yNode[4] = {
// 		floorY - 1,
// 		floorY,
// 		floorY + 1,
// 		floorY + 2,
// 	};
// 
// 	int zNode[4] = {
// 		floorZ - 1,
// 		floorZ,
// 		floorZ + 1,
// 		floorZ + 2,
// 	};
// 		
// 	elem_t dx = p[0] - floorX;
// 	elem_t dy = p[1] - floorY;
// 	elem_t dz = p[2] - floorZ;
// 	double t[3];
// 	
// 	for(int ordinate=0;ordinate<3;ordinate++){
// 		
// 		elem_t nx[4],ny[4],nz[4];
// 		
// 		for( int z=0;z<4;z++){
// 			for(int y=0;y<4;y++){
// 				
// 				for( int x=0;x<4;x++){
// 					int index[3] = {xNode[x], yNode[y], zNode[z]};												
// 					nx[x] = GetValueNorm( index, 0 );					
// 				}
// 
// 				t[0] = dx;
// 				t[1] = dx*dx;	
// 				t[2] = t[1]*dx;
// 
// 				ny[y] = CubicBSplineInterpolate( t, nx[0],nx[1],nx[2],nx[3] );
// 				//ny[1][y] = CubicBSplineInterpolate( t, nx[1][0],nx[1][1],nx[1][2],nx[1][3] );
// 				//ny[2][y] = CubicBSplineInterpolate( t, nx[2][0],nx[2][1],nx[2][2],nx[2][3] );
// 			}
// 			
// 			t[0] = dy;
// 			t[1] = dy*dy;	
// 			t[2] = t[1]*dy;
// 			
// 			nz[z] = CubicBSplineInterpolate( t, ny[0],ny[1],ny[2],ny[3] );
// 			//nz[1][z] = CubicBSplineInterpolate( t, ny[1][0],ny[1][1],ny[1][2],ny[1][3] );
// 			//nz[2][z] = CubicBSplineInterpolate( t, ny[2][0],ny[2][1],ny[2][2],ny[2][3] );
// 		}
// 		
// 		t[0] = dz;
// 		t[1] = dz*dz;	
// 		t[2] = t[1]*dz;
// 			
// 		N[ordinate] = CubicBSplineInterpolate( t, nz[0],nz[1],nz[2],nz[3] );
// 		//N[1] = CubicBSplineInterpolate( t, nz[1][0],nz[1][1],nz[1][2],nz[1][3] );
// 		//N[2] = CubicBSplineInterpolate( t, nz[2][0],nz[2][1],nz[2][2],nz[2][3] );
// 	}
// }


elem_t CCCData::GetTriCubicBSplineNew( const vector3_t p ) const
{
	int floorX = (int)p[0];
	int floorY = (int)p[1];
	int floorZ = (int)p[2];

	elem_t minusX = -p[0];
	elem_t minusY = -p[1];
	elem_t minusZ = -p[2];

	int xNode[4] = {
		floorX - 1,
		floorX,
		floorX + 1,
		floorX + 2,
	};

	int yNode[4] = {
		floorY - 1,
		floorY,
		floorY + 1,
		floorY + 2,
	};

	int zNode[4] = {
		floorZ - 1,
		floorZ,
		floorZ + 1,
		floorZ + 2,
	};

	elem_t dx[4] = {
		xNode[0] + minusX,
		xNode[1] + minusX,
		xNode[2] + minusX,
		xNode[3] + minusX,
	};

	elem_t dy[4] = {
		yNode[0] + minusY,
		yNode[1] + minusY,
		yNode[2] + minusY,
		yNode[3] + minusY,
	};

// 	elem_t dz[4] = {
// 		zNode[0] + minusZ,
// 		zNode[1] + minusZ,
// 		zNode[2] + minusZ,
// 		zNode[3] + minusZ,
// 	};

	//printf("New\n");
	elem_t sumX;
	elem_t sumY;
	elem_t sumResult = 0;

	for(int k=0;k<4;++k){
		sumY = 0;
		for(int j=0;j<4;++j){
			sumX = 0;
			for(int i=0;i<4;++i){
				sumX += Fetch( xNode[i],yNode[j],zNode[k] ) * CubicBSpline( dx[i] );
			}
			sumY += sumX * CubicBSpline( dy[j] );
		}
		sumResult += sumY * CubicBSpline( zNode[k] + minusZ );
	}
	return sumResult;
}

double CCCData::GetTriCubicBSpline( const vector3_t p ) const

{
	int   GridLines[3][4];
	double Offset[3][4];


	int floorX = (int)p[0];
	int floorY = (int)p[1];
	int floorZ = (int)p[2];

	GridLines[X][0] = floorX - 1;
	GridLines[X][1] = floorX;
	GridLines[X][2] = floorX + 1;
	GridLines[X][3] = floorX + 2;

	GridLines[Y][0] = floorY - 1;
	GridLines[Y][1] = floorY;
	GridLines[Y][2] = floorY + 1;
	GridLines[Y][3] = floorY + 2;

	GridLines[Z][0] = floorZ - 1;
	GridLines[Z][1] = floorZ;
	GridLines[Z][2] = floorZ + 1;
	GridLines[Z][3] = floorZ + 2;

	Offset[X][0] = GridLines[X][0] - p[0];
	Offset[X][1] = GridLines[X][1] - p[0];
	Offset[X][2] = GridLines[X][2] - p[0];
	Offset[X][3] = GridLines[X][3] - p[0];

	Offset[Y][0] = GridLines[Y][0] - p[1];
	Offset[Y][1] = GridLines[Y][1] - p[1];
	Offset[Y][2] = GridLines[Y][2] - p[1];
	Offset[Y][3] = GridLines[Y][3] - p[1];

	Offset[Z][0] = GridLines[Z][0] - p[2];
	Offset[Z][1] = GridLines[Z][1] - p[2];
	Offset[Z][2] = GridLines[Z][2] - p[2];
	Offset[Z][3] = GridLines[Z][3] - p[2];

	// Now be sane with the Gridlines; Messy part :( yuck

	for(int lines=0;lines < 4;++lines){

		if( GridLines[X][lines] < 0 ) GridLines[X][lines] = 0;
		else if( GridLines[X][lines] >= m_xSize ) GridLines[X][lines] = m_xSize - 1;

		if( GridLines[Y][lines] < 0 ) GridLines[Y][lines] = 0;
		else if( GridLines[Y][lines] >= m_ySize ) GridLines[Y][lines] = m_ySize - 1;

		if( GridLines[Z][lines] < 0 ) GridLines[Z][lines] = 0;
		else if( GridLines[Z][lines] >= m_zSize ) GridLines[Z][lines] = m_zSize - 1;

	}

	double sum =0;
	double BSplineWeightAlongY, BSplineWeightAlongZ;

	for(int k=0;k<4;++k){

		BSplineWeightAlongZ = CubicBSpline( Offset[Z][k] );

		for(int j=0;j<4;++j){

			BSplineWeightAlongY = CubicBSpline( Offset[Y][j] );

			for(int i=0;i<4;++i){

				sum += Fetch( GridLines[X][i], GridLines[Y][j],GridLines[Z][k] ) * CubicBSpline( Offset[X][i] ) * BSplineWeightAlongY * BSplineWeightAlongZ;

			}
		}
	}

	return sum;
}

void CCCData::GetNormalNearest( const vector3_t p, vector3_t N) const
{
	glmCopyVector3d(m_pppNormals[(int)(p[2]+0.5f)+GUARD_BAND][(int)(p[1]+0.5f)+GUARD_BAND][(int)(p[0]+0.5f)+GUARD_BAND],N);
}

void glmLerpVector3d( vector3_t output, const vector3_t a , const double ak, const vector3_t b , const double bk )
{
	output[X] = a[X] * ak + b[X] * bk;
	output[Y] = a[Y] * ak + b[Y] * bk;
	output[Z] = a[Z] * ak + b[Z] * bk;
}


void CCCData::GetNormalTriLinear( const vector3_t v, vector3_t N) const
{
	vector3_t p[8];

	memset(p,0,sizeof(vector3_t) * 8 );

	int floorX = ((int)v[0]);
	int nextX  = floorX + 1;
	int floorY = ((int)v[1]);
	int nextY  = floorY+1;
	int floorZ = ((int)v[2]);
	int nextZ  = floorZ+1;


	bool nextXWithinBound = nextX < m_xSize;
	bool nextYWithinBound = nextY < m_ySize;
	bool nextZWithinBound = nextZ < m_zSize;

	if( nextZWithinBound ){
		FetchNormal( floorX, floorY, nextZ, p[0] );

		if( nextXWithinBound ) {

			FetchNormal( nextX, floorY, nextZ, p[1] );

			if( nextYWithinBound ){
				FetchNormal(nextX, nextY,nextZ, p[2]);
			}
		}

		if( nextYWithinBound ){
			FetchNormal( floorX, nextY, nextZ, p[3] );
		}

	}

	FetchNormal( floorX, floorY, floorZ, p[4] );

	if( nextXWithinBound ){
		FetchNormal( nextX, floorY, floorZ, p[5] );

		if( nextYWithinBound ){
			FetchNormal( nextX, nextY, floorZ, p[6] );
		}
	}

	if( nextYWithinBound ){
		FetchNormal(floorX, nextY, floorZ, p[7] );
	}

	double dx = v[0] - floorX;
	double dy = v[1] - floorY;
	double dz = v[2] - floorZ;

	double invDx = 1 - dx;
	double invDy = 1 - dy;
	double invDz = 1 - dz;

	vector3_t interp[6];

	// Hell lot of linear interpolations going on
	glmLerpVector3d( interp[0], p[0], invDx, p[1], dx );
	glmLerpVector3d( interp[1], p[4], invDx, p[5], dx );
	glmLerpVector3d( interp[2], interp[0],dz, interp[1], invDz );
	glmLerpVector3d( interp[3], p[3], invDx, p[2], dx );
	glmLerpVector3d( interp[4], p[7], invDx, p[6], dx );
	glmLerpVector3d( interp[5], interp[3], dz, interp[4], invDz );
	glmLerpVector3d( N, interp[2], invDy, interp[5], dy );
}


void CCCData::GetNormalTriCubicBSplineNew( const vector3_t p, vector3_t N) const
{
	// This is obviously not the most optimal way to compute the gradient
	// but we wanted to be fair with how BCC interpolation was implemented
	// where all the indexings and everything would be recalculated.
	int floorX = (int)p[0];
		int floorY = (int)p[1];
		int floorZ = (int)p[2];

		elem_t minusX = -p[0];
		elem_t minusY = -p[1];
		elem_t minusZ = -p[2];

		int xNode[4] = {
			floorX - 1,
			floorX,
			floorX + 1,
			floorX + 2,
		};

		int yNode[4] = {
			floorY - 1,
			floorY,
			floorY + 1,
			floorY + 2,
		};

		int zNode[4] = {
			floorZ - 1,
			floorZ,
			floorZ + 1,
			floorZ + 2,
		};

		elem_t dx[4] = {
			xNode[0] + minusX,
			xNode[1] + minusX,
			xNode[2] + minusX,
			xNode[3] + minusX,
		};

		elem_t dy[4] = {
			yNode[0] + minusY,
			yNode[1] + minusY,
			yNode[2] + minusY,
			yNode[3] + minusY,
		};
		
	for(int ordinate=0;ordinate<3;++ordinate){


	// 	elem_t dz[4] = {
	// 		zNode[0] + minusZ,
	// 		zNode[1] + minusZ,
	// 		zNode[2] + minusZ,
	// 		zNode[3] + minusZ,
	// 	};

		//printf("New\n");
		elem_t sumX;
		elem_t sumY;
		elem_t sumResult = 0;

		for(int k=0;k<4;++k){
			sumY = 0;
			for(int j=0;j<4;++j){
				sumX = 0;
				for(int i=0;i<4;++i){
					int index[3] = {xNode[i],yNode[j],zNode[k]};
					sumX += GetValueNorm( index, ordinate ) * CubicBSpline( dx[i] );
				}
				sumY += sumX * CubicBSpline( dy[j] );
			}
			sumResult += sumY * CubicBSpline( zNode[k] + minusZ );
		}
		N[ordinate] =  sumResult;
	}
	//printf("%f %f %f\n",N[0],N[1],N[2]);
}

void CCCData::GetNormalTriCubicBSpline( const vector3_t p, vector3_t N) const
{
	int   GridLines[3][4];
	double Offset[3][4];


	int floorX = (int)p[0];
	int floorY = (int)p[1];
	int floorZ = (int)p[2];

	GridLines[X][0] = floorX - 1;
	GridLines[X][1] = floorX;
	GridLines[X][2] = floorX + 1;
	GridLines[X][3] = floorX + 2;

	GridLines[Y][0] = floorY - 1;
	GridLines[Y][1] = floorY;
	GridLines[Y][2] = floorY + 1;
	GridLines[Y][3] = floorY + 2;

	GridLines[Z][0] = floorZ - 1;
	GridLines[Z][1] = floorZ;
	GridLines[Z][2] = floorZ + 1;
	GridLines[Z][3] = floorZ + 2;

	Offset[X][0] = GridLines[X][0] - p[0];
	Offset[X][1] = GridLines[X][1] - p[0];
	Offset[X][2] = GridLines[X][2] - p[0];
	Offset[X][3] = GridLines[X][3] - p[0];

	Offset[Y][0] = GridLines[Y][0] - p[1];
	Offset[Y][1] = GridLines[Y][1] - p[1];
	Offset[Y][2] = GridLines[Y][2] - p[1];
	Offset[Y][3] = GridLines[Y][3] - p[1];

	Offset[Z][0] = GridLines[Z][0] - p[2];
	Offset[Z][1] = GridLines[Z][1] - p[2];
	Offset[Z][2] = GridLines[Z][2] - p[2];
	Offset[Z][3] = GridLines[Z][3] - p[2];

	// Now be sane with the Gridlines; Messy part :( yuck

	for(int lines=0;lines < 4;++lines){

		if( GridLines[X][lines] < 0 ) GridLines[X][lines] = 0;
		else if( GridLines[X][lines] >= m_xSize ) GridLines[X][lines] = m_xSize - 1;

		if( GridLines[Y][lines] < 0 ) GridLines[Y][lines] = 0;
		else if( GridLines[Y][lines] >= m_ySize ) GridLines[Y][lines] = m_ySize - 1;

		if( GridLines[Z][lines] < 0 ) GridLines[Z][lines] = 0;
		else if( GridLines[Z][lines] >= m_zSize ) GridLines[Z][lines] = m_zSize - 1;

	}

	vector3_t sum ={0,0,0};
	double BSplineWeightAlongY, BSplineWeightAlongZ;
	vector3_t normal;

	for(int k=0;k<4;++k){

		BSplineWeightAlongZ = CubicBSpline( Offset[Z][k] );

		for(int j=0;j<4;++j){

			BSplineWeightAlongY = CubicBSpline( Offset[Y][j] );

			for(int i=0;i<4;++i){

				FetchNormal( GridLines[X][i], GridLines[Y][j],GridLines[Z][k], normal );
				glmMultVectorByScalar3d(normal,CubicBSpline( Offset[X][i] ) * BSplineWeightAlongY * BSplineWeightAlongZ);
				glmAddVector3d( sum, sum, normal );

			}
		}
	}

	glmCopyVector3d( sum, N );
}

void CCCData::HookGradientComponentEstimators( )
{
	if( m_bEpsilonNormalMode ) return;
	
	if( !m_bOnTheFlyNormals ){
		printf("Setting up normal estimators for precomputed normals\n");
		m_gradEstimator[0] = new CCCGradientXEstimatorPrecomp( this );
		m_gradEstimator[1] = new CCCGradientYEstimatorPrecomp( this );
		m_gradEstimator[2] = new CCCGradientZEstimatorPrecomp( this );
	}else{
		printf("Setting up normal estimators for on-the-fly normals\n");
		switch( m_normalEst){
			case CC_NORMAL_CD_2ND:
				m_gradEstimator[0] = new CCESTIMATORCLASSNAME(X,CentralDifference2nd)( this );
				m_gradEstimator[1] = new CCESTIMATORCLASSNAME(Y,CentralDifference2nd)( this );
				m_gradEstimator[2] = new CCESTIMATORCLASSNAME(Z,CentralDifference2nd)( this );
				break;
			case CC_NORMAL_CD_4TH:
				m_gradEstimator[0] = new CCESTIMATORCLASSNAME(X,CentralDifference4th)( this );
				m_gradEstimator[1] = new CCESTIMATORCLASSNAME(Y,CentralDifference4th)( this );
				m_gradEstimator[2] = new CCESTIMATORCLASSNAME(Z,CentralDifference4th)( this );
				break;
			default:
				printf("The normal estimation filter is not supported for On-the-fly computation!\n");
				break;

		};

	}
}

