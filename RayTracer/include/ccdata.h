#include "data.h"

#include <map>
#include <string>

class CCCData : public CBaseData{
public:
	enum RECONSTRUCTION_FILTER{
		CC_NEAREST,
  		CC_TRILINEAR,
    		CC_TRICUBICBSPLINE
	};

	enum NORMAL_DISCRETE_FILTER{
		CC_NORMAL_CD_2ND,
		CC_NORMAL_CD_4TH,
		CC_NORMAL_CD_6TH,
		CC_NORMAL_CD_10TH,
	};

	RECONSTRUCTION_FILTER   m_reconFilter;
	RECONSTRUCTION_FILTER   m_normalReconFilter;

	NORMAL_DISCRETE_FILTER  m_normalEst;
	elem_t			m_invScale[3];

	std::string             m_gridPointDumpFile;

	static std::map< std::string, RECONSTRUCTION_FILTER >  m_reconstructionFilterNameMap;
	static std::map< std::string, NORMAL_DISCRETE_FILTER > m_normalFilterNameMap;
	static bool m_bMapInitialized;
	static void InitializeMaps();

	elem_t ComputeNormalXCentralDifference2nd( int x, int y, int z ) const;
	elem_t ComputeNormalYCentralDifference2nd( int x, int y, int z ) const;
	elem_t ComputeNormalZCentralDifference2nd( int x, int y, int z ) const;
	void   PrecompNormalCentralDifference2nd( vector3_t N, int x,int y, int z );

	elem_t ComputeNormalXCentralDifference4th( int x, int y, int z ) const;
	elem_t ComputeNormalYCentralDifference4th( int x, int y, int z ) const;
	elem_t ComputeNormalZCentralDifference4th( int x, int y, int z ) const;
	void   PrecompNormalCentralDifference4th( vector3_t N, int x,int y, int z );


	elem_t ComputeNormalXCentralDifference6th( int x, int y, int z ) const;
	elem_t ComputeNormalYCentralDifference6th( int x, int y, int z ) const;
	elem_t ComputeNormalZCentralDifference6th( int x, int y, int z ) const;
	void   PrecompNormalCentralDifference6th( vector3_t N, int x,int y, int z );

	elem_t ComputeNormalXCentralDifference10th( int x, int y, int z ) const;
	elem_t ComputeNormalYCentralDifference10th( int x, int y, int z ) const;
	elem_t ComputeNormalZCentralDifference10th( int x, int y, int z ) const;
	void   PrecompNormalCentralDifference10th( vector3_t N, int x,int y, int z );
public:
	CCCData( const char* name ) : CBaseData( std::string( name ) )
	{
		m_reconFilter = CC_TRICUBICBSPLINE;
		m_normalEst = CC_NORMAL_CD_4TH;
	}

	CCCData( const char* name, const char* fname, double xScale, double yScale, double zScale );
	CCCData( const char* name, int xSize, int ySize, int zSize, double xScale, double yScale, double zScale, atomic_t *pData, bool bOnTheFlyNormals = false ) : CBaseData( std::string(name) )
	{
		m_normalEst = CC_NORMAL_CD_4TH;
		m_reconFilter = m_normalReconFilter = CC_TRICUBICBSPLINE;
		CreateData( xSize, ySize, zSize, xScale, yScale, zScale, pData );

		m_invScale[X] = 1.0/GetVolume()->m_xScale;
		m_invScale[Y] = 1.0/GetVolume()->m_yScale;
		m_invScale[Z] = 1.0/GetVolume()->m_zScale;

		m_bOnTheFlyNormals = bOnTheFlyNormals;
	}

	CCCData( const char *name, const char *normalEstName, bool bOnTheFlyNormals = false, const char* gridPointDumpFile=NULL);
	CCCData( const char* name, const char* fname, DATA_FORMAT fmt, double xScale, double yScale, double zScale, const char* normalEstName, bool bOnTheFlyNormals = false, const char *gridPointDumpFile = NULL);
	CCCData( const char* name, const char* fname, DATA_FORMAT fmt, const char* xname, const char *yname, const char *zname, double xScale, double yScale, double zScale, const char *gridPointDumpFile );

	virtual ~CCCData()
	{
		if( m_sepNormals[0]  /*if one is full, others are full too*/){
			for(int z=0;z<(m_zSize+2*GUARD_BAND);z++){
				for(int y=0;y<(m_ySize+2*GUARD_BAND);++y){
					DELETE_SAFE_ARRAY(m_sepNormals[X][z][y]);
					DELETE_SAFE_ARRAY(m_sepNormals[Y][z][y]);
					DELETE_SAFE_ARRAY(m_sepNormals[Z][z][y]);
				}
				DELETE_SAFE_ARRAY(m_sepNormals[X][z]);
				DELETE_SAFE_ARRAY(m_sepNormals[Y][z]);
				DELETE_SAFE_ARRAY(m_sepNormals[Z][z]);
			}
			DELETE_SAFE_ARRAY(m_sepNormals[X]);
			DELETE_SAFE_ARRAY(m_sepNormals[Y]);
			DELETE_SAFE_ARRAY(m_sepNormals[Z]);
		}

		DELETE_SAFE( m_gradEstimator[0] );
		DELETE_SAFE( m_gradEstimator[1] );
		DELETE_SAFE( m_gradEstimator[2] );
	}

	virtual void Evaluate( const vector3_t min, const vector3_t max, const vector3_t step,
		elem_t (*fn)(const vector3_t p), void (*fnNormal)(const vector3_t p, vector3_t N) = 0 );

	virtual void  CalculateNormals( void (*fnNormal)(const vector3_t p, vector3_t N) = 0 );

	elem_t GetValueNorm( int index[], int ordinate ) const
	{
	/*	if( index[X] < 0 || index[X] >= m_xSize ||
		  index[Y] < 0 || index[Y] >= m_ySize ||
		  index[Z] < 0 || index[Z] >= m_zSize ) return 0;

	  	return m_sepNormals[ordinate][index[Z]][index[Y]][index[X]];*/

		return m_gradEstimator[ordinate]->GetGradientComponent( index[0], index[1], index[2] );
	}

	double GetNearest( const vector3_t p ) const ;
	double GetTriLinear( const vector3_t p ) const ;
	double GetTriCubicBSpline( const vector3_t p ) const;
	elem_t GetTriCubicBSplineNew( const vector3_t p ) const;
	elem_t GetTriCubicBSplineFast( const vector3_t p) const;
	elem_t GetTriCubicBSplineTurbo( const vector3_t p) const;
	
	void GetNormalNearest( const vector3_t p , vector3_t N) const ;
	void GetNormalTriLinear( const vector3_t p, vector3_t N) const ;
	void GetNormalTriCubicBSpline( const vector3_t p, vector3_t N) const;

	void GetNormalTriCubicBSplineNew( const vector3_t p, vector3_t N) const;
	void GetNormalTricubicBSplineFast( const vector3_t p, vector3_t N) const;
	void GetNormalTriCubicBSplineTurbo( const vector3_t p, vector3_t N) const;

	double GetFromVolumeSpace( const vector3_t p ) const
	{
		switch (m_reconFilter){
			case CC_NEAREST:
				return GetNearest( p );
			case CC_TRILINEAR:
				return GetTriLinear( p );
			case CC_TRICUBICBSPLINE:
#ifdef TRICUBIC_BSPLINE_TURBO
			    return GetTriCubicBSplineTurbo(p);
#else				
				return GetTriCubicBSplineFast( p );
#endif				
		}
		return 0;
	}

	void GetNormalFromVolumeSpace( const vector3_t p, vector3_t N ) const
	{
		switch (m_normalReconFilter){
			case CC_NEAREST:
				GetNormalNearest( p, N );
				 break;
			case CC_TRILINEAR:
				GetNormalTriLinear(p, N );
				break;
			case CC_TRICUBICBSPLINE:
#ifdef TRICUBIC_BSPLINE_TURBO
				GetNormalTriCubicBSplineTurbo(p,N);
#else				
				GetNormalTricubicBSplineFast( p, N );
#endif				
				break;
		}
	}

	virtual double Get( const vector3_t p ) const
	{
		vector3_t xyz;
		m_pVolume->FindXYZ( p,xyz );
		//return GetFromVolumeSpace( xyz[X], xyz[Y], xyz[Z] );
		return GetFromVolumeSpace( xyz );
	}

	virtual void GetNormal( const vector3_t p, vector3_t N ) const
	{
		if( m_bEpsilonNormalMode ){
			CBaseData::GetNormal( p, N );
			return;
		}

		vector3_t xyz;
		m_pVolume->FindXYZ( p,xyz );
		//GetNormalFromVolumeSpace( xyz[X], xyz[Y],xyz[Z], N );
		GetNormalFromVolumeSpace( xyz, N );

	}

	void SetReconstructionFilter( RECONSTRUCTION_FILTER mode )
	{
		m_reconFilter = mode;
	}

	RECONSTRUCTION_FILTER GetReconstructionFilter() { return m_reconFilter; }

	virtual bool SetNormalEstimatorFilterByName( const char* name );
	virtual bool SetReconstructionFilterByName( const char* name );
	virtual bool SetNormalReconstructionFilterByName( const char *filter );
	virtual void HookGradientComponentEstimators( );
};

class CCCGradientComponentEstimator :public CBaseGradientComponentEstimator{
	protected: CCCData *m_pData;
	public: CCCGradientComponentEstimator( CCCData *pData )
	: m_pData( pData ){}
};

class CCCGradientXEstimatorPrecomp : public CCCGradientComponentEstimator{
	public:
		CCCGradientXEstimatorPrecomp(CCCData *pData )
	: CCCGradientComponentEstimator( pData)
		{
			printf("\tCCCGradientXEstimatorPrecomp Initialized\n");
		}
		virtual elem_t GetGradientComponent( int x, int y, int z ){
//			if( x < 0 || x >= m_pData->m_xSize ||
//						 y < 0 || y >= m_pData->m_ySize ||
//						 z < 0 || z >= m_pData->m_zSize ) return 0;
			//return m_pData->m_sepNormals[0][z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND];
			return m_pData->m_pppNormals[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND][0];
		}
};

class CCCGradientYEstimatorPrecomp : public CCCGradientComponentEstimator{
	public:
		CCCGradientYEstimatorPrecomp(CCCData *pData )
	: CCCGradientComponentEstimator( pData)
		{
			printf("\tCCCGradientYEstimatorPrecomp Initialized\n");
		}
		virtual elem_t GetGradientComponent( int x, int y, int z ){
//			if( x < 0 || x >= m_pData->m_xSize ||
//						 y < 0 || y >= m_pData->m_ySize ||
//						 z < 0 || z >= m_pData->m_zSize ) return 0;
			//return m_pData->m_sepNormals[1][z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND];
			return m_pData->m_pppNormals[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND][1];
		}
};

class CCCGradientZEstimatorPrecomp : public CCCGradientComponentEstimator{
	public:
		CCCGradientZEstimatorPrecomp(CCCData *pData )
	: CCCGradientComponentEstimator( pData)
		{
			printf("\tCCCGradientZEstimatorPrecomp Initialized\n");
		}
		virtual elem_t GetGradientComponent( int x, int y, int z ){
//			if( x < 0 || x >= m_pData->m_xSize ||
//						 y < 0 || y >= m_pData->m_ySize ||
//						 z < 0 || z >= m_pData->m_zSize ) return 0;
			//return m_pData->m_sepNormals[2][z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND];
			return m_pData->m_pppNormals[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND][2];
		}
};

#define CCESTIMATORCLASSNAME(component,name) CCCGradient##component##Estimator##name
#define CCESTIMATOR(component,name) \
class CCESTIMATORCLASSNAME(component,name) : public CCCGradientComponentEstimator{ 	\
	public:	 								     	\
		CCESTIMATORCLASSNAME(component,name) (CCCData *pData ) 		\
	: CCCGradientComponentEstimator( pData)					\
		{									\
			printf("\t" STR(CCESTIMATORCLASSNAME(component,name)) " Initialized\n");\
		}									\
		virtual elem_t GetGradientComponent( int x, int y, int z )		\
		{									\
			return m_pData->ComputeNormal##component##name( x,y,z);	\
		}									\
};

// 2ND
CCESTIMATOR(X,CentralDifference2nd)
CCESTIMATOR(Y,CentralDifference2nd)
CCESTIMATOR(Z,CentralDifference2nd)

//4TH
CCESTIMATOR(X,CentralDifference4th)
CCESTIMATOR(Y,CentralDifference4th)
CCESTIMATOR(Z,CentralDifference4th)

