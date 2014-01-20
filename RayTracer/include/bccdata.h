#ifndef __BCC_DATA_H__
#define __BCC_DATA_H__

#include "data.h"
#include <map>


class CBCCData : public CBaseData{
public:
	typedef enum{
		BCC_NORMAL_BCD,
		BCC_NORMAL_SOCD,
  		BCC_NORMAL_ICD,
		BCC_NORMAL_SIZE_OPTIMAL16,
		BCC_NORMAL_SEMI_OPTIMAL18,
		BCC_NORMAL_OPTIMAL26
	} NORMAL_DISCRETE_FILTER;

	typedef enum{
		BCC_QUINTIC_BOX_SPLINE,
		BCC_LINEAR_BOX_SPLINE,
	} RECONSTRUCTION_FILTER;

public:

	static bool m_bMapInitialized;
	static std::map< std::string, NORMAL_DISCRETE_FILTER> m_normalFilterNameMap;
	static std::map< std::string, RECONSTRUCTION_FILTER>  m_reconstructionFilterNameMap;

	static void InitializeMaps();


	//elem_t ***m_sepNormals[3];

	NORMAL_DISCRETE_FILTER m_normalEst;
	RECONSTRUCTION_FILTER  m_reconFilter;
	RECONSTRUCTION_FILTER  m_normalReconFilter;

	char       m_gridNormalOutputFile[256];


	
	inline void BCCToIndex( int index[3], int x, int y, int z ) const
	{

		if ( z & 0x1 ) { // checking if z is odd
		  x--;
		  y--;
		}

		x >>= 1;
		y >>= 1;

		index[0] = x;
		index[1] = y;
		index[2] = z;

	}

	inline void IndexToBCC( int &x, int &y, int &z, int index[3] )
	{
		x = index[0] << 1;
		y = index[1] << 1;
		z = index[2];

		if( z & 0x1 ) {
			++x;
			++y;
		}
	}


	double GetValue( int index[]  ) const
	{

	  /*if( index[X] < 0 || index[X] >= m_xSize ||
		  index[Y] < 0 || index[Y] >= m_ySize ||
		  index[Z] < 0 || index[Z] >= m_zSize ) return 0;
	  */

	  return Fetch( index[X], index[Y], index[Z] );
	}

	elem_t GetValueNorm( int index[], int XYZ ) const
	{
// 	  if( index[X] < 0 || index[X] >= m_xSize ||
// 		  index[Y] < 0 || index[Y] >= m_ySize ||
// 		  index[Z] < 0 || index[Z] >= m_zSize ) return 0;
//
// 	  return m_sepNormals[XYZ][index[Z]][index[Y]][index[X]];

		/*switch( XYZ ){
			case 0:
				return ComputeNormalXCentralDifference(index[0],index[1],index[2]);
				break;
			case 1:
				return ComputeNormalYCentralDifference(index[0],index[1],index[2]);
				break;
			case 2:
				return ComputeNormalZCentralDifference(index[0],index[1],index[2]);
				break;

		}
		return 0;*/

		return m_gradEstimator[XYZ]->GetGradientComponent( index[0], index[1], index[2] );
	}

	double rho22(double alpha, double beta, double gamma) const
	{
		register double alpha2 = alpha * alpha;
		return
		  -alpha2*alpha*
		  (one_over_6*gamma*beta-two_over_24*alpha*
		  (gamma+beta)+six_over_120 *alpha2);
	}

	void Permute( elem_t o[8], const elem_t i[8], const int p[8]) const
	{
		for(int k=1;k<7;++k)
		{
			o[k] = i[ p[k] ];
		}
	}
				  
	inline void WorldToBCC( vector3_t bcc, const vector3_t world ) const
	{
	  //bcc[X] = (world[X] + m_pVolume->m_HXSize) / m_halfScale[X];
	  //bcc[Y] = (world[Y] + m_pVolume->m_HYSize) / m_halfScale[Y];
	  //bcc[Z] = (world[Z] + m_pVolume->m_HZSize) / m_halfScale[Z];
	  
	  bcc[X] = world[X] * m_invHalfScale[X] + m_translate[X];
	  bcc[Y] = world[Y] * m_invHalfScale[Y] + m_translate[Y];
	  bcc[Z] = world[Z] * m_invHalfScale[Z] + m_translate[Z];

	}

	inline void BCCToWorld( vector3_t world, const vector3_t bcc) const
	{
		world[X] = bcc[X] * m_halfScale[X] - m_pVolume->m_HXSize;
		world[Y] = bcc[Y] * m_halfScale[Y] - m_pVolume->m_HYSize;
		world[Z] = bcc[Z] * m_halfScale[Z] - m_pVolume->m_HZSize;
	}

	double Interpolate( const vector3_t bcc ) const;
	
	double InterpolateQuinticBoxSplineTurbo( const vector3_t bcc ) const;
	double InterpolateLinearBoxSplineTurbo( const vector3_t bcc ) const;
	double InterpolateQuinticBoxSplineSteroid( const vector3_t bcc) const;
	void   InterpolateNormalQuinticBoxSplineTurbo( const vector3_t bcc, vector3_t N ) const;
	void   InterpolateNormalLinearBoxSplineTurbo( const vector3_t bcc, vector3_t N ) const;
	
	

	//SOCD
	elem_t  ComputeNormalXCentralDifference( int x,int y, int z) const ;
	elem_t  ComputeNormalYCentralDifference( int x,int y, int z) const ;
	elem_t  ComputeNormalZCentralDifference( int x,int y, int z) const ;
	void   PrecompNormalCentralDifference( vector3_t N, int x,int y, int z );


	void   PrecompNormalInterpolatedCentralDifference( vector3_t N, int x,int y, int z );

	//BCD
	elem_t  ComputeNormalXBoxDifference( int x,int y, int z) const ;
	elem_t  ComputeNormalYBoxDifference( int x,int y, int z) const ;
	elem_t  ComputeNormalZBoxDifference( int x,int y, int z) const ;
	void   PrecompNormalBoxDifference( vector3_t N, int x, int y, int z );

	//OPT16
	elem_t  ComputeNormalXSizeOptimalFilter( int x,int y, int z) const ;
	elem_t  ComputeNormalYSizeOptimalFilter( int x,int y, int z) const ;
	elem_t  ComputeNormalZSizeOptimalFilter( int x,int y, int z) const ;
	void   PrecompNormalSizeOptimalFilter( vector3_t N, int x, int y, int z );

	void   PrecompNormalSemiOptimalFilter( vector3_t N, int x, int y, int z );

	//OPT26
	elem_t  ComputeNormalXOptimalFilter( int x,int y, int z) const ;
	elem_t  ComputeNormalYOptimalFilter( int x,int y, int z) const ;
	elem_t  ComputeNormalZOptimalFilter( int x,int y, int z) const ;
	void   PrecompNormalOptimalFilter( vector3_t N, int x, int y, int z );
	
	void   Fill_1_Even(int id[],double p1[],double p2[], double p3[], double p4[],double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const;
	void   Fill_2_Even(int id[],double p1[],double p2[], double p3[], double p4[],double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const ;
	void   Fill_3_Even(int id[],double p1[],double p2[], double p3[], double p4[],double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const;
	void   Fill_4_Even(int id[],double p1[],double p2[], double p3[], double p4[],double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const;
	void   Fill_5_Even(int id[],double p1[],double p2[], double p3[], double p4[],double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const;
	void   Fill_6_Even(int id[],double p1[],double p2[], double p3[], double p4[],double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const;
	void   Fill_1_Odd(int id[],double p1[],double p2[], double p3[], double p4[],double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const;
	void   Fill_2_Odd(int id[],double p1[],double p2[], double p3[], double p4[],double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const;
	void   Fill_3_Odd(int id[],double p1[],double p2[], double p3[], double p4[],double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const;
	void   Fill_4_Odd(int id[],double p1[],double p2[], double p3[], double p4[],double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const;
	void   Fill_5_Odd(int id[],double p1[],double p2[], double p3[], double p4[],double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const;
	void   Fill_6_Odd(int id[],double p1[],double p2[], double p3[], double p4[],double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const;
	
	typedef void (CBCCData::*FNFILL) ( int id[],double p1[],double p2[], double p3[], double p4[],double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const;
	
	static CBCCData::FNFILL m_fnFill[16];
public:


	//Spacing.This should be half of the scaling in all three axes, i.e. halfScale[X] = xScale / 2.0f and so on
	double m_halfScale[3];
	double m_invHalfScale[3];
	double m_translate[3];

	CBCCData(const char* name, NORMAL_DISCRETE_FILTER normalEst = BCC_NORMAL_BCD, bool bOnTheFlyNormals = false, const char* gridNormalOutputFile = NULL ) : CBaseData(std::string(name)), m_normalEst(normalEst)
	{
		m_reconFilter = BCC_QUINTIC_BOX_SPLINE;
		m_normalReconFilter = BCC_QUINTIC_BOX_SPLINE;
		m_bOnTheFlyNormals = bOnTheFlyNormals;

		if( gridNormalOutputFile){
			strcpy( m_gridNormalOutputFile, gridNormalOutputFile );
		}
		else
			m_gridNormalOutputFile[0] = 0;

	}

	CBCCData(const char* name, const char* filter, bool bOnTheFlyNormals=false, const char* gridNormalOutputFile = NULL ) : CBaseData(std::string(name))
	{
		m_reconFilter = BCC_QUINTIC_BOX_SPLINE;
		m_normalReconFilter = BCC_QUINTIC_BOX_SPLINE;
		m_bOnTheFlyNormals = bOnTheFlyNormals;

		SetNormalEstimatorFilterByName( filter );

		if( gridNormalOutputFile){
			strcpy( m_gridNormalOutputFile, gridNormalOutputFile );
		}
		else
			m_gridNormalOutputFile[0] = 0;
	}


	virtual void Evaluate( const vector3_t min, const vector3_t max, const vector3_t step,
		elem_t (*fn)(const vector3_t p), void (*fnNormal)(const vector3_t p, vector3_t N) = 0 )
		;
	CBCCData( const char* name, const char* fname, double xScale, double yScale, double zScale, NORMAL_DISCRETE_FILTER normalEst = BCC_NORMAL_BCD, bool bOnTheFlyNormals = false, const char* gridNormalOutputFile = NULL );
	CBCCData( const char* name, const char* fname, DATA_FORMAT fmt, double xScale, double yScale, double zScale, const char* filter, bool bOnTheFlyNormals=false, const char* gridNormalOutputFile = NULL );
	CBCCData( const char* name, const char* fname, DATA_FORMAT fmt, const char *xname, const char *yname, const char *zname, double xScale, double yScale, double zScale, const char* gridNormalOutputFile = NULL );

/*
	CBCCData( const char* name, int xSize, int ySize, int zSize, double xScale, double yScale, double zScale, atomic_t *pData, NORMAL_DISCRETE_FILTER normalEst = BCC_NORMAL_BCD, const char *gridNormalOutputFile = NULL ): CBaseData(std::string(name)), m_normalEst(normalEst)
	{
		m_reconFilter = BCC_QUINTIC_BOX_SPLINE;
		m_normalReconFilter = BCC_QUINTIC_BOX_SPLINE;

		if( gridNormalOutputFile){
			strcpy( m_gridNormalOutputFile, gridNormalOutputFile );
		}
		else
			m_gridNormalOutputFile[0] = 0;

		m_halfScale[X] = xScale / 2.0f;
		m_halfScale[Y] = yScale / 2.0f;
		m_halfScale[Z] = zScale / 1.0f;

		CreateData( xSize, ySize, zSize, xScale, yScale, zScale, pData );

	}
*/
	virtual ~CBCCData()
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

	virtual void  CalculateNormals( void (*fnNormal)(const vector3_t p, vector3_t N) = 0 );
	virtual double Get( const vector3_t p ) const
	{
		
		//printf("%f %f %f start\n", p[0],p[1],p[2]);
		vector3_t bcc;
		WorldToBCC( bcc, p );

		switch( m_reconFilter ){
			case BCC_QUINTIC_BOX_SPLINE:
				//return Interpolate(bcc);
				return InterpolateQuinticBoxSplineTurbo( bcc );
				//return InterpolateQuinticBoxSplineSteroid(bcc);
				break;
			case BCC_LINEAR_BOX_SPLINE:
				return InterpolateLinearBoxSplineTurbo( bcc );
				break;		  
		}

		return 0;

	}
	
	virtual void   GetNormal( const vector3_t p, vector3_t N ) const
	{		
		if( m_bEpsilonNormalMode ){
			CBaseData::GetNormal( p, N );
			return;
		}
		
		vector3_t bcc;
		WorldToBCC( bcc, p );
		
		switch(m_normalReconFilter){ 
			case BCC_LINEAR_BOX_SPLINE:
				return InterpolateNormalLinearBoxSplineTurbo( bcc, N );
				break;
			case BCC_QUINTIC_BOX_SPLINE:
				return InterpolateNormalQuinticBoxSplineTurbo( bcc, N );
				break;
		}
	}
	
	virtual bool   SetNormalEstimatorFilterByName( const char* filter );
	virtual bool   SetReconstructionFilterByName( const char* filter );
	virtual bool   SetNormalReconstructionFilterByName( const char *filter );

	virtual void HookGradientComponentEstimators( );
};

class CBCCGradientComponentEstimator :public CBaseGradientComponentEstimator{
	protected: CBCCData *m_pData;
	public: CBCCGradientComponentEstimator( CBCCData *pData )
	: m_pData( pData ){}
};

class CBCCGradientXEstimatorPrecomp : public CBCCGradientComponentEstimator{
public:
	CBCCGradientXEstimatorPrecomp(CBCCData *pData )
			: CBCCGradientComponentEstimator( pData)
	{
		printf("\tCBCCGradientXEstimatorPrecomp Initialized\n");
	}
	virtual elem_t GetGradientComponent( int x, int y, int z ){
//		if( x < 0 || x >= m_pData->m_xSize ||
// 		  y < 0 || y >= m_pData->m_ySize ||
// 		  z < 0 || z >= m_pData->m_zSize ) return 0;
		return m_pData->m_pppNormals[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND][0];
	}
};

class CBCCGradientYEstimatorPrecomp : public CBCCGradientComponentEstimator{
public:
	CBCCGradientYEstimatorPrecomp(CBCCData *pData )
		: CBCCGradientComponentEstimator( pData)
	{
		printf("\tCBCCGradientYEstimatorPrecomp Initialized\n");
	}
	virtual elem_t GetGradientComponent( int x, int y, int z ){
//		if( x < 0 || x >= m_pData->m_xSize ||
//		    y < 0 || y >= m_pData->m_ySize ||
//		    z < 0 || z >= m_pData->m_zSize ) return 0;
		return m_pData->m_pppNormals[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND][1];
	}
};

class CBCCGradientZEstimatorPrecomp : public CBCCGradientComponentEstimator{
public:
	CBCCGradientZEstimatorPrecomp(CBCCData *pData )
		: CBCCGradientComponentEstimator( pData)
	{
		printf("\tCBCCGradientZEstimatorPrecomp Initialized\n");
	}
	virtual elem_t GetGradientComponent( int x, int y, int z ){
//		if( x < 0 || x >= m_pData->m_xSize ||
//		    y < 0 || y >= m_pData->m_ySize ||
//		    z < 0 || z >= m_pData->m_zSize ) return 0;
		return m_pData->m_pppNormals[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND][2];
	}
};

#define BCCESTIMATORCLASSNAME(component,name) CBCCGradient##component##Estimator##name
#define STR(x) #x

#define BCCESTIMATOR(component,name) \
class BCCESTIMATORCLASSNAME(component,name) : public CBCCGradientComponentEstimator{ 	\
	public:	 								     	\
		BCCESTIMATORCLASSNAME(component,name) (CBCCData *pData ) 		\
	: CBCCGradientComponentEstimator( pData)					\
		{									\
			printf("\t" STR(BCCESTIMATORCLASSNAME(component,name)) " Initialized\n");\
		}									\
		virtual elem_t GetGradientComponent( int x, int y, int z )		\
		{									\
			return m_pData->ComputeNormal##component##name( x,y,z);	\
		}									\
};


// class CBCCGradientXEstimatorSOCD : public CBCCGradientComponentEstimator{
// public:
// 	CBCCGradientXEstimatorSOCD(CBCCData *pData )
// 	: CBCCGradientComponentEstimator( pData)
// 	{
// 		printf("\tCBCCGradientXEstimatorSOCD Initialized\n");
// 	}
// 	virtual elem_t GetGradientComponent( int x, int y, int z )
// 	{
// 		return m_pData->ComputeNormalXCentralDifference( x,y,z);
// 	}
// };

// class CBCCGradientYEstimatorSOCD : public CBCCGradientComponentEstimator{
// public:
// 	CBCCGradientYEstimatorSOCD(CBCCData *pData )
// 	: CBCCGradientComponentEstimator( pData)
// 	{
// 		printf("\tCBCCGradientYEstimatorSOCD Initialized\n");
// 	}
// 	virtual elem_t GetGradientComponent( int x, int y, int z )
// 	{
// 		return m_pData->ComputeNormalYCentralDifference( x,y,z);
// 	}
// };
//
// class CBCCGradientZEstimatorSOCD : public CBCCGradientComponentEstimator{
// public:
// 	CBCCGradientZEstimatorSOCD(CBCCData *pData )
// 	: CBCCGradientComponentEstimator( pData)
// 	{
// 		printf("\tCBCCGradientZEstimatorSOCD Initialized\n");
// 	}
// 	virtual elem_t GetGradientComponent( int x, int y, int z )
// 	{
// 		return m_pData->ComputeNormalZCentralDifference( x,y,z);
// 	}
// };

///// SOCD /////////////////////////////////////////////////////////////
BCCESTIMATOR(X,CentralDifference)
BCCESTIMATOR(Y,CentralDifference)
BCCESTIMATOR(Z,CentralDifference)


//BCD//////////////////////////////////////////////////////////////////
BCCESTIMATOR(X,BoxDifference)
BCCESTIMATOR(Y,BoxDifference)
BCCESTIMATOR(Z,BoxDifference)

//OPT16////////////////////////////////////////////////////////////////
BCCESTIMATOR(X,SizeOptimalFilter)
BCCESTIMATOR(Y,SizeOptimalFilter)
BCCESTIMATOR(Z,SizeOptimalFilter)

//OPT26////////////////////////////////////////////////////////////////
BCCESTIMATOR(X,OptimalFilter)
BCCESTIMATOR(Y,OptimalFilter)
BCCESTIMATOR(Z,OptimalFilter)
#endif
