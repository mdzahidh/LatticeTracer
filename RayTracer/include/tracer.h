#ifndef __TRACER_H__
#define __TRACER_H__

#include "utils.h"
#include "data.h"
#include "volume.h"

enum{
	RT_FALSE = 0x0,
	RT_TRUE  = 0x1
};

enum{
	RT_LIGHTING           = 0x1,
	RT_SHADOW     	      = 0x2,
	RT_REFLECTION 	      = 0x4,
	RT_TWO_SIDED_LIGHTING = 0x8,
} ;

enum{
	RT_DVR = 0x0,
	RT_ISR = 0x1,
};

#define MAX_LIGHTS 8
#define MAX_DATA   8

typedef double (*fnAlphaTransfer)( const double data, const vector3_t xyz);
typedef void  (*fnColorTransfer)( material_t *mat, const double data, const vector3_t xyz, const vector3_t normal);
typedef void  (*fnDataCallback)( int x, int y, const vector3_t p, int n, double *pValues, vector3_t *pNormals);

void rtSetSamplingRate( double rate );
int  rtAddData( CBaseData *pData );
int  rtTraceRay( int x, int y, vector4_t color, const Ray_t *pRay, vector4_t background, int bounce = 0 );
void rtSetColorTransfer( fnColorTransfer fn );
void rtSetAlphaTransfer( fnAlphaTransfer fn );
int  rtAddLight( vector3_t light, color_t color );
void rtSetAmbientLight( color_t color );
void rtBegin();
void rtEnd();
void rtSet( int state, int value );
void rtSetDataCallback( fnDataCallback fn );
void rtSetData( CBaseData *pData );
void rtSetShadowOffsetParam( double p );
void rtSetShadowAmbient( color_t color );
void rtSetAmbient( color_t color );
void rtSetIsoSurface( double v );
void rtSetIsoSurfaceThreshold( double t );
void rtSetMode( int mode );
  		
#endif
