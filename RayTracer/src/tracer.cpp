#include "tracer.h"
#include "volume.h"
#include "data.h"
#include "camera.h"



// index 0 is the mother data grid. All the other data grid will be 
// will be evaluated with respect to the position of the mother data
// grid.
static CBaseData *g_pData[MAX_DATA]; 
static double     g_Values[MAX_DATA];
static vector3_t  g_Normals[MAX_DATA];


//static int		 g_xBound[MAX_DATA] = 0;
//static int		 g_yBound[MAX_DATA] = 0;
//static int		 g_zBound = 0;
static double	 g_sampleRate = 1;

static vector4_t g_lights[MAX_LIGHTS];
static color_t	 g_lightColors[MAX_LIGHTS];
//static color_t 	 g_ambientColor = {0.1f,0.1f,0.1f,0};
//static color_t   g_shadowAmbientColor = {0.25,0.25f,0.25f,0};

static color_t 	 g_ambientColor = {0,0,0,0};
static color_t   g_shadowAmbientColor = {0.3,0.3,0.3,0};

static int 		 g_activeLights = 0;
static vector3_t g_camPos;
static double    g_sampleDistance = 1;
static double    g_shadowOffsetParam = 2.0;
static double	 g_isoSurface = 0;
static double 	 g_isoSurfaceThreshold = EPSILON;
static unsigned int	 g_state = RT_LIGHTING | RT_SHADOW | RT_TWO_SIDED_LIGHTING;
static unsigned int  g_activeData  = 0;
static unsigned int  g_mode = RT_DVR;

fnAlphaTransfer g_pfnAlphaTransfer = 0;
fnColorTransfer g_pfnColorTransfer = 0;
fnDataCallback  g_pfnDataCallback  = 0;


static int numOfOnes( int value )
{
	if( !value ) return 0;
	int count = 0;
	while( value ) {
		value = value & ( value - 1 );
		++count;
	}
	return count;
}

double rtOpacityCorrection( double alpha, double sampleSpacing )
{
	return sqrt(1 - pow( (1-alpha), sampleSpacing ));
}

void rtSetMode( int m )
{
	g_mode = m;
}

void rtSetIsoSurfaceThreshold( double t )
{	
	g_isoSurfaceThreshold = t;	
}

void rtSetIsoSurface( double v )
{
	g_isoSurface = v;
}

void rtSetShadowAmbient( color_t color )
{
	g_shadowAmbientColor[R] = color[R];
	g_shadowAmbientColor[G] = color[G];
	g_shadowAmbientColor[B] = color[B];
	g_shadowAmbientColor[A] = color[A];
}

void rtSetAmbient( color_t color )
{
	g_ambientColor[R] = color[R];
	g_ambientColor[G] = color[G];
	g_ambientColor[B] = color[B];
	g_ambientColor[A] = color[A];
}

void rtSet( int state, int value )
{
  g_state = (g_state & ~state) | (value <<  numOfOnes(state-1));  
}

void rtSetDataCallback( fnDataCallback fn )
{
	g_pfnDataCallback = fn;
}

void rtBegin()
{  
	camGetPosition( g_camPos );
	g_sampleDistance = g_shadowOffsetParam / g_sampleRate;
}

void rtSetShadowOffsetParam( double p )
{
	g_shadowOffsetParam = p;
}

void rtEnd()
{
}

void rtSetAmbientLight( color_t color )
{
	glmCopyVector4d( color, g_ambientColor );
} 
		       
int rtAddLight( vector4_t light, color_t color )
{
	if( light[3] == 0.0 ){
		vector4_t tmpLight;
		glmCopyVector4d( light, tmpLight );
		glmNormalizeVector3d( tmpLight );
		glmCopyVector4d( tmpLight, g_lights[ g_activeLights] );
	}
	else{
		glmCopyVector4d( light, g_lights[ g_activeLights] );
	}
	
	glmCopyVector4d( color, g_lightColors[ g_activeLights] );
	return g_activeLights++;	
}

void rtSetSamplingRate( double rate )
{
	g_sampleRate = rate;
}

void rtSetData( CBaseData *pData )
{
	g_pData[ 0 ] = pData;
	g_activeData = 1;
}

int rtAddData( CBaseData *pData )
{
	g_pData[ g_activeData ]  = pData;

	return g_activeData++;
}


void rtSetColorTransfer( fnColorTransfer fn )
{
	g_pfnColorTransfer = fn;
}

void rtSetAlphaTransfer( fnAlphaTransfer fn )
{
	g_pfnAlphaTransfer = fn;
}


void glmMultVectorAdd3d( double out[], const double a[], const double b[], const double k )
{
	out[0] = a[0] + b[0] * k;
	out[1] = a[1] + b[1] * k;
	out[2] = a[2] + b[2] * k;
}



void glmMultVectorAdd4d( double out[], const double a[], const double b[], const double k )
{
	out[0] = a[0] + b[0] * k;
	out[1] = a[1] + b[1] * k;
	out[2] = a[2] + b[2] * k;
	out[3] = a[3] + b[3] * k;
}


void glmMultVectorByVector3d( double a[], const double b[] )
{
	a[X] *= b[X];
	a[Y] *= b[Y];
	a[Z] *= b[Z];
}

void glmMultVectorByVector4d( double a[], const double b[] )
{
	a[X] *= b[X];
	a[Y] *= b[Y];
	a[Z] *= b[Z];
	a[W] *= b[W];
}

double rtTraceShadowRayDirectional( const vector3_t dir, const vector3_t point, const vector3_t normal, const vector3_t traceRay  )
{
	vector3_t entry,exit;
	Ray_t ray;
	vector3_t startPoint;
	double offset = 0;
	
	double NdotDir = glmDotProduct3d( dir, normal );
	double NdotR  =  glmDotProduct3d(  normal, traceRay );

	if( (NdotDir > EPSILON) || (NdotDir < -EPSILON) ) {
	  offset = g_sampleDistance * NdotR / NdotDir;
	}
	
	if ( offset < 0 ) offset = -offset;

	 

	glmMultVectorAdd3d( startPoint, point, dir, offset );
	
	glmCopyVector3d( startPoint, ray.start );
	glmCopyVector3d( dir, ray.dir );
			
	if( !g_pData[0]->GetVolume()->FindEntryExitPoints( &ray, entry, exit ) ){		
		return 1;
	}

			
	vector3_t entryExit;
	glmSubtractVector3d( entryExit, exit, entry );
	double l = glmComputeVectorLength3d( entryExit );
	
	vector3_t dv;
	vector3_t p;
	
	glmCopyVector3d(ray.dir, dv);	
	glmMultVectorByScalar3d(dv,1.0f/g_sampleRate);	
	glmCopyVector3d( entry, p );

	unsigned int numSamples = (unsigned int)(l * g_sampleRate + 0.5) + 1;

	double     incomingAlpha = 1.0;	
	double compositeAlpha;
		
	while( numSamples-- ){
		
		if( incomingAlpha < EPSILON && incomingAlpha > -EPSILON ) {
			break;
		}
		
		//double interpData = g_pData->Get(xyz[X],xyz[Y],xyz[Z]);
		double interpData = g_pData[0]->Get(p);

		double currentAlpha = g_pfnAlphaTransfer( interpData, p );
		
		if( currentAlpha > EPSILON ){ // Remember Alpha is never negative
						
			material_t mat;
			
			//TODO: For now normal is set to be 0 for 
			// Shadow ray computation !
			vector3_t normal = {0,0,0};
			
			g_pfnColorTransfer( &mat, interpData, p, normal );		
		
			compositeAlpha = currentAlpha * incomingAlpha * mat.color[3];			
		
			incomingAlpha -= compositeAlpha;
		}
				
		glmAddVector3d( p, p, dv );		
	}
	
	return incomingAlpha;
};

double rtTraceShadowRayDirectionalIso( const vector3_t dir, const vector3_t point, const vector3_t normal, const vector3_t traceRay  )
{
	vector3_t entry,exit;	
	Ray_t ray;
	vector3_t startPoint;
	double offset = 0;
	
	double NdotDir = glmDotProduct3d( dir, normal );
	double NdotR  =  glmDotProduct3d(  normal, traceRay );

	if( (NdotDir > EPSILON) || (NdotDir < -EPSILON) ) {
		offset = g_sampleDistance * NdotR / NdotDir;
	}
	
	if ( offset < 0 ) offset = -offset;

	 

	glmMultVectorAdd3d( startPoint, point, dir, offset );
	
	glmCopyVector3d( startPoint, ray.start );
	glmCopyVector3d( dir, ray.dir );
	
	if( !g_pData[0]->GetVolume()->FindEntryExitPoints( &ray, entry, exit ) ){				
		return 1;
	}

	vector3_t entryExit;
	glmSubtractVector3d( entryExit, exit, entry );
	double l = glmComputeVectorLength3d( entryExit );
	
	vector3_t dv;
	vector3_t nextStartPoint,p;
	
	glmCopyVector3d(ray.dir, dv);	
	glmMultVectorByScalar3d(dv,1.0f/g_sampleRate);	
	glmCopyVector3d( entry, nextStartPoint );

	unsigned int numSamples = (unsigned int)(l * g_sampleRate + 0.5);
	
	double     incomingAlpha = 1.0;	
	double compositeAlpha = 0;	
		
	double currentData = 0;
	double currentDiff = 0;
	double nextData = 0;
	double nextDiff = 0;
	int    windowFound = 0;
	vector3_t nextp,midp;
				
	while( numSamples ){
		
		if( incomingAlpha < EPSILON && incomingAlpha > -EPSILON ) {
			break;
		}	
		
		glmCopyVector3d( nextStartPoint, p );
		
		windowFound = 0;
		
		while( numSamples ){
					
			--numSamples;
			
			currentData = g_pData[0]->Get(p);
			currentDiff = currentData - g_isoSurface;
			
					
			glmAddVector3d( nextp, p, dv );
			
			nextData = g_pData[0]->Get(nextp);
			nextDiff = nextData - g_isoSurface;
			
			if( (nextDiff * currentDiff) <= 0 ){
				//We found the bounding window
				windowFound = 1;
				glmCopyVector3d( nextp, nextStartPoint );
				break;							
			}
			
			glmCopyVector3d( nextp, p );
						
		}
		
	

		double alpha = 0;	
		
		if( windowFound ){
			
			double midData;
			double midDiff;
			
			if( ABS(currentDiff) < g_isoSurfaceThreshold ){
				// If the current position is the position
				// of the iso-surface, then just use this one.
				glmCopyVector3d( p, midp );
				midData = currentData;
			}
			else{
				glmAddVector3d( midp, p, nextp );
				glmMultVectorByScalar3d( midp, 0.5 );
				midData = g_pData[0]->Get(midp);
				midDiff =  midData - g_isoSurface;
					
				while( ABS(midDiff) > g_isoSurfaceThreshold ){ 
													
					if( (midDiff * nextDiff) <= 0 ){
						glmCopyVector3d( midp, p );
					}
					else{
						glmCopyVector3d( midp, nextp );			
					}
					
					glmAddVector3d( midp, p, nextp );
					glmMultVectorByScalar3d( midp, 0.5 );
					midData = g_pData[0]->Get(midp);
					midDiff = midData - g_isoSurface;
				}
			}
			
			material_t mat;
							
			glmAddVector3d( nextStartPoint, midp, dv );
			
			//vector3_t normal;						
			//g_pData[0]->GetNormal(midp, normal );
			
			alpha = g_pfnAlphaTransfer( midData, midp );
			g_pfnColorTransfer( &mat, g_isoSurface, midp, normal );		
			
																
	
			
			compositeAlpha = alpha * incomingAlpha * mat.color[3];
			incomingAlpha -= compositeAlpha;
			//glmMultVectorAdd3d( color, color, mat.color,  compositeAlpha );
		}
	}
	
	return incomingAlpha;
};

double rtTraceShadowRayPositional( const vector3_t light, const vector3_t point, const vector3_t normal, const vector3_t traceRay )
{
	vector3_t entry,exit;
	Ray_t ray;
	vector3_t endPoint;
	double offset = 0;

	
	glmCopyVector3d( light, ray.start );
	glmSubtractVector3d( ray.dir, point, light );
	glmNormalizeVector3d( ray.dir );
	
	if( !g_pData[0]->GetVolume()->FindEntryExitPoints( &ray, entry, exit ) ){		
		return 0;
	}

	double NdotDir = glmDotProduct3d( ray.dir, normal );
	double NdotR  =  glmDotProduct3d(  normal, traceRay );

	if( (NdotDir > EPSILON) || (NdotDir < -EPSILON) ) {
	  offset = g_sampleDistance * NdotR / NdotDir;
	}
	
	if ( offset < 0 ) offset = -offset;

	glmCopyVector3d( point, endPoint);
	glmMultVectorAdd3d( endPoint, endPoint, ray.dir, -offset ); 
			   
	vector3_t entryExit;
	glmSubtractVector3d( entryExit, endPoint, entry );
	double l = glmComputeVectorLength3d( entryExit );
	
	vector3_t dv;
	vector3_t p;
	
	glmCopyVector3d(ray.dir, dv);	
	glmMultVectorByScalar3d(dv,1.0f/g_sampleRate);	
	glmCopyVector3d( entry, p );

	unsigned int numSamples = (unsigned int)(l * g_sampleRate + 0.5) + 1;

	double     incomingAlpha = 1.0;
	
	double compositeAlpha;
		
	while( numSamples-- ){
		
		if( incomingAlpha < EPSILON && incomingAlpha > -EPSILON ) {
			break;
		}
		
		//double interpData = g_pData->Get(xyz[X],xyz[Y],xyz[Z]);
		double interpData = g_pData[0]->Get(p);

		double currentAlpha = g_pfnAlphaTransfer( interpData, p );
		
		if( currentAlpha > EPSILON ){ // Remember Alpha is never negative
						
			material_t mat;
			
			//TODO: For now normal is set to be 0 for 
			// Shadow ray computation !
			vector3_t normal = {0,0,0};
			
			g_pfnColorTransfer( &mat, interpData, p, normal );	

			compositeAlpha = currentAlpha * incomingAlpha * mat.color[3];			
		
			incomingAlpha -= compositeAlpha;
		}
				
		glmAddVector3d( p, p, dv );		
	}
		
	return incomingAlpha;
};

double rtTraceShadowRayPositionalIso( const vector3_t light, const vector3_t point, const vector3_t normal, const vector3_t traceRay  )
{
// 	vector3_t entry,exit;	
// 	Ray_t ray;
// 	vector3_t startPoint;
// 	double offset = 0;
// 	
// 	double NdotDir = glmDotProduct3d( dir, normal );
// 	double NdotR  =  glmDotProduct3d(  normal, traceRay );
// 
// 	if( (NdotDir > EPSILON) || (NdotDir < -EPSILON) ) {
// 		offset = g_sampleDistance * NdotR / NdotDir;
// 	}
// 	
// 	if ( offset < 0 ) offset = -offset;
// 
// 	 
// 
// 	glmMultVectorAdd3d( startPoint, point, dir, offset );
// 	
// 	glmCopyVector3d( startPoint, ray.start );
// 	glmCopyVector3d( dir, ray.dir );
// 	
// 	if( !g_pData[0]->GetVolume()->FindEntryExitPoints( &ray, entry, exit ) ){				
// 		return 1;
// 	}

	vector3_t entry,exit;
	Ray_t ray;
	vector3_t endPoint;
	double offset = 0;

	
	glmCopyVector3d( light, ray.start );
	glmSubtractVector3d( ray.dir, point, light );
	glmNormalizeVector3d( ray.dir );
	
	if( !g_pData[0]->GetVolume()->FindEntryExitPoints( &ray, entry, exit ) ){		
		return 0;
	}

	double NdotDir = glmDotProduct3d( ray.dir, normal );
	double NdotR  =  glmDotProduct3d(  normal, traceRay );

	if( (NdotDir > EPSILON) || (NdotDir < -EPSILON) ) {
		offset = g_sampleDistance * NdotR / NdotDir;
	}
	
	if ( offset < 0 ) offset = -offset;

	glmCopyVector3d( point, endPoint);
	glmMultVectorAdd3d( endPoint, endPoint, ray.dir, -offset );
	
	vector3_t entryExit;
	glmSubtractVector3d( entryExit, endPoint, entry );
	double l = glmComputeVectorLength3d( entryExit );
	
	vector3_t dv;
	vector3_t nextStartPoint,p;
	
	glmCopyVector3d(ray.dir, dv);	
	glmMultVectorByScalar3d(dv,1.0f/g_sampleRate);	
	glmCopyVector3d( entry, nextStartPoint );

	unsigned int numSamples = (unsigned int)(l * g_sampleRate + 0.5);
	
	double     incomingAlpha = 1.0;	
	double compositeAlpha = 0;	
	
	//double lastData = g_pData[0]->Get(p);
//	double lastData = g_pData[0]->Get(nextStartPoint);
//	double lastDiff = lastData - g_isoSurface; 
	
//	if( ABS(lastDiff) <= g_isoSurfaceThreshold ){
			
//	};
	
	double currentData = 0;
	double currentDiff = 0;
	double nextData = 0;
	double nextDiff = 0;
	int    windowFound = 0;
	vector3_t nextp,midp;
	
			
	while( numSamples ){
		
		if( incomingAlpha < EPSILON && incomingAlpha > -EPSILON ) {
			break;
		}	
		
		glmCopyVector3d( nextStartPoint, p );
		
		windowFound = 0;
		
		while( numSamples ){
					
			--numSamples;
			
			currentData = g_pData[0]->Get(p);
			currentDiff = currentData - g_isoSurface;
			
					
			glmAddVector3d( nextp, p, dv );
			
			nextData = g_pData[0]->Get(nextp);
			nextDiff = nextData - g_isoSurface;
			
			if( (nextDiff * currentDiff) <= 0 ){
				//We found the bounding window
				windowFound = 1;
				glmCopyVector3d( nextp, nextStartPoint );
				break;							
			}
			
			glmCopyVector3d( nextp, p );
						
		}
		
	

		double alpha = 0;	
		
		if( windowFound ){
			
			double midData;
			double midDiff;
			
			if( ABS(currentDiff) < g_isoSurfaceThreshold ){
				// If the current position is the position
				// of the iso-surface, then just use this one.
				glmCopyVector3d( p, midp );
				midData = currentData;
			}
			else{
				glmAddVector3d( midp, p, nextp );
				glmMultVectorByScalar3d( midp, 0.5 );
				midData = g_pData[0]->Get(midp);
				midDiff =  midData - g_isoSurface;
					
				while( ABS(midDiff) > g_isoSurfaceThreshold ){ 
													
					if( (midDiff * nextDiff) <= 0 ){
						glmCopyVector3d( midp, p );
					}
					else{
						glmCopyVector3d( midp, nextp );			
					}
					
					glmAddVector3d( midp, p, nextp );
					glmMultVectorByScalar3d( midp, 0.5 );
					midData = g_pData[0]->Get(midp);
					midDiff = midData - g_isoSurface;
				}
			}
			
			glmAddVector3d( nextStartPoint, midp, dv );
			
			material_t mat;
							
			//vector3_t normal;						
			//g_pData[0]->GetNormal(midp, normal );
			
			alpha = g_pfnAlphaTransfer( midData, midp );
			g_pfnColorTransfer( &mat, g_isoSurface, midp, normal );		
			
																
	
			
			compositeAlpha = alpha * incomingAlpha * mat.color[3];
			incomingAlpha -= compositeAlpha;
			//glmMultVectorAdd3d( color, color, mat.color,  compositeAlpha );
		}
	}
	
	return incomingAlpha;
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Deals with only lighting
// Note: color is only operated on RGB, and not on A. Light shouldn't be able to
// change opacity of a voxel.
//
// Lighting Model: 
// OutputColor = Material.Color ( sum( diffuse[i] * (Light[i].Color * Shadow[i] + ShadowAmbient) ) + Ambient )
// i = counting over the number of lights   
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void rtShade( material_t *mat, const vector3_t pos, const vector3_t normal, const vector3_t ray )
{
	
	color_t totalDiffuse = {0,0,0,0};
	color_t totalSpecular = {0,0,0,0};
	color_t totalAmbience = {0,0,0,0};
	color_t ambience;

	glmCopyVector3d( mat->color, ambience );
	//glmMultVectorByVector3d( ambience, g_ambientColor );
	glmMultVectorByVector3d( ambience, g_shadowAmbientColor );
	
	//glmMultAddVector3d( totalAmbience, totalAmbience, ambience);
	
	//Instead of starting from zero, we biased it
	totalAmbience[R] = mat->color[R] * g_ambientColor[R];
	totalAmbience[G] = mat->color[G] * g_ambientColor[G];
	totalAmbience[B] = mat->color[B] * g_ambientColor[B];
	totalAmbience[A] = mat->color[A] * g_ambientColor[A];
			
	for(int i=0;i<g_activeLights;++i){
		
		if( g_lights[i][3] > 0 ) {
			
			
			double shadow =1;
			if( g_state & RT_SHADOW)
			  shadow = rtTraceShadowRayPositional( g_lights[i], pos, normal, ray );

			double dSquare;
			vector3_t lightToP;

			glmSubtractVector3d( lightToP, pos, g_lights[i] );
			dSquare = glmDotProduct3d( lightToP, lightToP );
			glmMultVectorByScalar3d( lightToP, 1.0f/sqrt(dSquare) );
			
				
			//Positional light
			double diffuse = -glmDotProduct3d( normal, lightToP ) * 100000.0f / (dSquare + 1.0f);
			
			//NOTE: Doing diffuse for both backfacing and frontfacing normals !!!			
			if( g_state & RT_TWO_SIDED_LIGHTING ){
				if( diffuse < 0 ) diffuse = -diffuse;
			}
						

			if( diffuse > 0 ){
				vector4_t lightColor;
				
				glmCopyVector3d( g_lightColors[i], lightColor );
				glmMultVectorByScalar3d( lightColor, diffuse * shadow );
				glmAddVector3d( totalDiffuse, totalDiffuse, lightColor );

				glmMultVectorAdd3d( totalAmbience, totalAmbience, ambience, diffuse );
				
				if( shadow > EPSILON && mat->shininess > 0 ){
					
					//calculate specular component
					vector3_t posToEye;
					vector3_t reflectedVector;
					
					glmSubtractVector3d( posToEye, g_camPos, pos );
					glmNormalizeVector3d( posToEye );
					glmReflectVector3d( reflectedVector, lightToP, normal );
					
					double spec = glmDotProduct3d( reflectedVector, posToEye );
					
					if( spec > 0 ){
						color_t specColor = {1,1,1,0};
						spec = pow( spec, mat->shinePower ) * shadow * mat->shininess;
						glmMultVectorAdd3d( totalSpecular, totalSpecular, specColor, spec );
						
					}
				}

			}			
		}
		else{
			vector3_t shadowDir = {-g_lights[i][0], -g_lights[i][1], -g_lights[i][2] };
			
			double shadow = 1;
			
			if( g_state & RT_SHADOW ){			  
			  shadow = rtTraceShadowRayDirectional( shadowDir, pos, normal, ray );			  
			}
			
				
			//Directional light			
			double diffuse = -glmDotProduct3d( normal, g_lights[i] );			
			
			//NOTE: Doing diffuse for both backfacing and frontfacing normals !!!
			if( g_state & RT_TWO_SIDED_LIGHTING ){
				if( diffuse < 0 ) diffuse = -diffuse;
			}
			
			if( diffuse > 0 ){
				vector4_t lightColor;
				
				glmCopyVector3d( g_lightColors[i], lightColor );
				glmMultVectorByScalar3d( lightColor, diffuse * shadow );
				glmAddVector3d( totalDiffuse, totalDiffuse, lightColor );

				glmMultVectorAdd3d( totalAmbience, totalAmbience, ambience, diffuse );
				
				//glmAddVector3d( totalAmbience, totalAmbience, ambience );
				
				if( shadow > EPSILON && mat->shininess > 0 ){
					
					//calculate specular component
					vector3_t posToEye;
					vector3_t reflectedVector;
					
					glmSubtractVector3d( posToEye, g_camPos, pos );
					glmNormalizeVector3d( posToEye );
					glmReflectVector3d( reflectedVector, g_lights[i], normal );
					
					double spec = glmDotProduct3d( reflectedVector, posToEye );
					
					if( spec > 0 ){
						color_t specColor = {1,1,1,0};
						spec = pow( spec, mat->shinePower ) * shadow * mat->shininess;
						glmMultVectorAdd3d( totalSpecular, totalSpecular, specColor, spec );
						
					}
				}

			}			
		}
	}
		
			       
	glmMultVectorByVector3d( mat->color, totalDiffuse );
	glmAddVector3d( mat->color, mat->color, totalSpecular );
	glmAddVector3d( mat->color, mat->color, totalAmbience );
}



static int rtTraceRayDirectVolume( int x, int y, vector4_t color, const Ray_t *pRay, vector4_t background, int bounce )
{
	vector3_t entry,exit;
	if( !g_pData[0]->GetVolume()->FindEntryExitPoints( pRay, entry, exit ) ){
		color[X] = 0; color[Y] = 0; color[Z] = 0; color[W] = 0;		
		return 0;
	}

	vector3_t entryExit;
	glmSubtractVector3d( entryExit, exit, entry );
	double l = glmComputeVectorLength3d( entryExit );
	
	vector3_t dv;
	vector3_t p;
	double    sampleSpacing = 1.0f/g_sampleRate;
	
	glmCopyVector3d(pRay->dir, dv);	
	glmMultVectorByScalar3d(dv,sampleSpacing);	
	glmCopyVector3d( entry, p );

	unsigned int numSamples = (unsigned int)(l * g_sampleRate + 0.5) + 1;

	int 	   hitCount = 0;
	double     incomingAlpha = 1.0;

	color[0] = color[1] = color[2] = color[3] = 0;
	

	double compositeAlpha;	
	
	while( numSamples-- ){
		
		if( incomingAlpha < EPSILON && incomingAlpha > -EPSILON ) {
			break;
		}
				
		double interpData = g_pData[0]->Get(p);

		double currentAlpha = g_pfnAlphaTransfer( interpData, p );
		
		
		
		if( currentAlpha > EPSILON ){ // Remember Alpha is never negative
			
			//printf("%f\n",sampleSpacing);
			// Opacity correction
			currentAlpha = rtOpacityCorrection( currentAlpha, sampleSpacing );
			//////
			
			++hitCount;
			
			material_t mat;
			
			vector3_t normal;						
			g_pData[0]->GetNormal(p, normal );
			
			g_pfnColorTransfer( &mat, interpData, p, normal );		
		
							
			// We only iterate through the other grids when we have data callback
			if( g_pfnDataCallback ){
				
				// Copying the interpolated data and normal of the mother grid
				// to index 0
				g_Values[0] = interpData;
				glmCopyVector3d( normal, g_Normals[0] );


				// Here we will iterate through all the other data grids
				for(unsigned int k = 1; k < g_activeData; ++k ){
					
					g_Values[k] = g_pData[k]->Get(p);
					g_pData[k]->GetNormal( p, g_Normals[k] );
				}

				g_pfnDataCallback( x,y, p, g_activeData, g_Values, g_Normals );
			}

			glmNormalizeVector3d( normal );
			
			

			//////
			// Shading part comes here on cColor
			//////
			
			if(g_state & RT_LIGHTING)
			  rtShade( &mat, p, normal, pRay->dir );

			// TODO: Bounce is not working
			if( bounce < 0){
				if( (1 - currentAlpha) < EPSILON ){
					//Bounce rays
					Ray_t reflectedRay;
					glmCopyVector3d(p, reflectedRay.start);

					if( glmDotProduct3d(normal,pRay->dir) < 0)
						glmMultVectorByScalar3d(normal,-1);

					glmReflectVector3d( reflectedRay.dir, pRay->dir, normal );

					color_t refColor;
					
					rtTraceRay( x,y,refColor, &reflectedRay, background, bounce+1 );

					glmInterpolateVector3d( mat.color, refColor, mat.color, 0.9f );

				}
			}	
			
			compositeAlpha = currentAlpha * incomingAlpha * mat.color[3];	
			glmMultVectorAdd3d( color, color, mat.color, compositeAlpha );
		
			incomingAlpha -= compositeAlpha;
		}
		
		glmAddVector3d( p, p, dv );
	}

	if( incomingAlpha > 0 ){
		glmMultVectorAdd3d( color, color, background, incomingAlpha );
	}

	return hitCount;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Deals with only lighting
// Note: color is only operated on RGB, and not on A. Light shouldn't be able to
// change opacity of a voxel.
//
// Lighting Model: 
// OutputColor = Material.Color ( sum( diffuse[i] * (Light[i].Color * Shadow[i] + ShadowAmbient) ) + Ambient )
// i = counting over the number of lights   
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void rtShadeIso( material_t *mat, const vector3_t pos, const vector3_t normal, const vector3_t ray )
{
	
	color_t totalDiffuse = {0,0,0,0};
	color_t totalSpecular = {0,0,0,0};
	color_t totalAmbience = {0,0,0,0};
	color_t ambience;

	glmCopyVector3d( mat->color, ambience );
	//glmMultVectorByVector3d( ambience, g_ambientColor );
	glmMultVectorByVector3d( ambience, g_shadowAmbientColor );
	
	//glmMultAddVector3d( totalAmbience, totalAmbience, ambience);
	
	//Instead of starting from zero, we biased it
	totalAmbience[R] = mat->color[R] * g_ambientColor[R];
	totalAmbience[G] = mat->color[G] * g_ambientColor[G];
	totalAmbience[B] = mat->color[B] * g_ambientColor[B];
	totalAmbience[A] = mat->color[A] * g_ambientColor[A];
			
	for(int i=0;i<g_activeLights;++i){
		
		if( g_lights[i][3] > 0 ) {
			
			
			double shadow =1;
			if( g_state & RT_SHADOW)
				shadow = rtTraceShadowRayPositionalIso( g_lights[i], pos, normal, ray );

			double dSquare;
			vector3_t lightToP;

			glmSubtractVector3d( lightToP, pos, g_lights[i] );
			dSquare = glmDotProduct3d( lightToP, lightToP );
			glmMultVectorByScalar3d( lightToP, 1.0f/sqrt(dSquare) );
			
				
			//Positional light
			double diffuse = -glmDotProduct3d( normal, lightToP ) * 100000.0f / (dSquare + 1.0f);
			
			//NOTE: Doing diffuse for both backfacing and frontfacing normals !!!			
			if( g_state & RT_TWO_SIDED_LIGHTING ){
				if( diffuse < 0 ) diffuse = -diffuse;
			}
						
			
			if( diffuse > 0 ){
				vector4_t lightColor;
				
				glmCopyVector3d( g_lightColors[i], lightColor );
				glmMultVectorByScalar3d( lightColor, diffuse * shadow );
				glmAddVector3d( totalDiffuse, totalDiffuse, lightColor );

				glmMultVectorAdd3d( totalAmbience, totalAmbience, ambience, diffuse );
				
				if( shadow > EPSILON && mat->shininess > 0 ){
					
					//calculate specular component
					vector3_t posToEye;
					vector3_t reflectedVector;
					
					glmSubtractVector3d( posToEye, g_camPos, pos );
					glmNormalizeVector3d( posToEye );
					glmReflectVector3d( reflectedVector, lightToP, normal );
					
					double spec = glmDotProduct3d( reflectedVector, posToEye );
					
					if( spec > 0 ){
						color_t specColor = {1,1,1,0};
						spec = pow( spec, mat->shinePower ) * shadow * mat->shininess;
						glmMultVectorAdd3d( totalSpecular, totalSpecular, specColor, spec );
						
					}
				}

			}			
		}
		else{
			vector3_t shadowDir = {-g_lights[i][0], -g_lights[i][1], -g_lights[i][2] };
			
			double shadow = 1;
			
			if( g_state & RT_SHADOW ){			  
				shadow = rtTraceShadowRayDirectionalIso( shadowDir, pos, normal, ray );			  
			}
			
				
			//Directional light			
			double diffuse = -glmDotProduct3d( normal, g_lights[i] );			
			
			//NOTE: Doing diffuse for both backfacing and frontfacing normals !!!
			if( g_state & RT_TWO_SIDED_LIGHTING ){
				if( diffuse < 0 ) diffuse = -diffuse;
			}
			
			
			if( diffuse > 0 ){
				vector4_t lightColor;
				
				glmCopyVector3d( g_lightColors[i], lightColor );
				glmMultVectorByScalar3d( lightColor, diffuse * shadow );
				glmAddVector3d( totalDiffuse, totalDiffuse, lightColor );

				glmMultVectorAdd3d( totalAmbience, totalAmbience, ambience, diffuse );
				
				//glmAddVector3d( totalAmbience, totalAmbience, ambience );
				
				if( shadow > EPSILON && mat->shininess > 0 ){
					
					//calculate specular component
					vector3_t posToEye;
					vector3_t reflectedVector;
					
					glmSubtractVector3d( posToEye, g_camPos, pos );
					glmNormalizeVector3d( posToEye );
					glmReflectVector3d( reflectedVector, g_lights[i], normal );
					
					double spec = glmDotProduct3d( reflectedVector, posToEye );
					
					if( spec > 0 ){
						color_t specColor = {1,1,1,0};
						spec = pow( spec, mat->shinePower ) * shadow * mat->shininess;
						glmMultVectorAdd3d( totalSpecular, totalSpecular, specColor, spec );
						
					}
				}

			}			
		}
	}		
			       
	glmMultVectorByVector3d( mat->color, totalDiffuse );
	glmAddVector3d( mat->color, mat->color, totalSpecular );
	glmAddVector3d( mat->color, mat->color, totalAmbience );
}

static int rtTraceRayIsoSurface( int x, int y, vector4_t color, const Ray_t *pRay, vector4_t background, int bounce )
{
	
	vector3_t entry,exit;
	if( !g_pData[0]->GetVolume()->FindEntryExitPoints( pRay, entry, exit ) ){
		color[X] = 0; color[Y] = 0; color[Z] = 0; color[W] = 0;		
		return 0;
	}

	vector3_t entryExit;
	glmSubtractVector3d( entryExit, exit, entry );
	double l = glmComputeVectorLength3d( entryExit );
	
	vector3_t dv;
	vector3_t nextStartPoint,p;
	
	glmCopyVector3d(pRay->dir, dv);	
	glmMultVectorByScalar3d(dv,1.0f/g_sampleRate);	
	glmCopyVector3d( entry, nextStartPoint );

	unsigned int numSamples = (unsigned int)(l * g_sampleRate + 0.5);

	int 	   hitCount = 0;
	double     incomingAlpha = 1.0;

	color[0] = color[1] = color[2] = color[3] = 0;
	

	double compositeAlpha = 0;	
	
	//double lastData = g_pData[0]->Get(p);
//	double lastData = g_pData[0]->Get(nextStartPoint);
//	double lastDiff = lastData - g_isoSurface; 
	
//	if( ABS(lastDiff) <= g_isoSurfaceThreshold ){
			
//	};
	
	double currentData = 0;
	double currentDiff = 0;
	double nextData = 0;
	double nextDiff = 0;
	int    windowFound = 0;
	vector3_t nextp,midp;
	
	
	
	
	while( numSamples ){
		
		if( incomingAlpha < EPSILON && incomingAlpha > -EPSILON ) {
			break;
		}	
		
		glmCopyVector3d( nextStartPoint, p );
		
		windowFound = 0;
		
		while( numSamples ){
					
			--numSamples;
			
			currentData = g_pData[0]->Get(p);
			currentDiff = currentData - g_isoSurface;
			
					
			glmAddVector3d( nextp, p, dv );
			
			nextData = g_pData[0]->Get(nextp);
			nextDiff = nextData - g_isoSurface;
			
			if( (nextDiff * currentDiff) <= 0 ){
				//We found the bounding window
				windowFound = 1;
				glmCopyVector3d( nextp, nextStartPoint );
				break;							
			}
			
			glmCopyVector3d( nextp, p );
						
		}
		
	

		double alpha = 0;	
		
		if( windowFound ){
			double midData;
			double midDiff;
			
			if( ABS(currentDiff) < g_isoSurfaceThreshold ){
				// If the current position is the position
				// of the iso-surface, then just use this one.
				glmCopyVector3d( p, midp );
				midData = currentData;
			}
			else{	
				glmAddVector3d( midp, p, nextp );
				glmMultVectorByScalar3d( midp, 0.5 );
				midData = g_pData[0]->Get(midp);
				midDiff =  midData - g_isoSurface;
					
				while( ABS(midDiff) > g_isoSurfaceThreshold ){ 
													
					if( (midDiff * nextDiff) <= 0 ){
						glmCopyVector3d( midp, p );
					}
					else{
						glmCopyVector3d( midp, nextp );			
					}
					
					glmAddVector3d( midp, p, nextp );
					glmMultVectorByScalar3d( midp, 0.5 );
					midData = g_pData[0]->Get(midp);
					midDiff = midData - g_isoSurface;
				}
			}
						
			glmAddVector3d( nextStartPoint, midp, dv );
			
			++hitCount;
			
			material_t mat;
				
			vector3_t normal;						
			g_pData[0]->GetNormal(midp, normal );
			
			alpha = g_pfnAlphaTransfer( midData, midp );
			g_pfnColorTransfer( &mat, g_isoSurface, midp, normal );		
			
								
				// We only iterate through the other grids when we have data callback
			if( g_pfnDataCallback ){
					
					// Copying the interpolated data and normal of the mother grid
					// to index 0
				g_Values[0] = midData;
				glmCopyVector3d( normal, g_Normals[0] );
	
	
					// Here we will iterate through all the other data grids
				for(unsigned int k = 1; k < g_activeData; ++k ){
						
					g_Values[k] = g_pData[k]->Get(midp);
					g_pData[k]->GetNormal( midp, g_Normals[k] );
				}
	
				g_pfnDataCallback( x,y, midp, g_activeData, g_Values, g_Normals );
			}
	
			glmNormalizeVector3d( normal );
				
				
	
			//////
				// Shading part comes here on cColor
			//////
				
			if(g_state & RT_LIGHTING)
				rtShadeIso( &mat, midp, normal, pRay->dir );
			
			compositeAlpha = alpha * incomingAlpha * mat.color[3];
			incomingAlpha -= compositeAlpha;
			glmMultVectorAdd3d( color, color, mat.color,  compositeAlpha );
		}
	}
	
	if( incomingAlpha > 0 )
		glmMultVectorAdd3d( color, color, background, incomingAlpha );
	
	return hitCount;
};

int rtTraceRay( int x, int y, vector4_t color, const Ray_t *pRay, vector4_t background, int bounce )
{
	switch( g_mode ){
		case RT_DVR:
			return rtTraceRayDirectVolume( x,y,color,pRay,background,bounce );			
		case RT_ISR:
			return rtTraceRayIsoSurface( x,y,color,pRay,background,bounce );					
	}
	return 0;		
}
