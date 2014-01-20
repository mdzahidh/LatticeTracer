#include "transferfunc.h"
#include <stdio.h>

std::map<std::string,fnAlphaTransfer> CTransferFunction::m_alphaMap;
std::map<std::string,fnColorTransfer> CTransferFunction::m_colorMap;

double Tri( double t )
{
	if( (t>=0) && (t < 0.5f) ) return 1 - 2*t;
	else if( t > -0.5 && t < 0 ) return 1 + 2*t;
	return 0;
}

double Rect( double t )
{
	if( (t > -0.5f) && (t < 0.5f) ) return 1;
	return 0;
}

double Step(double t )
{
	if( t >= 0 ) return 1;
	else return 0;
}

double Ramp(double t )
{
    if( (t > 0) && (t <= 1)){
        return t;
    }
    return 0;
}

void ColorTransferCommon( material_t *mat, const double data, const vector3_t xyz, const vector3_t norm )
{

  mat->color[R] = 0.607843137;
  mat->color[G] = 0.607843137;
  mat->color[B] = 0.901960784;


  mat->color[A] = 1;
  mat->shininess = 1;
  mat->shinePower = 20;

}

double AlphaHAM( const double t, const vector3_t p)
{

    if( p[X] < 0 && p[Y] > 0 ) return 0;

    double diff = t - 0.4;

	//printf("diff = %f\n",diff);

	if( (diff > -0.01) && ( diff < 0.01 ) ){
		//printf("Passed\n");
		return 1;
	}

	return 0;
}
void ColorHAM( material_t *mat, const double data, const vector3_t xyz, const vector3_t norm )
{

  	mat->color[R] = 0.5176470588235294118;
	mat->color[G] = 0.3372549019607843137;
	mat->color[B] = 0.1725490196078431373;

	mat->color[A] = 1;
	mat->shininess = 1;
	mat->shinePower = 20;

}

double AlphaHAMMulti( const double t, const vector3_t p)
{

	double alpha =
		        0.9  * (Step(t - (0.4 - 0.01)) - Step( t - (0.4 + 0.01))) +
				0.5  * (Step(t - (0.5 - 0.01)) - Step( t - (0.5 + 0.01))) +
				0.2 * (Step(t - (0.6 - 0.02)) - Step( t - (0.6 + 0.02)));
	return alpha;
}
void ColorHAMMulti( material_t *mat, const double data, const vector3_t xyz, const vector3_t norm )
{
	color_t colorTable[3] = {
		{0.5176470588235294118,0.3372549019607843137,0.1725490196078431373,1},
		{0.211764706,0.756862745,0.454901961,1},
		{0.749019608,0.51372549,0.929411765,1}
	};

	int idx;// = ((int) ((data - 0.39) / (0.1)) ) % 3;

	if( (Step(data - (0.4 - 0.01)) - Step( data - (0.4 + 0.01)))  > 0 ){
		idx = 0;
	}
	else if( (Step(data - (0.5 - 0.01)) - Step( data - (0.5 + 0.01)) ) > 0 ) {
		idx = 1;
	}
	else idx = 2;

	mat->color[0] = colorTable[idx][0];
	mat->color[1] = colorTable[idx][1];
	mat->color[2] = colorTable[idx][2];
	mat->color[3] = colorTable[idx][3];

	mat->shininess = 1;
	mat->shinePower = 20;

}

double AlphaML( const double t, const vector3_t p)
{

	//printf("%f\n",t);

  	if( fabs(p[X]) > 0.95 ||
  		fabs(p[Y]) > 0.95 ||
  		fabs(p[Z]) > 0.95 )
  		return 0;
 // //
 // //
   	double diff = t - 0.5;
 // //
 // // 	//printf("diff = %f\n",diff);
 // //
   	if( (diff > -0.01) && ( diff < 0.01 ) ){
 // // 		//printf("Passed\n");
   		return 0.1;
   	}

 	return 0.001;

//	return 1;
}

void ColorML( material_t *mat, const double data, const vector3_t xyz, const vector3_t norm )
{

	ColorTransferCommon( mat, data, xyz, norm );

}

double AlphaModel( const double t, const vector3_t p)
{
	return Step( t - 100 ) - Step( t - 110 );
	//return Step(t - 1200 );
}


double AlphaBunny( const double t, const vector3_t p) {
  return Step(t - 2000);
}


void ColorModel( material_t *mat, const double data, const vector3_t xyz, const vector3_t norm )
{

	ColorTransferCommon( mat, data, xyz, norm );

}

double AlphaMouse( const double t, const vector3_t p)
{
	return Step( t - 1000 );
}

void ColorMouse( material_t *mat, const double data, const vector3_t xyz, const vector3_t norm )
{
	mat->color[R] = 193.0/255.0;
	mat->color[G] = 206.0/255.0;
	mat->color[B] = 218.0/255.0;

	mat->color[A]   = 1;
	mat->shininess  = 0.2;
	mat->shinePower = 100;
}

void ColorBunny( material_t *mat, const double data, const vector3_t xyz, const vector3_t norm ) {
  mat->color[R] = 162.0/255.0;
  mat->color[G] = 151.0/255.0;
  mat->color[B] = 139.0/255.0;

  mat->color[A]   = 1;
  mat->shininess  = 0.3;
  mat->shinePower = 100;
}

double AlphaMLPDE( const double t, const vector3_t p )
{
	if( fabs(p[X]) > 0.95 ||
		   fabs(p[Y]) > 0.95 ||
		   fabs(p[Z]) > 0.95 )
		return 0;


	double diff = t - 0.5;

	//printf("diff = %f\n",diff);

	if( (diff > -0.01) && ( diff < 0.01 ) ){
		//printf("Passed\n");
		return 1;
	}

	return 0;
}

void ColorMLPDE( material_t *mat, const double data, const vector3_t p,
		 const vector3_t norm)
{
	ColorTransferCommon( mat, data, p, norm );
	double len = sqrt(norm[X] * norm[X] +
		     norm[Y] * norm[Y] +
				norm[Z] * norm[Z]);
	if( (len > 4*0.7) && (len < 4*1.22) ){
		mat->color[R] = 1;
		mat->color[G] = 0;
		mat->color[B] = 0;
	}
}


void ColorFish( material_t *mat, const double data, const vector3_t xyz, const vector3_t norm )
{
	mat->color[R] = 1.0;
	mat->color[G] = 1.0;
	mat->color[B] = 1.0;

	mat->color[A]   = 1;
	mat->shininess = 0.4;
	//mat->shininess  = 0.6;
	mat->shinePower = 100;
}


double AlphaUnity( const double t, const vector3_t p) { return 1.0; }
double AlphaRamp( const double t, const vector3_t p ) {return Ramp(t);}

double AlphaKarp( const double t, const vector3_t p)
{
    if( (t > 92) && (t < 94) ){ //Skin
	    return Tri( (t - 93)/(2) ) * 0.001;
    }
    else if( (t > 100) && (t<115) ){ //Bones
        return Tri( (t - 107.5)/15 ) * 0.05;
    }
    else if( (t > 125) && (t<135) ){ //Inner bones
        return Ramp( (t - 125) / 10.0 ) * 0.1;
    }

    return 0.0;
}

void ColorKarp( material_t *mat, const double t, const vector3_t xyz, const vector3_t norm )
{
    mat->color[A]   = 1;

    if( (t > 92) && (t < 94) ){ //Skin
        mat->color[R] = 1.0;
        mat->color[G] = 1.0;
        mat->color[B] = 1.0;
        mat->shininess = 1.0;
        mat->shinePower = 200;
    }
    else if( (t > 100) && (t<115) ){ //Bones
        mat->color[R] = 1.0;
        mat->color[G] = 0.6667;
        mat->color[B] = 0.498;
        mat->shininess = 1.0;
        mat->shinePower = 200;
    }
    else if( (t > 125) && (t<135) ){ //Inner bones
        mat->color[R] = 1.0;
        mat->color[G] = 0.3333;
        mat->color[B] = 1.0;
        mat->shininess = 1.0;
        mat->shinePower = 400;
    }
}

double AlphaAngio( const double t, const vector3_t p)
{
	if( t >= 75 ) return 1;
	return 0;
} 

double AlphaTooth( const double t, const vector3_t p)
{
	if( p[Y] > 0 ) return 0;
	return 1;
}

void ColorTooth( material_t *mat, const double t, const vector3_t xyz, const vector3_t norm )
{
	ColorTransferCommon( mat, t, xyz, norm );	
}

double AlphaSheep( const double t, const vector3_t p)
{
	if( p[Z] > 0.3 && p[X] < 0.2) return 0;
	
	//if( p[Z] > 0 ) return 0;
	
	//return Rect( (t - 153.0 )/10.0 );
	
	
	return 1;
	//return 1;
	//return Tri(t - 80.0);
	//if( t > 80 ) return 1;
	//return 0;
	//return Tri((t-100)*0.1) * 1.3e-4;
	
}

void ColorSheepGrad( material_t *mat, const double t, const vector3_t xyz, const vector3_t norm )
{
	double L =  sqrt(norm[X] * norm[X] + 
				norm[Y] * norm[Y] + 
				norm[Z] * norm[Z]) *  0.011428571;
				
	double lum = (L / 130.0);
	if (lum > 1.0) lum = 1.0;
	
	mat->color[A] = 1;	
	mat->color[R] = lum;
    mat->color[G] = lum;
    mat->color[B] = lum;
			
	mat->shininess = 1.0;
    mat->shinePower = 100;		
	
}

void ColorSheep( material_t *mat, const double t, const vector3_t xyz, const vector3_t norm )
{
	//ColorTransferCommon( mat, t, xyz, norm );
	
	mat->color[A] = 1;
	//if( xyz[Z] <  0.300 ) {
	mat->color[R] = 1.0;
    mat->color[G] = 0.5;
    mat->color[B] = 0.5;
    
	
	mat->shininess = 1.0;
    mat->shinePower = 100;		
	
// 	}
// 	else{
// 		mat->color[R] = 0.4;
//         mat->color[G] = 0.4;
//         mat->color[B] = 0.4;
//         mat->shininess = 1.0;
//         mat->shinePower = 100;		
// 	}
}

double AlphaPhantom( const double t, const vector3_t p)
{
	if( p[Y] < -0.2 ) return 0;	
	//return 1.0;
	//return Ramp( t );
	return Tri( (t - 0.3)/0.1 )*0.75 + Tri( (t - 0.2)/0.1 )*0.25;	
}

void ColorPhantom( material_t *mat, const double t, const vector3_t xyz, const vector3_t norm )
{
	mat->color[A] = 1;	
	mat->color[R] = 0.0;
    mat->color[G] = 0.0;
    mat->color[B] = 0.0;    	
	mat->shininess = 1.0;
    mat->shinePower = 100;		
	
	if( Rect( (t-1.0) / 0.1 ) > 0 ){
		mat->color[B] = 1.0;		
	}
	else if(Rect( (t-0.2) / 0.1 ) > 0 ){
		mat->color[G] = 1.0;
	}
	else if(Rect( (t-0.3) / 0.1 ) > 0 ){
		mat->color[R] = 1.0;
	}
}

double AlphaPhantomSphere( const double t, const vector3_t p)
{
	//return Tri( (t - 27.5)/5.0 );
	return Ramp( (t - 250)/5.0 ) + Step(t-255.0);
}

void CTransferFunction::InitTransferFunctionSystem()
{
	m_alphaMap["AlphaModel"] = AlphaModel;
	m_colorMap["ColorModel"] = ColorModel;


	m_alphaMap["AlphaML"]    = AlphaML;
	m_colorMap["ColorML"]    = ColorML;

	m_alphaMap["AlphaUnity"] = AlphaUnity;
	m_alphaMap["AlphaRamp"] = AlphaRamp;

	m_alphaMap["AlphaHAM"]   = AlphaHAM;
	m_colorMap["ColorHAM"]   = ColorHAM;

	m_alphaMap["AlphaMouse"] = AlphaMouse;
	m_colorMap["ColorMouse"] = ColorMouse;

	m_colorMap["ColorFish"] = ColorFish;

	m_alphaMap["AlphaBunny"] = AlphaBunny;
	m_colorMap["ColorBunny"] = ColorBunny;

	m_alphaMap["AlphaHAMMulti"] = AlphaHAMMulti;
	m_colorMap["ColorHAMMulti"] = ColorHAMMulti;

	m_alphaMap["AlphaMLPDE"] = AlphaMLPDE;
	m_colorMap["ColorMLPDE"] = ColorMLPDE;


    m_alphaMap["AlphaKarp"] = AlphaKarp;
    m_colorMap["ColorKarp"] = ColorKarp;
	
	m_alphaMap["AlphaAngio"] = AlphaAngio;
	
	m_alphaMap["AlphaSheep"] = AlphaSheep;
	m_colorMap["ColorSheep"] = ColorSheep;
	m_colorMap["ColorSheepGrad"] = ColorSheepGrad;
	
	m_alphaMap["AlphaTooth"] = AlphaTooth;
	m_colorMap["ColorTooth"] = ColorTooth;
	
	m_alphaMap["AlphaPhantom"] = AlphaPhantom;
	m_colorMap["ColorPhantom"] = ColorPhantom;
	
	m_alphaMap["AlphaPhantomSphere"] = AlphaPhantomSphere;
}

bool CTransferFunction::SetAlphaTransferByName( std::string name )
{
	fnAlphaTransfer fnAlpha = m_alphaMap[name];
	if( !fnAlpha ){
		printf("Alpha transfer by the name %s doesn't exist!\n",name.c_str());
		return false;
	}

	rtSetAlphaTransfer( fnAlpha );

	return true;
}

bool CTransferFunction::SetColorTransferByName( std::string name )
{
	fnColorTransfer fnColor = m_colorMap[name];
	if( !fnColor){
		printf("Color transfer by the name %s doesn't exist!\n",name.c_str());
		return false;
	}

	rtSetColorTransfer( fnColor );

	return false;
}


