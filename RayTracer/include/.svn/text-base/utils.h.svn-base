#ifndef __UTILS_H__
#define __UTILS_H__

enum{
	X = 0,
	Y,
 	Z,
  	W
};

enum{
	R = 0,
	G,
	B,
	A
};


typedef double  elem_t;
typedef elem_t vector2_t[2];
typedef elem_t vector3_t[3];
typedef elem_t vector4_t[4];
typedef elem_t matrix44_t[16];
typedef vector4_t color_t;

typedef struct{
	color_t color;
	double	shininess;  // 0..1
	double   shinePower; // exponent	
} material_t;


#define glmInitVector4d(out,x,y,z,w) do{ out[X] = x; out[Y] = y; out[Z] = z; out[W] = w;}while(0)
#define glmInitVector3d(out,x,y,z) do{ out[X] = x; out[Y] = y; out[Z] = z; }while(0)
		     		
#define DELETE_SAFE(x) if( x ) { delete x; x = NULL; };
#define DELETE_SAFE_ARRAY(x) if( x ) { delete [] x; x = NULL; };

#define PI			3.1415926535897932384626433832795
#define PI_OVER_2	1.5707963267948966192313216916398
#define PI_2		6.283185307179586476925286766559
#define PI_SQUARE   9.8696044010893586188344909998762
#define PI_CUBE     31.006276680299820175476315067101

#define DEG2RAD(x) (((x) / 180.0f) * PI)

#define EPSILON 5e-4
#define DELTA   1e-3

#define ABS(x) ( (x) < 0 ? -(x) : (x) )
#define STR(x) #x

#endif

