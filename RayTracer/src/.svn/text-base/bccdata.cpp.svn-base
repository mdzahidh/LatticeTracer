 #include <string.h>
#include "bccdata.h"
#include "camera.h"
#include <stdlib.h>

using namespace std;



//#define LINEAR
#define ODD(x) ((x)&1)
#define EVEN(x) !ODD(x)

#define BCC2CC(index3D,x,y,z) {BCCToIndex(index3D,x,y,z);}

#define FILL_PPIPED_STEROID(p,Q){ \
	BCC2CC(index3D,x0,y0,z0); \
	p[0] = Fetch( index3D[0], index3D[1], index3D[2] ); \
	BCC2CC(index3D,(x0-(Q[0])),(y0+(Q[1])),(z0+(Q[2]))) \
	p[1] = Fetch( index3D[0], index3D[1], index3D[2] ); \
	BCC2CC(index3D,(x0+(Q[0])),(y0-(Q[1])),(z0+(Q[2]))) \
	p[2] = Fetch( index3D[0], index3D[1], index3D[2] ); \
	BCC2CC(index3D,(x0+(Q[0])),(y0+(Q[1])),(z0-(Q[2]))) \
	p[3] = Fetch( index3D[0], index3D[1], index3D[2] ); \
	BCC2CC(index3D,(x0+2*(Q[0])),y0,z0) \
	p[4] = Fetch( index3D[0], index3D[1], index3D[2] ); \
	BCC2CC(index3D,x0,(y0+2*(Q[1])),z0) \
	p[5] = Fetch( index3D[0], index3D[1], index3D[2] ); \
	BCC2CC(index3D,x0,y0,(z0+2*(Q[2]))) \
	p[6] = Fetch( index3D[0], index3D[1], index3D[2] ); \
	BCC2CC(index3D,(x0+(Q[0])),(y0+(Q[1])),(z0+(Q[2]))) \
	p[7] = Fetch( index3D[0], index3D[1], index3D[2] ); \
}

#define FILL_PPIPED_NORM(p,a,b,c,maxdir,middir,mindir,ordinate) { \
                   BCCToIndex(index3D,x0,y0,z0); \
                   p[0]=GetValueNorm(index3D, ordinate);\
                   BCCToIndex(index3D,x0-(a),y0+(b),z0+(c)); \
                   p[maxdir]=GetValueNorm(index3D, ordinate);\
                   BCCToIndex(index3D,x0+(a),y0-(b),z0+(c)); \
                   p[middir]=GetValueNorm(index3D,ordinate);\
                   BCCToIndex(index3D,x0+(a),y0+(b),z0-(c)); \
                   p[mindir]=GetValueNorm(index3D,ordinate);\
                   BCCToIndex(index3D,x0+2*(a),y0,z0); \
                   p[3+maxdir]=GetValueNorm(index3D,ordinate);\
                   BCCToIndex(index3D,x0,y0+2*(b),z0); \
                   p[3+middir]=GetValueNorm(index3D,ordinate);\
                   BCCToIndex(index3D,x0,y0,z0+2*(c)); \
                   p[3+mindir]=GetValueNorm(index3D,ordinate);\
                   BCCToIndex(index3D,x0+(a),y0+(b),z0+(c)); \
                   p[7]=GetValueNorm(index3D,ordinate);\
                  }

#define FILL_PPIPED(p,a,b,c,maxdir,middir,mindir) { \
                   BCCToIndex(index3D,x0,y0,z0); \
                   p[0]=Fetch( index3D[0], index3D[1], index3D[2] );\
                   BCCToIndex(index3D,x0-(a),y0+(b),z0+(c)); \
                   p[maxdir]=Fetch( index3D[0], index3D[1], index3D[2] );\
                   BCCToIndex(index3D,x0+(a),y0-(b),z0+(c)); \
                   p[middir]=Fetch( index3D[0], index3D[1], index3D[2] );\
                   BCCToIndex(index3D,x0+(a),y0+(b),z0-(c)); \
                   p[mindir]=Fetch( index3D[0], index3D[1], index3D[2] );\
                   BCCToIndex(index3D,x0+2*(a),y0,z0); \
                   p[3+maxdir]=Fetch( index3D[0], index3D[1], index3D[2] );\
                   BCCToIndex(index3D,x0,y0+2*(b),z0); \
                   p[3+middir]=Fetch( index3D[0], index3D[1], index3D[2] );\
                   BCCToIndex(index3D,x0,y0,z0+2*(c)); \
                   p[3+mindir]=Fetch( index3D[0], index3D[1], index3D[2] );\
                   BCCToIndex(index3D,x0+(a),y0+(b),z0+(c)); \
                   p[7]=Fetch( index3D[0], index3D[1], index3D[2] );\
                  }


#define CONV_PPIPED_OLD(value, p,alpha,beta,gamma) \
   value += (-10.0f*p[0]+4.0f*(p[1]+p[2]+p[3])-2.0f*(p[4]+p[5]+p[6])+p[7]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*p[0]-2.0f*(p[2]+p[3])+p[4]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*p[0]-2.0f*(p[1]+p[3])+p[5]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*p[0]-2.0f*(p[1]+p[2])+p[6]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*p[0]+p[3]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*p[0]+p[2]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*p[0]+p[1]) * rho22(gamma, alpha-1, beta-1) + \
   (p[0]) * rho22(alpha-1, beta-1, gamma-1);

#define RHO_FAST(alpha,beta,gamma) -alpha3*(gamma*beta-0.5*alpha*(gamma+beta)+0.3 *alpha2);
#define CONV_PPIPED_OPT(value,p,alpha,beta,gamma)  \
    p123 = p[1] + p[2] + p[3]; \
    p0 = 4.0 * p[0]; \
    alpha2 = alpha * alpha;\
    alpha3 = one_over_6 * alpha2 * alpha;\
    base_rho = RHO_FAST(alpha,beta,gamma);\
    base_rho2 = base_rho - 0.5 * alpha3 * alpha; \
    value += (-2.5*p0+4.0*(p123)-2.0*(p[4]+p[5]+p[6])+p[7]) \
    * base_rho; \
    alpha2 = base_rho2 + alpha3 * beta; \
    base_rho2 += alpha3 * gamma; \
    value += \
    (p0-2.0*(p123-p[1])+p[4]) \
    * (alpha2) + \
    (p0-2.0*(p123-p[2])+p[5]) \
    * (base_rho2) + \
    (-.5*p0+p[3])  \
    * (base_rho2 + alpha2 - alpha3 - base_rho); \
    alpha -= 1.0; \
    alpha2 = beta * beta; \
    alpha3 = one_over_6 * alpha2 * beta; \
    base_rho =  RHO_FAST(beta, gamma, alpha); \
    value += (p0-2.0*(p123-p[3])+p[6]) \
    * (base_rho); \
    p0 = -.5 * p0; \
    value += (p0+p[2])\
    * (base_rho + alpha3 * (alpha - .5 * beta)) + \
    (p0+p[1]) \
    * rho22(gamma, alpha, beta-1.0) + \
    (-.5 * p0)\
    * rho22(alpha, beta-1.0, gamma-1.0);\


#define CONV_PPIPED(value,p,alpha,beta,gamma) { \
    register double p123 = p[1] + p[2] + p[3]; \
    register double p0 = 4.0f * p[0]; \
    register double alpha2 = alpha * alpha;\
    register double alpha3 = one_over_6 * alpha2 * alpha;\
    register double base_rho = RHO_FAST(alpha,beta,gamma);\
    register double base_rho2 = base_rho - 0.5f * alpha3 * alpha; \
    value += (-2.5f*p0+4.0f*(p123)-2.0f*(p[4]+p[5]+p[6])+p[7]) \
    * base_rho; \
    alpha2 = base_rho2 + alpha3 * beta; \
    base_rho2 += alpha3 * gamma; \
    value += \
    (p0-2.0f*(p123-p[1])+p[4]) \
    * (alpha2) + \
    (p0-2.0f*(p123-p[2])+p[5]) \
    * (base_rho2) + \
    (-.5f*p0+p[3])  \
    * (base_rho2 + alpha2 - alpha3 - base_rho); \
    alpha -= 1.0f; \
    alpha2 = beta * beta; \
    alpha3 = one_over_6 * alpha2 * beta; \
    base_rho =  RHO_FAST(beta, gamma, alpha); \
    value += (p0-2.0f*(p123-p[3])+p[6]) \
    * (base_rho); \
    p0 = -.5f * p0; \
    value += (p0+p[2])\
    * (base_rho + alpha3 * (alpha - .5f * beta)) + \
    (p0+p[1]) \
    * rho22(gamma, alpha, beta-1.0f) + \
    (-.5f * p0)\
    * rho22(alpha, beta-1.0f, gamma-1.0f);\
}


#define FILL_NORM(v,n) {v##x[n] = GetValueNorm(id,0);\
	v##y[n] = GetValueNorm(id,1);\
	v##z[n] = GetValueNorm(id,2);}


#define Norm_RegionOne_FillData_Even \
id[2]-=3;       FILL_NORM(p3,4);\
id[2]++;        FILL_NORM(p1,6);\
id[0]++;        FILL_NORM(p3,2);\
id[0]--;id[1]++;        FILL_NORM(p3,3);\
id[0]++;        FILL_NORM(p3,7);\
id[0]-=2;id[1]-=2;id[2]++;      FILL_NORM(p1,7);\
id[0]++;        FILL_NORM(p1,1);\
id[0]++;        FILL_NORM(p4,3);\
id[0]-=2;id[1]++;       FILL_NORM(p1,2);\
id[0]++;        FILL_NORM(p3,0);\
id[0]++;        FILL_NORM(p3,6);\
id[0]--;id[1]++;        FILL_NORM(p3,5);\
id[1]-=2;id[2]++;       FILL_NORM(p1,5);\
id[0]++;        FILL_NORM(p4,5);\
id[0]-=2;id[1]++;       FILL_NORM(p1,4);\
id[0]++;        FILL_NORM(p1,0);\
id[0]++;        FILL_NORM(p4,0);\
id[0]++;        FILL_NORM(p4,4);\
id[0]-=2;id[1]++;       FILL_NORM(p2,1);\
id[0]++;        FILL_NORM(p3,1);\
id[0]-=2;id[1]-=2;id[2]++;      FILL_NORM(p1,3);\
id[0]++;        FILL_NORM(p4,1);\
id[0]++;        FILL_NORM(p4,7);\
id[0]-=2;id[1]++;       FILL_NORM(p2,6);\
id[0]++;        FILL_NORM(p2,0);\
id[0]++;        FILL_NORM(p4,2);\
id[0]--;id[1]++;        FILL_NORM(p2,5);\
id[1]--;id[2]++;        FILL_NORM(p2,2);\
id[0]++;        FILL_NORM(p4,6);\
id[0]--;id[1]++;        FILL_NORM(p2,7);\
id[0]++;        FILL_NORM(p2,3);\
id[0]--;id[1]--;id[2]++;        FILL_NORM(p2,4);\

#define Norm_RegionOne_FillData_Odd \
id[0]++;id[1]++;id[2]-=3;       FILL_NORM(p3,4);\
id[0]--;id[1]--;id[2]++;        FILL_NORM(p1,6);\
id[0]++;        FILL_NORM(p3,2);\
id[0]--;id[1]++;        FILL_NORM(p3,3);\
id[0]++;        FILL_NORM(p3,7);\
id[0]--;id[1]--;id[2]++;        FILL_NORM(p1,7);\
id[0]++;        FILL_NORM(p1,1);\
id[0]++;        FILL_NORM(p4,3);\
id[0]-=2;id[1]++;       FILL_NORM(p1,2);\
id[0]++;        FILL_NORM(p3,0);\
id[0]++;        FILL_NORM(p3,6);\
id[0]--;id[1]++;        FILL_NORM(p3,5);\
id[0]--;id[1]-=3;id[2]++;       FILL_NORM(p1,5);\
id[0]++;        FILL_NORM(p4,5);\
id[0]-=2;id[1]++;       FILL_NORM(p1,4);\
id[0]++;        FILL_NORM(p1,0);\
id[0]++;        FILL_NORM(p4,0);\
id[0]++;        FILL_NORM(p4,4);\
id[0]-=2;id[1]++;       FILL_NORM(p2,1);\
id[0]++;        FILL_NORM(p3,1);\
id[0]--;id[1]--;id[2]++;        FILL_NORM(p1,3);\
id[0]++;        FILL_NORM(p4,1);\
id[0]++;        FILL_NORM(p4,7);\
id[0]-=2;id[1]++;       FILL_NORM(p2,6);\
id[0]++;        FILL_NORM(p2,0);\
id[0]++;        FILL_NORM(p4,2);\
id[0]--;id[1]++;        FILL_NORM(p2,5);\
id[0]--;id[1]-=2;id[2]++;       FILL_NORM(p2,2);\
id[0]++;        FILL_NORM(p4,6);\
id[0]--;id[1]++;        FILL_NORM(p2,7);\
id[0]++;        FILL_NORM(p2,3);\
id[2]++;        FILL_NORM(p2,4);\

/////////////////////////////////


/////////////////////////////////

#define Norm_RegionTwo_FillData_Even \
id[2]-=3;       FILL_NORM(p3,4);\
id[2]++;        FILL_NORM(p1,6);\
id[0]++;        FILL_NORM(p3,3);\
id[0]--;id[1]++;        FILL_NORM(p3,2);\
id[0]++;        FILL_NORM(p3,7);\
id[0]-=2;id[1]-=2;id[2]++;      FILL_NORM(p1,7);\
id[0]++;        FILL_NORM(p1,2);\
id[0]--;id[1]++;        FILL_NORM(p1,1);\
id[0]++;        FILL_NORM(p3,0);\
id[0]++;        FILL_NORM(p3,5);\
id[0]-=2;id[1]++;       FILL_NORM(p4,3);\
id[0]++;        FILL_NORM(p3,6);\
id[1]-=2;id[2]++;       FILL_NORM(p1,4);\
id[0]--;id[1]++;        FILL_NORM(p1,5);\
id[0]++;        FILL_NORM(p1,0);\
id[0]++;        FILL_NORM(p2,1);\
id[0]-=2;id[1]++;       FILL_NORM(p4,5);\
id[0]++;        FILL_NORM(p4,0);\
id[0]++;        FILL_NORM(p3,1);\
id[0]--;id[1]++;        FILL_NORM(p4,4);\
id[0]--;id[1]-=3;id[2]++;       FILL_NORM(p1,3);\
id[0]++;        FILL_NORM(p2,6);\
id[0]--;id[1]++;        FILL_NORM(p4,1);\
id[0]++;        FILL_NORM(p2,0);\
id[0]++;        FILL_NORM(p2,5);\
id[0]-=2;id[1]++;       FILL_NORM(p4,7);\
id[0]++;        FILL_NORM(p4,2);\
id[1]--;id[2]++;        FILL_NORM(p2,2);\
id[0]++;        FILL_NORM(p2,7);\
id[0]--;id[1]++;        FILL_NORM(p4,6);\
id[0]++;        FILL_NORM(p2,3);\
id[0]--;id[1]--;id[2]++;        FILL_NORM(p2,4);\

#define Norm_RegionTwo_FillData_Odd \
id[0]++;id[1]++;id[2]-=3;       FILL_NORM(p3,4);\
id[0]--;id[1]--;id[2]++;        FILL_NORM(p1,6);\
id[0]++;        FILL_NORM(p3,3);\
id[0]--;id[1]++;        FILL_NORM(p3,2);\
id[0]++;        FILL_NORM(p3,7);\
id[0]--;id[1]--;id[2]++;        FILL_NORM(p1,7);\
id[0]++;        FILL_NORM(p1,2);\
id[0]--;id[1]++;        FILL_NORM(p1,1);\
id[0]++;        FILL_NORM(p3,0);\
id[0]++;        FILL_NORM(p3,5);\
id[0]-=2;id[1]++;       FILL_NORM(p4,3);\
id[0]++;        FILL_NORM(p3,6);\
id[0]--;id[1]-=3;id[2]++;       FILL_NORM(p1,4);\
id[0]--;id[1]++;        FILL_NORM(p1,5);\
id[0]++;        FILL_NORM(p1,0);\
id[0]++;        FILL_NORM(p2,1);\
id[0]-=2;id[1]++;       FILL_NORM(p4,5);\
id[0]++;        FILL_NORM(p4,0);\
id[0]++;        FILL_NORM(p3,1);\
id[0]--;id[1]++;        FILL_NORM(p4,4);\
id[1]-=2;id[2]++;       FILL_NORM(p1,3);\
id[0]++;        FILL_NORM(p2,6);\
id[0]--;id[1]++;        FILL_NORM(p4,1);\
id[0]++;        FILL_NORM(p2,0);\
id[0]++;        FILL_NORM(p2,5);\
id[0]-=2;id[1]++;       FILL_NORM(p4,7);\
id[0]++;        FILL_NORM(p4,2);\
id[0]--;id[1]-=2;id[2]++;       FILL_NORM(p2,2);\
id[0]++;        FILL_NORM(p2,7);\
id[0]--;id[1]++;        FILL_NORM(p4,6);\
id[0]++;        FILL_NORM(p2,3);\
id[2]++;        FILL_NORM(p2,4);\

/////////////////////////////////


/////////////////////////////////

#define Norm_RegionThree_FillData_Even \
id[2]-=2;       FILL_NORM(p1,5);\
id[1]++;        FILL_NORM(p4,5);\
id[0]--;id[1]-=2;id[2]++;       FILL_NORM(p1,7);\
id[0]++;        FILL_NORM(p1,3);\
id[0]--;id[1]++;        FILL_NORM(p1,1);\
id[0]++;        FILL_NORM(p4,1);\
id[0]--;id[1]++;        FILL_NORM(p4,3);\
id[0]++;        FILL_NORM(p4,7);\
id[1]-=2;id[2]++;       FILL_NORM(p1,4);\
id[0]--;id[1]++;        FILL_NORM(p1,6);\
id[0]++;        FILL_NORM(p1,0);\
id[0]++;        FILL_NORM(p2,2);\
id[0]-=2;id[1]++;       FILL_NORM(p3,2);\
id[0]++;        FILL_NORM(p4,0);\
id[0]++;        FILL_NORM(p4,6);\
id[0]--;id[1]++;        FILL_NORM(p4,4);\
id[0]--;id[1]-=3;id[2]++;       FILL_NORM(p1,2);\
id[0]++;        FILL_NORM(p2,6);\
id[0]-=2;id[1]++;       FILL_NORM(p3,4);\
id[0]++;        FILL_NORM(p3,0);\
id[0]++;        FILL_NORM(p2,0);\
id[0]++;        FILL_NORM(p2,4);\
id[0]-=2;id[1]++;       FILL_NORM(p3,6);\
id[0]++;        FILL_NORM(p4,2);\
id[0]--;id[1]--;id[2]++;        FILL_NORM(p3,3);\
id[0]++;        FILL_NORM(p2,1);\
id[0]++;        FILL_NORM(p2,7);\
id[0]-=2;id[1]++;       FILL_NORM(p3,7);\
id[0]++;        FILL_NORM(p3,1);\
id[0]++;        FILL_NORM(p2,3);\
id[0]-=2;id[1]--;id[2]++;       FILL_NORM(p3,5);\
id[0]++;        FILL_NORM(p2,5);\

#define Norm_RegionThree_FillData_Odd \
id[2]-=2;       FILL_NORM(p1,5);\
id[1]++;        FILL_NORM(p4,5);\
id[1]--;id[2]++;        FILL_NORM(p1,7);\
id[0]++;        FILL_NORM(p1,3);\
id[0]--;id[1]++;        FILL_NORM(p1,1);\
id[0]++;        FILL_NORM(p4,1);\
id[0]--;id[1]++;        FILL_NORM(p4,3);\
id[0]++;        FILL_NORM(p4,7);\
id[0]--;id[1]-=3;id[2]++;       FILL_NORM(p1,4);\
id[0]--;id[1]++;        FILL_NORM(p1,6);\
id[0]++;        FILL_NORM(p1,0);\
id[0]++;        FILL_NORM(p2,2);\
id[0]-=2;id[1]++;       FILL_NORM(p3,2);\
id[0]++;        FILL_NORM(p4,0);\
id[0]++;        FILL_NORM(p4,6);\
id[0]--;id[1]++;        FILL_NORM(p4,4);\
id[1]-=2;id[2]++;       FILL_NORM(p1,2);\
id[0]++;        FILL_NORM(p2,6);\
id[0]-=2;id[1]++;       FILL_NORM(p3,4);\
id[0]++;        FILL_NORM(p3,0);\
id[0]++;        FILL_NORM(p2,0);\
id[0]++;        FILL_NORM(p2,4);\
id[0]-=2;id[1]++;       FILL_NORM(p3,6);\
id[0]++;        FILL_NORM(p4,2);\
id[0]-=2;id[1]-=2;id[2]++;      FILL_NORM(p3,3);\
id[0]++;        FILL_NORM(p2,1);\
id[0]++;        FILL_NORM(p2,7);\
id[0]-=2;id[1]++;       FILL_NORM(p3,7);\
id[0]++;        FILL_NORM(p3,1);\
id[0]++;        FILL_NORM(p2,3);\
id[0]--;id[2]++;        FILL_NORM(p3,5);\
id[0]++;        FILL_NORM(p2,5);\

/////////////////////////////////


/////////////////////////////////

#define Norm_RegionFour_FillData_Even \
id[2]-=2;       FILL_NORM(p1,5);\
id[0]++;        FILL_NORM(p4,5);\
id[0]-=2;id[1]--;id[2]++;       FILL_NORM(p1,7);\
id[0]++;        FILL_NORM(p1,1);\
id[0]++;        FILL_NORM(p4,3);\
id[0]-=2;id[1]++;       FILL_NORM(p1,3);\
id[0]++;        FILL_NORM(p4,1);\
id[0]++;        FILL_NORM(p4,7);\
id[0]--;id[1]--;id[2]++;        FILL_NORM(p1,6);\
id[0]++;        FILL_NORM(p3,2);\
id[0]-=2;id[1]++;       FILL_NORM(p1,4);\
id[0]++;        FILL_NORM(p1,0);\
id[0]++;        FILL_NORM(p4,0);\
id[0]++;        FILL_NORM(p4,4);\
id[0]-=2;id[1]++;       FILL_NORM(p2,2);\
id[0]++;        FILL_NORM(p4,6);\
id[0]--;id[1]-=3;id[2]++;       FILL_NORM(p3,4);\
id[0]--;id[1]++;        FILL_NORM(p1,2);\
id[0]++;        FILL_NORM(p3,0);\
id[0]++;        FILL_NORM(p3,6);\
id[0]-=2;id[1]++;       FILL_NORM(p2,6);\
id[0]++;        FILL_NORM(p2,0);\
id[0]++;        FILL_NORM(p4,2);\
id[0]--;id[1]++;        FILL_NORM(p2,4);\
id[1]-=2;id[2]++;       FILL_NORM(p3,3);\
id[0]++;        FILL_NORM(p3,7);\
id[0]--;id[1]++;        FILL_NORM(p2,1);\
id[0]++;        FILL_NORM(p3,1);\
id[0]--;id[1]++;        FILL_NORM(p2,7);\
id[0]++;        FILL_NORM(p2,3);\
id[0]--;id[1]-=2;id[2]++;       FILL_NORM(p3,5);\
id[1]++;        FILL_NORM(p2,5);\

#define Norm_RegionFour_FillData_Odd \
id[2]-=2;       FILL_NORM(p1,5);\
id[0]++;        FILL_NORM(p4,5);\
id[0]--;id[2]++;        FILL_NORM(p1,7);\
id[0]++;        FILL_NORM(p1,1);\
id[0]++;        FILL_NORM(p4,3);\
id[0]-=2;id[1]++;       FILL_NORM(p1,3);\
id[0]++;        FILL_NORM(p4,1);\
id[0]++;        FILL_NORM(p4,7);\
id[0]-=2;id[1]-=2;id[2]++;      FILL_NORM(p1,6);\
id[0]++;        FILL_NORM(p3,2);\
id[0]-=2;id[1]++;       FILL_NORM(p1,4);\
id[0]++;        FILL_NORM(p1,0);\
id[0]++;        FILL_NORM(p4,0);\
id[0]++;        FILL_NORM(p4,4);\
id[0]-=2;id[1]++;       FILL_NORM(p2,2);\
id[0]++;        FILL_NORM(p4,6);\
id[1]-=2;id[2]++;       FILL_NORM(p3,4);\
id[0]--;id[1]++;        FILL_NORM(p1,2);\
id[0]++;        FILL_NORM(p3,0);\
id[0]++;        FILL_NORM(p3,6);\
id[0]-=2;id[1]++;       FILL_NORM(p2,6);\
id[0]++;        FILL_NORM(p2,0);\
id[0]++;        FILL_NORM(p4,2);\
id[0]--;id[1]++;        FILL_NORM(p2,4);\
id[0]--;id[1]-=3;id[2]++;       FILL_NORM(p3,3);\
id[0]++;        FILL_NORM(p3,7);\
id[0]--;id[1]++;        FILL_NORM(p2,1);\
id[0]++;        FILL_NORM(p3,1);\
id[0]--;id[1]++;        FILL_NORM(p2,7);\
id[0]++;        FILL_NORM(p2,3);\
id[1]--;id[2]++;        FILL_NORM(p3,5);\
id[1]++;        FILL_NORM(p2,5);\

/////////////////////////////////


/////////////////////////////////

#define Norm_RegionFive_FillData_Even \
id[2]-=2;       FILL_NORM(p1,4);\
id[0]--;id[1]--;id[2]++;        FILL_NORM(p1,7);\
id[0]++;        FILL_NORM(p1,2);\
id[0]--;id[1]++;        FILL_NORM(p1,3);\
id[0]++;        FILL_NORM(p2,6);\
id[1]--;id[2]++;        FILL_NORM(p1,6);\
id[0]++;        FILL_NORM(p3,3);\
id[0]-=2;id[1]++;       FILL_NORM(p1,5);\
id[0]++;        FILL_NORM(p1,0);\
id[0]++;        FILL_NORM(p2,1);\
id[0]--;id[1]++;        FILL_NORM(p2,2);\
id[0]++;        FILL_NORM(p2,7);\
id[0]--;id[1]-=3;id[2]++;       FILL_NORM(p3,4);\
id[0]--;id[1]++;        FILL_NORM(p1,1);\
id[0]++;        FILL_NORM(p3,0);\
id[0]++;        FILL_NORM(p3,5);\
id[0]-=2;id[1]++;       FILL_NORM(p4,1);\
id[0]++;        FILL_NORM(p2,0);\
id[0]++;        FILL_NORM(p2,5);\
id[0]--;id[1]++;        FILL_NORM(p2,4);\
id[1]-=2;id[2]++;       FILL_NORM(p3,2);\
id[0]++;        FILL_NORM(p3,7);\
id[0]-=2;id[1]++;       FILL_NORM(p4,5);\
id[0]++;        FILL_NORM(p4,0);\
id[0]++;        FILL_NORM(p3,1);\
id[0]--;id[1]++;        FILL_NORM(p4,6);\
id[0]++;        FILL_NORM(p2,3);\
id[0]-=2;id[1]-=2;id[2]++;      FILL_NORM(p4,3);\
id[0]++;        FILL_NORM(p3,6);\
id[0]--;id[1]++;        FILL_NORM(p4,7);\
id[0]++;        FILL_NORM(p4,2);\
id[2]++;        FILL_NORM(p4,4);\

#define Norm_RegionFive_FillData_Odd \
id[2]-=2;       FILL_NORM(p1,4);\
id[2]++;        FILL_NORM(p1,7);\
id[0]++;        FILL_NORM(p1,2);\
id[0]--;id[1]++;        FILL_NORM(p1,3);\
id[0]++;        FILL_NORM(p2,6);\
id[0]--;id[1]-=2;id[2]++;       FILL_NORM(p1,6);\
id[0]++;        FILL_NORM(p3,3);\
id[0]-=2;id[1]++;       FILL_NORM(p1,5);\
id[0]++;        FILL_NORM(p1,0);\
id[0]++;        FILL_NORM(p2,1);\
id[0]--;id[1]++;        FILL_NORM(p2,2);\
id[0]++;        FILL_NORM(p2,7);\
id[1]-=2;id[2]++;       FILL_NORM(p3,4);\
id[0]--;id[1]++;        FILL_NORM(p1,1);\
id[0]++;        FILL_NORM(p3,0);\
id[0]++;        FILL_NORM(p3,5);\
id[0]-=2;id[1]++;       FILL_NORM(p4,1);\
id[0]++;        FILL_NORM(p2,0);\
id[0]++;        FILL_NORM(p2,5);\
id[0]--;id[1]++;        FILL_NORM(p2,4);\
id[0]--;id[1]-=3;id[2]++;       FILL_NORM(p3,2);\
id[0]++;        FILL_NORM(p3,7);\
id[0]-=2;id[1]++;       FILL_NORM(p4,5);\
id[0]++;        FILL_NORM(p4,0);\
id[0]++;        FILL_NORM(p3,1);\
id[0]--;id[1]++;        FILL_NORM(p4,6);\
id[0]++;        FILL_NORM(p2,3);\
id[0]--;id[1]--;id[2]++;        FILL_NORM(p4,3);\
id[0]++;        FILL_NORM(p3,6);\
id[0]--;id[1]++;        FILL_NORM(p4,7);\
id[0]++;        FILL_NORM(p4,2);\
id[0]--;id[1]--;id[2]++;        FILL_NORM(p4,4);\

/////////////////////////////////



/////////////////////////////////

#define Norm_RegionSix_FillData_Even \
id[2]-=2;       FILL_NORM(p1,4);\
id[0]--;id[1]--;id[2]++;        FILL_NORM(p1,7);\
id[0]++;        FILL_NORM(p1,3);\
id[0]--;id[1]++;        FILL_NORM(p1,2);\
id[0]++;        FILL_NORM(p2,6);\
id[1]--;id[2]++;        FILL_NORM(p1,5);\
id[0]--;id[1]++;        FILL_NORM(p1,6);\
id[0]++;        FILL_NORM(p1,0);\
id[0]++;        FILL_NORM(p2,2);\
id[0]-=2;id[1]++;       FILL_NORM(p3,3);\
id[0]++;        FILL_NORM(p2,1);\
id[0]++;        FILL_NORM(p2,7);\
id[0]-=2;id[1]-=2;id[2]++;      FILL_NORM(p1,1);\
id[0]++;        FILL_NORM(p4,1);\
id[0]-=2;id[1]++;       FILL_NORM(p3,4);\
id[0]++;        FILL_NORM(p3,0);\
id[0]++;        FILL_NORM(p2,0);\
id[0]++;        FILL_NORM(p2,4);\
id[0]-=2;id[1]++;       FILL_NORM(p3,5);\
id[0]++;        FILL_NORM(p2,5);\
id[1]-=2;id[2]++;       FILL_NORM(p4,5);\
id[0]--;id[1]++;        FILL_NORM(p3,2);\
id[0]++;        FILL_NORM(p4,0);\
id[0]++;        FILL_NORM(p4,6);\
id[0]-=2;id[1]++;       FILL_NORM(p3,7);\
id[0]++;        FILL_NORM(p3,1);\
id[0]++;        FILL_NORM(p2,3);\
id[0]-=2;id[1]-=2;id[2]++;      FILL_NORM(p4,3);\
id[0]++;        FILL_NORM(p4,7);\
id[0]--;id[1]++;        FILL_NORM(p3,6);\
id[0]++;        FILL_NORM(p4,2);\
id[2]++;        FILL_NORM(p4,4);\

#define Norm_RegionSix_FillData_Odd \
id[2]-=2;       FILL_NORM(p1,4);\
id[2]++;        FILL_NORM(p1,7);\
id[0]++;        FILL_NORM(p1,3);\
id[0]--;id[1]++;        FILL_NORM(p1,2);\
id[0]++;        FILL_NORM(p2,6);\
id[0]--;id[1]-=2;id[2]++;       FILL_NORM(p1,5);\
id[0]--;id[1]++;        FILL_NORM(p1,6);\
id[0]++;        FILL_NORM(p1,0);\
id[0]++;        FILL_NORM(p2,2);\
id[0]-=2;id[1]++;       FILL_NORM(p3,3);\
id[0]++;        FILL_NORM(p2,1);\
id[0]++;        FILL_NORM(p2,7);\
id[0]--;id[1]--;id[2]++;        FILL_NORM(p1,1);\
id[0]++;        FILL_NORM(p4,1);\
id[0]-=2;id[1]++;       FILL_NORM(p3,4);\
id[0]++;        FILL_NORM(p3,0);\
id[0]++;        FILL_NORM(p2,0);\
id[0]++;        FILL_NORM(p2,4);\
id[0]-=2;id[1]++;       FILL_NORM(p3,5);\
id[0]++;        FILL_NORM(p2,5);\
id[0]--;id[1]-=3;id[2]++;       FILL_NORM(p4,5);\
id[0]--;id[1]++;        FILL_NORM(p3,2);\
id[0]++;        FILL_NORM(p4,0);\
id[0]++;        FILL_NORM(p4,6);\
id[0]-=2;id[1]++;       FILL_NORM(p3,7);\
id[0]++;        FILL_NORM(p3,1);\
id[0]++;        FILL_NORM(p2,3);\
id[0]--;id[1]--;id[2]++;        FILL_NORM(p4,3);\
id[0]++;        FILL_NORM(p4,7);\
id[0]--;id[1]++;        FILL_NORM(p3,6);\
id[0]++;        FILL_NORM(p4,2);\
id[0]--;id[1]--;id[2]++;        FILL_NORM(p4,4);\


#define RegionOne_FillData_Even \
id[2]-=3;       p[20]=Fetch( id[0], id[1], id[2] );\
id[2]++;        p[6]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[18]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[19]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[23]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]-=2;id[2]++;      p[7]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[1]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[27]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[2]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[16]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[22]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[21]=Fetch( id[0], id[1], id[2] );\
id[1]-=2;id[2]++;       p[5]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[29]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[4]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[0]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[24]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[28]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[9]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[17]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]-=2;id[2]++;      p[3]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[25]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[31]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[14]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[8]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[26]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[13]=Fetch( id[0], id[1], id[2] );\
id[1]--;id[2]++;        p[10]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[30]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[15]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[11]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]--;id[2]++;        p[12]=Fetch( id[0], id[1], id[2] );\

#define RegionOne_FillData_Odd \
id[0]++;id[1]++;id[2]-=3;       p[20]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]--;id[2]++;        p[6]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[18]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[19]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[23]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]--;id[2]++;        p[7]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[1]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[27]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[2]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[16]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[22]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[21]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]-=3;id[2]++;       p[5]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[29]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[4]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[0]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[24]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[28]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[9]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[17]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]--;id[2]++;        p[3]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[25]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[31]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[14]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[8]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[26]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[13]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]-=2;id[2]++;       p[10]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[30]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[15]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[11]=Fetch( id[0], id[1], id[2] );\
id[2]++;        p[12]=Fetch( id[0], id[1], id[2] );\

/////////////////////////////////


/////////////////////////////////

#define RegionTwo_FillData_Even \
id[2]-=3;       p[20]=Fetch( id[0], id[1], id[2] );\
id[2]++;        p[6]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[19]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[18]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[23]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]-=2;id[2]++;      p[7]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[2]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[1]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[16]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[21]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[27]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[22]=Fetch( id[0], id[1], id[2] );\
id[1]-=2;id[2]++;       p[4]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[5]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[0]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[9]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[29]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[24]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[17]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[28]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]-=3;id[2]++;       p[3]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[14]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[25]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[8]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[13]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[31]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[26]=Fetch( id[0], id[1], id[2] );\
id[1]--;id[2]++;        p[10]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[15]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[30]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[11]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]--;id[2]++;        p[12]=Fetch( id[0], id[1], id[2] );\

#define RegionTwo_FillData_Odd \
id[0]++;id[1]++;id[2]-=3;       p[20]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]--;id[2]++;        p[6]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[19]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[18]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[23]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]--;id[2]++;        p[7]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[2]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[1]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[16]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[21]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[27]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[22]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]-=3;id[2]++;       p[4]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[5]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[0]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[9]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[29]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[24]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[17]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[28]=Fetch( id[0], id[1], id[2] );\
id[1]-=2;id[2]++;       p[3]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[14]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[25]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[8]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[13]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[31]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[26]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]-=2;id[2]++;       p[10]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[15]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[30]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[11]=Fetch( id[0], id[1], id[2] );\
id[2]++;        p[12]=Fetch( id[0], id[1], id[2] );\

/////////////////////////////////


/////////////////////////////////

#define RegionThree_FillData_Even \
id[2]-=2;       p[5]=Fetch( id[0], id[1], id[2] );\
id[1]++;        p[29]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]-=2;id[2]++;       p[7]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[3]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[1]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[25]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[27]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[31]=Fetch( id[0], id[1], id[2] );\
id[1]-=2;id[2]++;       p[4]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[6]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[0]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[10]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[18]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[24]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[30]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[28]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]-=3;id[2]++;       p[2]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[14]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[20]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[16]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[8]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[12]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[22]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[26]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]--;id[2]++;        p[19]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[9]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[15]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[23]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[17]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[11]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]--;id[2]++;       p[21]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[13]=Fetch( id[0], id[1], id[2] );\

#define RegionThree_FillData_Odd \
id[2]-=2;       p[5]=Fetch( id[0], id[1], id[2] );\
id[1]++;        p[29]=Fetch( id[0], id[1], id[2] );\
id[1]--;id[2]++;        p[7]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[3]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[1]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[25]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[27]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[31]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]-=3;id[2]++;       p[4]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[6]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[0]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[10]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[18]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[24]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[30]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[28]=Fetch( id[0], id[1], id[2] );\
id[1]-=2;id[2]++;       p[2]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[14]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[20]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[16]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[8]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[12]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[22]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[26]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]-=2;id[2]++;      p[19]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[9]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[15]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[23]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[17]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[11]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[2]++;        p[21]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[13]=Fetch( id[0], id[1], id[2] );\

/////////////////////////////////


/////////////////////////////////

#define RegionFour_FillData_Even \
id[2]-=2;       p[5]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[29]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]--;id[2]++;       p[7]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[1]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[27]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[3]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[25]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[31]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]--;id[2]++;        p[6]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[18]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[4]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[0]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[24]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[28]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[10]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[30]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]-=3;id[2]++;       p[20]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[2]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[16]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[22]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[14]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[8]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[26]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[12]=Fetch( id[0], id[1], id[2] );\
id[1]-=2;id[2]++;       p[19]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[23]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[9]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[17]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[15]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[11]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]-=2;id[2]++;       p[21]=Fetch( id[0], id[1], id[2] );\
id[1]++;        p[13]=Fetch( id[0], id[1], id[2] );\

#define RegionFour_FillData_Odd \
id[2]-=2;       p[5]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[29]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[2]++;        p[7]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[1]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[27]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[3]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[25]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[31]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]-=2;id[2]++;      p[6]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[18]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[4]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[0]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[24]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[28]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[10]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[30]=Fetch( id[0], id[1], id[2] );\
id[1]-=2;id[2]++;       p[20]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[2]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[16]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[22]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[14]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[8]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[26]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[12]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]-=3;id[2]++;       p[19]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[23]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[9]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[17]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[15]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[11]=Fetch( id[0], id[1], id[2] );\
id[1]--;id[2]++;        p[21]=Fetch( id[0], id[1], id[2] );\
id[1]++;        p[13]=Fetch( id[0], id[1], id[2] );\

/////////////////////////////////


/////////////////////////////////

#define RegionFive_FillData_Even \
id[2]-=2;       p[4]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]--;id[2]++;        p[7]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[2]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[3]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[14]=Fetch( id[0], id[1], id[2] );\
id[1]--;id[2]++;        p[6]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[19]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[5]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[0]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[9]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[10]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[15]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]-=3;id[2]++;       p[20]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[1]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[16]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[21]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[25]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[8]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[13]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[12]=Fetch( id[0], id[1], id[2] );\
id[1]-=2;id[2]++;       p[18]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[23]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[29]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[24]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[17]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[30]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[11]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]-=2;id[2]++;      p[27]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[22]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[31]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[26]=Fetch( id[0], id[1], id[2] );\
id[2]++;        p[28]=Fetch( id[0], id[1], id[2] );\

#define RegionFive_FillData_Odd \
id[2]-=2;       p[4]=Fetch( id[0], id[1], id[2] );\
id[2]++;        p[7]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[2]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[3]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[14]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]-=2;id[2]++;       p[6]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[19]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[5]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[0]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[9]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[10]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[15]=Fetch( id[0], id[1], id[2] );\
id[1]-=2;id[2]++;       p[20]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[1]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[16]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[21]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[25]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[8]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[13]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[12]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]-=3;id[2]++;       p[18]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[23]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[29]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[24]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[17]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[30]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[11]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]--;id[2]++;        p[27]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[22]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[31]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[26]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]--;id[2]++;        p[28]=Fetch( id[0], id[1], id[2] );\

/////////////////////////////////



/////////////////////////////////

#define RegionSix_FillData_Even \
id[2]-=2;       p[4]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]--;id[2]++;        p[7]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[3]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[2]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[14]=Fetch( id[0], id[1], id[2] );\
id[1]--;id[2]++;        p[5]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[6]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[0]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[10]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[19]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[9]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[15]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]-=2;id[2]++;      p[1]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[25]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[20]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[16]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[8]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[12]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[21]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[13]=Fetch( id[0], id[1], id[2] );\
id[1]-=2;id[2]++;       p[29]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[18]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[24]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[30]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[23]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[17]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[11]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]-=2;id[2]++;      p[27]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[31]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[22]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[26]=Fetch( id[0], id[1], id[2] );\
id[2]++;        p[28]=Fetch( id[0], id[1], id[2] );\

#define RegionSix_FillData_Odd \
id[2]-=2;       p[4]=Fetch( id[0], id[1], id[2] );\
id[2]++;        p[7]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[3]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[2]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[14]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]-=2;id[2]++;       p[5]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[6]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[0]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[10]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[19]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[9]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[15]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]--;id[2]++;        p[1]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[25]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[20]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[16]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[8]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[12]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[21]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[13]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]-=3;id[2]++;       p[29]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[18]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[24]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[30]=Fetch( id[0], id[1], id[2] );\
id[0]-=2;id[1]++;       p[23]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[17]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[11]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]--;id[2]++;        p[27]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[31]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]++;        p[22]=Fetch( id[0], id[1], id[2] );\
id[0]++;        p[26]=Fetch( id[0], id[1], id[2] );\
id[0]--;id[1]--;id[2]++;        p[28]=Fetch( id[0], id[1], id[2] );\

#include "conv_all_macros.h"



bool CBCCData::m_bMapInitialized = false;
std::map< std::string, CBCCData::NORMAL_DISCRETE_FILTER> CBCCData::m_normalFilterNameMap;
std::map< std::string, CBCCData::RECONSTRUCTION_FILTER>  CBCCData::m_reconstructionFilterNameMap;




void CBCCData::InitializeMaps()
{
	if( m_bMapInitialized ) return;

	MAP_ENTRY_BY_FLAG( m_reconstructionFilterNameMap, BCC_QUINTIC_BOX_SPLINE );
	MAP_ENTRY_BY_FLAG( m_reconstructionFilterNameMap, BCC_LINEAR_BOX_SPLINE );

	MAP_ENTRY_BY_FLAG( m_normalFilterNameMap, BCC_NORMAL_BCD );
	MAP_ENTRY_BY_FLAG( m_normalFilterNameMap, BCC_NORMAL_SOCD );
	MAP_ENTRY_BY_FLAG( m_normalFilterNameMap, BCC_NORMAL_ICD );
	MAP_ENTRY_BY_FLAG( m_normalFilterNameMap, BCC_NORMAL_SIZE_OPTIMAL16 );
	MAP_ENTRY_BY_FLAG( m_normalFilterNameMap, BCC_NORMAL_SEMI_OPTIMAL18 );
	MAP_ENTRY_BY_FLAG( m_normalFilterNameMap, BCC_NORMAL_OPTIMAL26 );

	m_bMapInitialized = true;
}

CBCCData::CBCCData( const char* name, const char *fname, double xScale, double yScale, double zScale, NORMAL_DISCRETE_FILTER normalEst, bool bOnTheFlyNormals, const char *gridNormalOutputFile ) : CBaseData(std::string(name)),m_normalEst(normalEst)
{
	m_reconFilter = BCC_QUINTIC_BOX_SPLINE;
	m_normalReconFilter = BCC_QUINTIC_BOX_SPLINE;
	m_bOnTheFlyNormals = bOnTheFlyNormals;

	if( gridNormalOutputFile){
		strcpy( m_gridNormalOutputFile, gridNormalOutputFile );
	}
	else
		m_gridNormalOutputFile[0] = 0;

	m_halfScale[X] = xScale / 2.0f;
	m_halfScale[Y] = yScale / 2.0f;
	m_halfScale[Z] = zScale / 1.0f;

	m_invHalfScale[X] = 1.0/m_halfScale[X];
	m_invHalfScale[Y] = 1.0/m_halfScale[Y];
	m_invHalfScale[Z] = 1.0/m_halfScale[Z];
			
	

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
		s = fgets(line, 256, fp ); // reads the first commented line
	}

	char dummy[20];
	//line should have dimension information now
	sscanf(line,"%s %d %d %d", dummy, &xSize, &ySize, &zSize );

	//now just read off the next 7 lines
	for(int i=0;i<7;++i){
		s = fgets(line, 256, fp ); // reads the first commented line
	}



	printf("Reading %s (%dx%dx%d).... ", fname, xSize, ySize, zSize );

	int totalAtomic = xSize * ySize * zSize;
	atomic_t *pMem = new atomic_t [ totalAtomic ];

	readsize = fread(pMem, totalAtomic * sizeof(atomic_t),1,fp);

	fclose(fp);

	CreateData( xSize, ySize, zSize, xScale, yScale, zScale, pMem );

	//Nasty hack to get the BCC grid right :(
	DELETE_SAFE(m_pVolume);
	m_pVolume = new CVolume( (m_xSize - 0.5) * xScale, (m_ySize - 0.5) * yScale, (m_zSize - 1) * zScale, xScale, yScale, zScale );

	DELETE_SAFE_ARRAY( pMem );


	m_translate[X] = m_pVolume->m_HXSize / m_halfScale[X];
	m_translate[Y] = m_pVolume->m_HYSize / m_halfScale[Y];
	m_translate[Z] = m_pVolume->m_HZSize / m_halfScale[Z];
	
	printf("Done!\n");
}

CBCCData::CBCCData( const char* name,
				    const char* fname,
					DATA_FORMAT fmt,
					double xScale,
					double yScale,
					double zScale,
					const char* filter,
					bool  bOnTheFlyNormals,
					const char* gridNormalOutputFile ) : CBaseData(std::string(name))
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

	m_halfScale[X] = xScale / 2.0f;
	m_halfScale[Y] = yScale / 2.0f;
	m_halfScale[Z] = zScale / 1.0f;

	m_invHalfScale[X] = 1.0/m_halfScale[X];
	m_invHalfScale[Y] = 1.0/m_halfScale[Y];
	m_invHalfScale[Z] = 1.0/m_halfScale[Z];

	m_xSize = m_ySize = m_zSize = 0;
	m_pppData = 0;

	CreateDataFromFile( fname,fmt,xScale,yScale,zScale );

	//Nasty hack to get the BCC grid right :(
	DELETE_SAFE(m_pVolume);
	m_pVolume = new CVolume( (m_xSize - 0.5) * xScale, (m_ySize - 0.5) * yScale, (m_zSize - 1) * zScale, xScale, yScale, zScale );

	m_translate[X] = m_pVolume->m_HXSize / m_halfScale[X];
	m_translate[Y] = m_pVolume->m_HYSize / m_halfScale[Y];
	m_translate[Z] = m_pVolume->m_HZSize / m_halfScale[Z];
}

/* This constructor assumes that we have float data tyep in the VUD file */


/* This constructor assumes that we have float data tyep in the VUD file */

CBCCData::CBCCData( const char* name,
				    const char* fname,
					DATA_FORMAT fmt,
					const char *xname,
					const char *yname,
					const char *zname,
					double xScale,
					double yScale,
					double zScale,
					const char* gridNormalOutputFile) : CBaseData( std::string( name ) )
{

	m_reconFilter = BCC_QUINTIC_BOX_SPLINE;
	m_normalReconFilter = BCC_QUINTIC_BOX_SPLINE;
	m_bOnTheFlyNormals = false;


	if( gridNormalOutputFile){
			strcpy( m_gridNormalOutputFile, gridNormalOutputFile );
	}
	else
		m_gridNormalOutputFile[0] = 0;

	int xsize,ysize,zsize;

	m_halfScale[X] = xScale / 2.0f;
	m_halfScale[Y] = yScale / 2.0f;
	m_halfScale[Z] = zScale / 1.0f;

	m_invHalfScale[X] = 1.0/m_halfScale[X];
	m_invHalfScale[Y] = 1.0/m_halfScale[Y];
	m_invHalfScale[Z] = 1.0/m_halfScale[Z];

	switch( fmt ){
			case DATA_FORMAT_BYTE:
				m_pppData       = LoadVUDFileByte( fname, m_xSize, m_ySize, m_zSize );
				m_sepNormals[X] = LoadVUDFileByte( xname, xsize,ysize,zsize );
				m_sepNormals[Y] = LoadVUDFileByte( yname, xsize,ysize,zsize  );
				m_sepNormals[Z] = LoadVUDFileByte( zname, xsize,ysize,zsize  );
				break;
			case DATA_FORMAT_FLOAT:
				m_pppData       = LoadVUDFileFloat( fname, m_xSize, m_ySize, m_zSize );
				m_sepNormals[X] = LoadVUDFileFloat( xname, xsize,ysize,zsize );
				m_sepNormals[Y] = LoadVUDFileFloat( yname, xsize,ysize,zsize  );
				m_sepNormals[Z] = LoadVUDFileFloat( zname, xsize,ysize,zsize  );
				break;
		}

	// Nasty hack :(
	m_pVolume = new CVolume( (m_xSize - 0.5) * xScale, (m_ySize - 0.5) * yScale, (m_zSize - 1) * zScale, xScale, yScale, zScale );

	m_translate[X] = m_pVolume->m_HXSize / m_halfScale[X];
	m_translate[Y] = m_pVolume->m_HYSize / m_halfScale[Y];
	m_translate[Z] = m_pVolume->m_HZSize / m_halfScale[Z];
	
	m_pppNormals = new vector3_t** [ m_zSize + 2*GUARD_BAND ];

	for(int z = 0; z< (m_zSize+2*GUARD_BAND);++z){

		m_pppNormals[z] = new vector3_t* [ m_ySize +  2*GUARD_BAND];

		for(int y=0;y<(m_ySize+2*GUARD_BAND);++y){
			m_pppNormals[z][y] = new vector3_t [ m_xSize + 2*GUARD_BAND ];

            memset(m_pppNormals[z][y], 0, sizeof( vector3_t) * (m_xSize +2*GUARD_BAND));

		}
	}
	
	for(int z = 0; z< m_zSize;++z){
		for(int y=0;y<m_ySize;++y){
			for(int x=0;x<m_xSize;++x){
				m_pppNormals[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND][X] = m_sepNormals[X][z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND];								
				m_pppNormals[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND][Y] = m_sepNormals[Y][z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND];
				m_pppNormals[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND][Z] = m_sepNormals[Z][z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND];
			}
		}
	}
	
	
}

bool   CBCCData::SetNormalEstimatorFilterByName( const char *filter )
{
	if( CBaseData::SetNormalEstimatorFilterByName( filter ) ){
		// This means we have filter in the base class with this name
		return true;
	}

	bool ret = true;

	InitializeMaps();

	if( m_normalFilterNameMap.find( std::string(filter) ) == m_normalFilterNameMap.end() ){
		printf(" BCC Lattice ( %s )\n", GetName() );
		printf("\tThe normal discrete filter %s doesn't exists\n", filter);
		printf("\tNormal discrete filter will be set to default\n");

		ret = false;
	}

	m_normalEst = m_normalFilterNameMap[std::string(filter)];

	return ret;
}

bool   CBCCData::SetReconstructionFilterByName( const char *filter )
{
	if( CBaseData::SetReconstructionFilterByName( filter ) ){
		return true;
	}

	bool ret = true;

	InitializeMaps();

	if( m_reconstructionFilterNameMap.find( std::string(filter) ) == m_reconstructionFilterNameMap.end() ){
		printf(" BCC Lattice ( %s )\n", GetName() );
		printf("\tThe reconstruction filter %s doesn't exists\n", filter);
		printf("\tReconstruction filter will be set to default\n");

		ret = false;
	}

	m_reconFilter = m_reconstructionFilterNameMap[std::string(filter)];

	return ret;

}

bool CBCCData::SetNormalReconstructionFilterByName( const char* filter )
{

	if( CBaseData::SetNormalReconstructionFilterByName( filter )) {
		return true;
	}

	bool ret = true;

	InitializeMaps();

	if( m_reconstructionFilterNameMap.find( std::string(filter) ) == m_reconstructionFilterNameMap.end() ){
		printf(" BCC Lattice ( %s )\n", GetName() );
		printf("\tThe normal reconstruction filter %s doesn't exists\n", filter);
		printf("\tNormal reconstruction filter will be set to default\n");

		ret = false;
	}

	m_normalReconFilter = m_reconstructionFilterNameMap[std::string(filter)];

	return ret;
}

void CBCCData::HookGradientComponentEstimators( )
{
	if( m_bEpsilonNormalMode ) return;
	
	if( !m_bOnTheFlyNormals ){
		printf("Setting up normal estimators for precomputed normals\n");
 		m_gradEstimator[0] = new CBCCGradientXEstimatorPrecomp( this );
 		m_gradEstimator[1] = new CBCCGradientYEstimatorPrecomp( this );
 		m_gradEstimator[2] = new CBCCGradientZEstimatorPrecomp( this );
	}
	else{
		printf("Setting up normal estimators for on-the-fly normals\n");

		switch( m_normalEst ){
			case BCC_NORMAL_SOCD:
				m_gradEstimator[0] = new BCCESTIMATORCLASSNAME(X,CentralDifference)( this );
				m_gradEstimator[1] = new BCCESTIMATORCLASSNAME(Y,CentralDifference)( this );
				m_gradEstimator[2] = new BCCESTIMATORCLASSNAME(Z,CentralDifference)( this );
				break;
			case BCC_NORMAL_BCD:
				m_gradEstimator[0] = new BCCESTIMATORCLASSNAME(X,BoxDifference)( this );
				m_gradEstimator[1] = new BCCESTIMATORCLASSNAME(Y,BoxDifference)( this );
				m_gradEstimator[2] = new BCCESTIMATORCLASSNAME(Z,BoxDifference)( this );
				break;
			case BCC_NORMAL_SIZE_OPTIMAL16:
				m_gradEstimator[0] = new BCCESTIMATORCLASSNAME(X,SizeOptimalFilter)( this );
				m_gradEstimator[1] = new BCCESTIMATORCLASSNAME(Y,SizeOptimalFilter)( this );
				m_gradEstimator[2] = new BCCESTIMATORCLASSNAME(Z,SizeOptimalFilter)( this );
				break;
			case BCC_NORMAL_OPTIMAL26:
				m_gradEstimator[0] = new BCCESTIMATORCLASSNAME(X,OptimalFilter)( this );
				m_gradEstimator[1] = new BCCESTIMATORCLASSNAME(Y,OptimalFilter)( this );
				m_gradEstimator[2] = new BCCESTIMATORCLASSNAME(Z,OptimalFilter)( this );
				break;
			default:
				printf("The normal estimation filter is not supported for On-the-fly computation!\n");
				break;

		}
	}
}

void CBCCData::PrecompNormalInterpolatedCentralDifference(vector3_t N, int x, int y, int z)
{
	int index3D[3] = {x,y,z};
	int bccX,bccY,bccZ;

	IndexToBCC( bccX,bccY,bccZ, index3D );
	double q[6];

	vector3_t bcc;

	bcc[X] = bccX + 1;
	bcc[Y] = bccY;
	bcc[Z] = bccZ;

	q[0] = Interpolate( bcc );

	bcc[X] = bccX - 1;
	bcc[Y] = bccY;
	bcc[Z] = bccZ;

	q[1] = Interpolate( bcc );

	bcc[X] = bccX;
	bcc[Y] = bccY + 1;
	bcc[Z] = bccZ;

	q[2] = Interpolate( bcc );

	bcc[X] = bccX;
	bcc[Y] = bccY - 1;
	bcc[Z] = bccZ;

	q[3] = Interpolate( bcc );

	bcc[X] = bccX;
	bcc[Y] = bccY;
	bcc[Z] = bccZ+1;

	q[4] = Interpolate( bcc );

	bcc[X] = bccX;
	bcc[Y] = bccY;
	bcc[Z] = bccZ-1;

	q[5] = Interpolate( bcc );

	N[X] = (q[0] - q[1]) / (2.0f * m_halfScale[X]);
	N[Y] = (q[2] - q[3]) / (2.0f * m_halfScale[Y]);
	N[Z] = (q[4] - q[5]) / (2.0f * m_halfScale[Z]);

}

elem_t  CBCCData::ComputeNormalXBoxDifference( int x,int y, int z) const
{

	int zp1 = z+1;
	int zm1 = z-1;

	if( z & 0x01 ){ //ODD
		int xp1 = x+1;
		int yp1 = y+1;

		return 0.125 * m_invHalfScale[X] * ( (Fetch( xp1,y,zp1 ) + Fetch( xp1, yp1, zp1 ) +
				Fetch( xp1,y,zm1) + Fetch( xp1,yp1,zm1))
				-
			   (Fetch(x,y,zp1) + Fetch(x,yp1,zp1) +
				Fetch(x,y,zm1) + Fetch(x,yp1,zm1) )
			);
	}
	else{
		int xm1 = x-1;
		int ym1 = y-1;

		return 0.125 * m_invHalfScale[X] * ((Fetch(x,ym1,zp1) + Fetch(x,y,zp1) + Fetch( x,ym1,zm1) + Fetch(x,y,zm1 ))
				-
				(Fetch(xm1,ym1,zp1) + Fetch(xm1,y,zp1) + Fetch(xm1,ym1,zm1) + Fetch(xm1,y,zm1)));
	}

	return 0;
}

elem_t  CBCCData::ComputeNormalYBoxDifference( int x,int y, int z) const
{

	int zp1 = z+1;
	int zm1 = z-1;

	if( z & 0x01 ){ //ODD
		int xp1 = x+1;
		int yp1 = y+1;

		return 0.125 * m_invHalfScale[Y] * ((Fetch(xp1,yp1,zp1) + Fetch(x,yp1,zp1) + Fetch(xp1,yp1,zm1) + Fetch(x,yp1,zm1))
				-
				(Fetch(xp1,y,zp1) + Fetch(x,y,zp1) + Fetch(xp1,y,zm1) + Fetch(x,y,zm1)));
	}
	else{
		int xm1 = x-1;
		int ym1 = y-1;

		return 0.125 * m_invHalfScale[Y] * ((Fetch(x,y,zp1) + Fetch(xm1,y,zp1) + Fetch(x,y,zm1) + Fetch(xm1,y,zm1))
				-
				(Fetch(x,ym1,zp1) + Fetch(xm1,ym1,zp1) + Fetch(x,ym1,zm1) + Fetch(xm1,ym1,zm1)));

	}


	return 0;
}

elem_t  CBCCData::ComputeNormalZBoxDifference( int x,int y, int z) const
{
	int zp1 = z+1;
	int zm1 = z-1;

	if( z & 0x01 ){ //ODD
		int xp1 = x+1;
		int yp1 = y+1;

		return 0.125 * m_invHalfScale[Z] * ((Fetch(x,y,zp1) + Fetch(x,yp1,zp1) + Fetch(xp1,y,zp1) + Fetch(xp1,yp1,zp1))
				-
				(Fetch(x,y,zm1) + Fetch(x,yp1,zm1) + Fetch(xp1,y,zm1) + Fetch(xp1,yp1,zm1)));
	}
	else{
		int xm1 = x-1;
		int ym1 = y-1;

		return 0.125 * m_invHalfScale[Z] * ((Fetch(xm1,ym1,zp1) + Fetch(xm1,y,zp1) + Fetch(x,ym1,zp1) + Fetch(x,y,zp1))
				-
				(Fetch(xm1,ym1,zm1) + Fetch(xm1,y,zm1) + Fetch(x,ym1,zm1) + Fetch(x,y,zm1)));

	}

	return 0;
}

#define BCC_DELTA(out,x,y,z,delX,delY,delZ) \
	do{ \
		int bcc[3];				\
		int idx[3];				\
		idx[0] = x;				\
		idx[1] = y;				\
		idx[2] = z;				\
		IndexToBCC(bcc[0],bcc[1],bcc[2],idx);	\
		bcc[X] += delX;				\
		bcc[Y] += delY;				\
		bcc[Z] += delZ;				\
		BCCToIndex(idx,bcc[0],bcc[1],bcc[2]);	\
		if( (idx[0] >= 0) && (idx[0] < m_xSize) && 	\
		    (idx[1] >= 0) && (idx[1] < m_ySize) && 	\
		    (idx[2] >= 0) && (idx[2] < m_zSize) ) 	\
			(out) = Fetch(idx[0],idx[1],idx[2]);	\
	}while(0)

void CBCCData::PrecompNormalBoxDifference( vector3_t N, int x, int y, int z )
{
	elem_t p[8] = {0};

	BCC_DELTA( p[0], x,y,z, -1,-1, 1 );
	BCC_DELTA( p[1], x,y,z,  1,-1, 1 );
	BCC_DELTA( p[2], x,y,z,  1, 1, 1 );
	BCC_DELTA( p[3], x,y,z, -1, 1, 1 );
	BCC_DELTA( p[4], x,y,z, -1,-1,-1 );
	BCC_DELTA( p[5], x,y,z,  1,-1,-1 );
	BCC_DELTA( p[6], x,y,z,  1, 1,-1 );
	BCC_DELTA( p[7], x,y,z, -1, 1,-1 );

	N[X] = (p[1] - p[0] + p[2] - p[3] + p[6] - p[7] + p[5] - p[4]) / (8.0f * m_halfScale[X]);
	N[Y] = (p[3] - p[0] + p[2] - p[1] + p[7] - p[4] + p[6] - p[5]) / (8.0f * m_halfScale[Y]);
	N[Z] = (p[3] - p[7] + p[2] - p[6] + p[1] - p[5] + p[0] - p[4]) / (8.0f * m_halfScale[Z]);
}


void   CBCCData::PrecompNormalSemiOptimalFilter( vector3_t N, int x, int y, int z )
{
}

elem_t  CBCCData::ComputeNormalXOptimalFilter( int x,int y, int z) const
{
	int xm1 = x-1;
	int xp1 = x+1;
	int ym1 = y-1;
	int yp1 = y+1;
	int zm1 = z-1;
	int zp1 = z+1;
	int zm2 = z-2;
	int zp2 = z+2;

	if( z & 0x01 ) { //ODD
		return m_invHalfScale[X] * ((      0.005208333   ) *
				((Fetch(xm1,ym1,zp2)+Fetch(xm1,yp1,zp2)+Fetch(xm1,ym1,zm2)+Fetch(xm1,yp1,zm2))-
				(Fetch(xp1,ym1,zp2)+Fetch(xp1,yp1,zp2)+Fetch(xp1,ym1,zm2)+Fetch(xp1,yp1,zm2)))+
				(      0.03125   ) * ((Fetch(xm1,y,zm2)+Fetch(xm1,y,zp2)+Fetch(xm1,ym1,z)+Fetch(xm1,yp1,z))-
				(Fetch(xp1,y,zp2)+Fetch(xp1,y,zm2)+Fetch(xp1,ym1,z)+Fetch(xp1,yp1,z)))+
				(      0.0625    ) * ((Fetch(xp1,y,z))- (Fetch(xm1,y,z)))+
				(      0.166666667     ) *
				((Fetch(xp1,y,zp1)+Fetch(xp1,yp1,zp1)+Fetch(xp1,y,zm1)+Fetch(xp1,yp1,zm1))-
				(Fetch(x,y,zp1)+Fetch(x,yp1,zp1)+Fetch(x,y,zm1)+Fetch(x,yp1,zm1))));
	}
	else{
		return m_invHalfScale[X] * ((      0.005208333   ) *
				((Fetch(xm1,ym1,zp2)+Fetch(xm1,yp1,zp2)+Fetch(xm1,ym1,zm2)+Fetch(xm1,yp1,zm2))-
				(Fetch(xp1,ym1,zp2)+Fetch(xp1,yp1,zp2)+Fetch(xp1,ym1,zm2)+Fetch(xp1,yp1,zm2)))+
				(      0.03125    ) * ((Fetch(xm1,y,zm2)+Fetch(xm1,y,zp2)+Fetch(xm1,ym1,z)+Fetch(xm1,yp1,z))
				-(Fetch(xp1,y,zp2)+Fetch(xp1,y,zm2)+Fetch(xp1,ym1,z)+Fetch(xp1,yp1,z)))+
				(      0.0625    ) * ((Fetch(xp1,y,z)) - (Fetch(xm1,y,z)))+
				(      0.166666667     ) * (
				(Fetch(x,ym1,zp1)+Fetch(x,y,zp1)+Fetch(x,ym1,zm1)+Fetch(x,y,zm1)) -
				(Fetch(xm1,ym1,zp1)+Fetch(xm1,y,zp1)+Fetch(xm1,ym1,zm1)+Fetch(xm1,y,zm1))));
	}

	return 0;
}

elem_t  CBCCData::ComputeNormalYOptimalFilter( int x,int y, int z) const
{
	int xm1 = x-1;
	int xp1 = x+1;
	int ym1 = y-1;
	int yp1 = y+1;
	int zm1 = z-1;
	int zp1 = z+1;
	int zm2 = z-2;
	int zp2 = z+2;

	if( z & 0x01 ) { //ODD
		return m_invHalfScale[Y] * ((      0.005208333   ) *
				((Fetch(xp1,ym1,zp2)+Fetch(xm1,ym1,zp2)+Fetch(xp1,ym1,zm2)+Fetch(xm1,ym1,zm2))-(Fetch(xp1,yp1,
				  zp2)+Fetch(xm1,yp1,zp2)+Fetch(xp1,yp1,zm2)+Fetch(xm1,yp1,zm2)))+
				(      0.03125    ) *
				((Fetch(x,ym1,zm2)+Fetch(x,ym1,zp2)+Fetch(xp1,ym1,z)+Fetch(xm1,ym1,z))-(Fetch(x,yp1,zp2)+Fetch
				(x,yp1,zm2)+Fetch(xp1,yp1,z)+Fetch(xm1,yp1,z)))+
				(      0.0625    ) * ((Fetch(x,yp1,z))-(Fetch(x,ym1,z)))+
				(      0.166666667     ) *
				((Fetch(xp1,yp1,zp1)+Fetch(x,yp1,zp1)+Fetch(xp1,yp1,zm1)+Fetch(x,yp1,zm1))-(Fetch(xp1,y,zp1)+
				Fetch(x,y,zp1)+Fetch(xp1,y,zm1)+Fetch(x,y,zm1))));
	}
	else{
		return m_invHalfScale[Y] *((      0.005208333   ) *
				((Fetch(xp1,ym1,zp2)+Fetch(xm1,ym1,zp2)+Fetch(xp1,ym1,zm2)+Fetch(xm1,ym1,zm2))-(Fetch(xp1,yp1,
				  zp2)+Fetch(xm1,yp1,zp2)+Fetch(xp1,yp1,zm2)+Fetch(xm1,yp1,zm2)))+
				(      0.03125    ) *
				((Fetch(x,ym1,zm2)+Fetch(x,ym1,zp2)+Fetch(xp1,ym1,z)+Fetch(xm1,ym1,z))-(Fetch(x,yp1,zp2)+Fetch
				(x,yp1,zm2)+Fetch(xp1,yp1,z)+Fetch(xm1,yp1,z)))+
				(      0.0625    ) * ((Fetch(x,yp1,z))-(Fetch(x,ym1,z)))+
				(      0.166666667     ) *
				((Fetch(x,y,zp1)+Fetch(xm1,y,zp1)+Fetch(x,y,zm1)+Fetch(xm1,y,zm1))-(Fetch(x,ym1,zp1)+
				Fetch(xm1,ym1,zp1)+Fetch(x,ym1,zm1)+Fetch(xm1,ym1,zm1))));
	}

	return 0;
}

elem_t  CBCCData::ComputeNormalZOptimalFilter( int x,int y, int z) const
{
	int xm1 = x-1;
	int xp1 = x+1;
	int ym1 = y-1;
	int yp1 = y+1;
	int zm1 = z-1;
	int zp1 = z+1;
	int zm2 = z-2;
	int zp2 = z+2;

	if( z & 0x01 ) { //ODD
		return m_invHalfScale[Z] *((      0.005208333   ) *
				((Fetch(xm1,ym1,zm2)+Fetch(xm1,yp1,zm2)+Fetch(xp1,ym1,zm2)+Fetch(xp1,yp1,zm2))-(Fetch(xm1,ym1,
				  zp2)+Fetch(xm1,yp1,zp2)+Fetch(xp1,ym1,zp2)+Fetch(xp1,yp1,zp2)))+
				(      0.03125    ) *
				((Fetch(xp1,y,zm2)+Fetch(xm1,y,zm2)+Fetch(x,ym1,zm2)+Fetch(x,yp1,zm2))-(Fetch(xm1,y,zp2)+Fetch
				(xp1,y,zp2)+Fetch(x,ym1,zp2)+Fetch(x,yp1,zp2)))+
				(      0.0625    ) * ((Fetch(x,y,zp2))-(Fetch(x,y,zm2)))+
				(      0.166666667     ) * (
				(Fetch(x,y,zp1)+Fetch(x,yp1,zp1)+Fetch(xp1,y,zp1)+Fetch(xp1,yp1,zp1))-
				(Fetch(x,y,zm1)+Fetch(x,yp1,zm1)+Fetch(xp1,y,zm1)+Fetch(xp1,yp1,zm1))));
	}
	else{
		return m_invHalfScale[Z] *((      0.005208333   ) *
				((Fetch(xm1,ym1,zm2)+Fetch(xm1,yp1,zm2)+Fetch(xp1,ym1,zm2)+Fetch(xp1,yp1,zm2))-(Fetch(xm1,ym1,
				  zp2)+Fetch(xm1,yp1,zp2)+Fetch(xp1,ym1,zp2)+Fetch(xp1,yp1,zp2)))+
				(      0.03125    ) *
				((Fetch(xp1,y,zm2)+Fetch(xm1,y,zm2)+Fetch(x,ym1,zm2)+Fetch(x,yp1,zm2))-(Fetch(xm1,y,zp2)+Fetch
				(xp1,y,zp2)+Fetch(x,ym1,zp2)+Fetch(x,yp1,zp2)))+
				(      0.0625    ) * ((Fetch(x,y,zp2))-(Fetch(x,y,zm2)))+
				(      0.166666667     ) *
				((Fetch(xm1,ym1,zp1)+Fetch(xm1,y,zp1)+Fetch(x,ym1,zp1)+Fetch(x,y,zp1))-(Fetch(xm1,ym1,zm1)+
				Fetch(xm1,y,zm1)+Fetch(x,ym1,zm1)+Fetch(x,y,zm1))));
	}

	return 0;
}

void   CBCCData::PrecompNormalOptimalFilter( vector3_t N, int x, int y, int z )
{
	int domainX[26][3] = {
        //first ring
        { -1, -1,  1}, {1, -1, 1}, {1, 1, 1}, { -1, 1,  1},
        { -1, -1, -1}, {1, -1,-1}, {1, 1,-1}, { -1, 1, -1},

        //second ring
        { 2, 0, 0 }, { -2,  0,  0 },

        //third ring
		{  2,  0,  2 }, {  2, 0, -2 },
		{ -2,  0, -2 }, { -2, 0,  2 },

        { -2, -2,  0 }, {  2, -2, 0 },
        {  2,  2,  0 }, { -2,  2, 0 },

        //fourth ring
        { -2, -2,  2 }, { 2, -2, 2 }, { 2, 2, 2 }, {-2, 2, 2},
        { -2, -2, -2 }, { 2, -2,-2 }, { 2, 2,-2 }, {-2, 2,-2},
	};

	int domainY[26][3] = {
        //first ring
        { -1, -1,  1}, {1, -1, 1}, {1, 1, 1}, { -1, 1,  1},
        { -1, -1, -1}, {1, -1,-1}, {1, 1,-1}, { -1, 1, -1},

        //second ring
        { 0, 2, 0 }, { 0,  -2,  0 },

        //third ring
		{  0, -2,  2 }, { 0, -2, -2 },
		{  0,  2, -2 }, { 0,  2,  2 },

        { -2, -2,  0 }, {  2, -2, 0 },
        {  2,  2,  0 }, { -2,  2, 0 },

        //fourth ring
        { -2, -2,  2 }, { 2, -2, 2 }, { 2, 2, 2 }, {-2, 2, 2},
        { -2, -2, -2 }, { 2, -2,-2 }, { 2, 2,-2 }, {-2, 2,-2},
	};

	int domainZ[26][3] = {
        //first ring
        { -1, -1,  1}, {1, -1, 1}, {1, 1, 1}, { -1, 1,  1},
        { -1, -1, -1}, {1, -1,-1}, {1, 1,-1}, { -1, 1, -1},

        //second ring
        { 0, 0, 2 }, { 0,  0,  -2 },

        //third ring
        {  0, -2,  2 }, { 0, -2, -2 },
        {  0,  2, -2 }, { 0,  2,  2 },

        {  2,  0,  2 }, {  2, 0, -2 },
        { -2,  0, -2 }, { -2, 0,  2 },


        //fourth ring
        { -2, -2,  2 }, { 2, -2, 2 }, { 2, 2, 2 }, {-2, 2, 2},
        { -2, -2, -2 }, { 2, -2,-2 }, { 2, 2,-2 }, {-2, 2,-2},
	};

	double filterX[26] = {
	  -1/6.0,
       1/6.0,
       1/6.0,
      -1/6.0,
      -1/6.0,
       1/6.0,
       1/6.0,
      -1/6.0,
       1/16.0,
      -1/16.0,
      -1/32.0,
      -1/32.0,
       1/32.0,
       1/32.0,
       1/32.0,
      -1/32.0,
      -1/32.0,
       1/32.0,
       1/192.0,
      -1/192.0,
      -1/192.0,
       1/192.0,
       1/192.0,
      -1/192.0,
      -1/192.0,
       1/192.0
	};

	double filterY[26] = {
		-1/6.0,
		-1/6.0,
		 1/6.0,
		 1/6.0,
		-1/6.0,
		-1/6.0,
		 1/6.0,
		 1/6.0,

		 1/16.0,
		-1/16.0,

		 1/32.0,
		 1/32.0,
		-1/32.0,
		-1/32.0,
		 1/32.0,
		 1/32.0,
		-1/32.0,
		-1/32.0,

		 1/192.0,
		 1/192.0,
		-1/192.0,
		-1/192.0,
		 1/192.0,
		 1/192.0,
		-1/192.0,
		-1/192.0,
	};

	double filterZ[26] = {
		 1/6.0,
		 1/6.0,
		 1/6.0,
		 1/6.0,
		-1/6.0,
		-1/6.0,
		-1/6.0,
		-1/6.0,

		 1/16.0,
		-1/16.0,

		-1/32.0,
		 1/32.0,
		 1/32.0,
		-1/32.0,

		-1/32.0,
		 1/32.0,
		 1/32.0,
		-1/32.0,

		-1/192.0,
		-1/192.0,
		-1/192.0,
		-1/192.0,
		 1/192.0,
		 1/192.0,
		 1/192.0,
		 1/192.0,
	};

	elem_t sumX = 0, sumY = 0, sumZ = 0;
	for( int i=0;i<26;++i){
		elem_t pX=0,pY=0,pZ = 0;

		BCC_DELTA( pX, x,y,z, domainX[i][0], domainX[i][1], domainX[i][2] );
		BCC_DELTA( pY, x,y,z, domainY[i][0], domainY[i][1], domainY[i][2] );
		BCC_DELTA( pZ, x,y,z, domainZ[i][0], domainZ[i][1], domainZ[i][2] );

		sumX += pX * filterX[i];
		sumY += pY * filterY[i];
		sumZ += pZ * filterZ[i];
	}

	N[X] = sumX / m_halfScale[X];
	N[Y] = sumY / m_halfScale[Y];
	N[Z] = sumZ / m_halfScale[Z];

	elem_t NF[3];
	NF[0] = ComputeNormalXOptimalFilter(x,y,z);
	NF[1] = ComputeNormalYOptimalFilter(x,y,z);
	NF[2] = ComputeNormalZOptimalFilter(x,y,z);

	if( ABS(NF[0] - N[0]) > EPSILON ) printf("OPT26: Error in X: Precomp=%lf OnTheFly=%lf\n", N[0],NF[0]);
	if( ABS(NF[1] - N[1]) > EPSILON ) printf("OPT26: Error in Y: Precomp=%lf OnTheFly=%lf\n", N[1],NF[1]);
	if( ABS(NF[2] - N[2]) > EPSILON ) printf("OPT26: Error in Z: Precomp=%lf OnTheFly=%lf\n", N[2],NF[2]);

}

elem_t  CBCCData::ComputeNormalXSizeOptimalFilter( int x,int y, int z) const
{

	int xm1 = x-1;
	int xp1 = x+1;
	int ym1 = y-1;
	int yp1 = y+1;
	int zm1 = z-1;
	int zp1 = z+1;
	int zm2 = z-2;
	int zp2 = z+2;

	if( z & 0x01) { //ODD
		return m_invHalfScale[X] * (
				(      0.020833333    ) * (
				(Fetch(xm1,ym1,zp2)+Fetch(xm1,yp1,zp2)+Fetch(xm1,ym1,zm2)+Fetch(xm1,yp1,zm2))
				-
				(Fetch(xp1,ym1,zp2)+Fetch(xp1,yp1,zp2)+Fetch(xp1,ym1,zm2)+Fetch(xp1,yp1,zm2))
				)
				+
				(      0.166666667    )  * (
				(Fetch(xp1,y,zp1)+Fetch(xp1,yp1,zp1)+Fetch(xp1,y,zm1)+Fetch(xp1,yp1,zm1))
				-
				(Fetch(x,y,zp1)+Fetch(x,yp1,zp1)+Fetch(x,y,zm1)+Fetch(x,yp1,zm1))
					)
					   );
	}
	else{
		return m_invHalfScale[X] * (
				(      0.020833333    ) * (
				(Fetch(xm1,ym1,zp2)+Fetch(xm1,yp1,zp2)+Fetch(xm1,ym1,zm2)+Fetch(xm1,yp1,zm2))
				    -
				(Fetch(xp1,ym1,zp2)+Fetch(xp1,yp1,zp2)+Fetch(xp1,ym1,zm2)+Fetch(xp1,yp1,zm2))
					  )
				+

				(      0.166666667     ) * (
				(Fetch(x,ym1,zp1)+Fetch(x,y,zp1)+Fetch(x,ym1,zm1)+Fetch(x,y,zm1))
				-
				(Fetch(xm1,ym1,zp1)+Fetch(xm1,y,zp1)+Fetch(xm1,ym1,zm1)+Fetch(xm1,y,zm1))
					   )
				    );
	}
}

elem_t  CBCCData::ComputeNormalYSizeOptimalFilter( int x,int y, int z) const
{
	int xm1 = x-1;
	int xp1 = x+1;
	int ym1 = y-1;
	int yp1 = y+1;
	int zm1 = z-1;
	int zp1 = z+1;
	int zm2 = z-2;
	int zp2 = z+2;

	if( z & 0x01) { //ODD

		return m_invHalfScale[Y] * ((      0.020833333    ) *
				((Fetch(xp1,ym1,zp2)+Fetch(xm1,ym1,zp2)+Fetch(xp1,ym1,zm2)+Fetch(xm1,ym1,zm2)) -
				(Fetch(xp1,yp1,zp2)+Fetch(xm1,yp1,zp2)+Fetch(xp1,yp1,zm2)+Fetch(xm1,yp1,zm2)))
				+
					   (     0.166666667    ) * (
				(Fetch(xp1,yp1,zp1)+Fetch(x,yp1,zp1)+Fetch(xp1,yp1,zm1)+Fetch(x,yp1,zm1)) -
				(Fetch(xp1,y,zp1)+Fetch(x,y,zp1)+Fetch(xp1,y,zm1)+Fetch(x,y,zm1))) );


	}
	else{

		return m_invHalfScale[Y] * ((      0.020833333    ) *
				((Fetch(xp1,ym1,zp2)+Fetch(xm1,ym1,zp2)+Fetch(xp1,ym1,zm2)+Fetch(xm1,ym1,zm2)) -
				(Fetch(xp1,yp1,zp2)+Fetch(xm1,yp1,zp2)+Fetch(xp1,yp1,zm2)+Fetch(xm1,yp1,zm2)))
				+
				(      0.166666667     ) * (
				(Fetch(x,y,zp1)+Fetch(xm1,y,zp1)+Fetch(x,y,zm1)+Fetch(xm1,y,zm1))
				- (Fetch(x,ym1,zp1)+Fetch(xm1,ym1,zp1)+Fetch(x,ym1,zm1)+Fetch(xm1,ym1,zm1))));
	}
}

elem_t  CBCCData::ComputeNormalZSizeOptimalFilter( int x,int y, int z) const
{
	int xm1 = x-1;
	int xp1 = x+1;
	int ym1 = y-1;
	int yp1 = y+1;
	int zm1 = z-1;
	int zp1 = z+1;
	int zm2 = z-2;
	int zp2 = z+2;

	if( z & 0x01) { //ODD
		return m_invHalfScale[Z] * ((      0.020833333    ) *
				((Fetch(xm1,ym1,zm2)+Fetch(xm1,yp1,zm2)+Fetch(xp1,ym1,zm2)+Fetch(xp1,yp1,zm2)) -
				(Fetch(xm1,ym1,zp2)+Fetch(xm1,yp1,zp2)+Fetch(xp1,ym1,zp2)+Fetch(xp1,yp1,zp2)))
				+
				(      0.166666667     ) *
				((Fetch(x,y,zp1)+Fetch(x,yp1,zp1)+Fetch(xp1,y,zp1)+Fetch(xp1,yp1,zp1))-
				(Fetch(x,y,zm1)+Fetch(x,yp1,zm1)+Fetch(xp1,y,zm1)+Fetch(xp1,yp1,zm1))));
	}
	else{
		return m_invHalfScale[Z] * ((      0.020833333    ) *
				((Fetch(xm1,ym1,zm2)+Fetch(xm1,yp1,zm2)+Fetch(xp1,ym1,zm2)+Fetch(xp1,yp1,zm2))-(Fetch(xm1,ym1,
				  zp2)+Fetch(xm1,yp1,zp2)+Fetch(xp1,ym1,zp2)+Fetch(xp1,yp1,zp2)))
				+
				(      0.166666667     ) * ((Fetch(xm1,ym1,zp1)+Fetch(xm1,y,zp1)+Fetch(x,ym1,zp1)+Fetch(x,y,zp1))-
				(Fetch(xm1,ym1,zm1)+Fetch(xm1,y,zm1)+Fetch(x,ym1,zm1)+Fetch(x,y,zm1))));
	}

}

void   CBCCData::PrecompNormalSizeOptimalFilter( vector3_t N, int x, int y, int z )
{
	int domain[16][3] = {
		{ -1, -1, 1,}, {1,-1,1}, {1,1,1}, {-1,1,1},
        { -1, -1, -1}, {1,-1,-1}, {1,1,-1}, {-1,1,-1},
        { -2, -2, 2,}, {2,-2,2}, {2,2,2}, {-2,2,2},
        { -2, -2, -2,}, {2,-2,-2}, {2,2,-2}, {-2,2,-2},
    };



	double filterX[16] = {-1.0/6.0, 1.0/6.0, 1.0/6.0, -1.0/6.0, -1.0/6.0, 1.0/6.0, 1.0/6.0, -1.0/6.0, 1.0/48.0, -1.0/48.0, -1.0/48.0, 1.0/48.0, 1.0/48.0, -1.0/48.0, -1.0/48.0, 1.0/48.0,};
	double filterY[16] = {-1.0/6.0, -1.0/6.0, 1.0/6.0, 1.0/6.0, -1.0/6.0, -1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/48.0, 1.0/48.0, -1.0/48.0, -1.0/48.0, 1.0/48.0, 1.0/48.0, -1.0/48.0, -1.0/48.0,};
	double filterZ[16] = {1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, -1.0/6.0, -1.0/6.0, -1.0/6.0, -1.0/6.0, -1.0/48.0, -1.0/48.0, -1.0/48.0, -1.0/48.0, 1.0/48.0, 1.0/48.0, 1.0/48.0, 1.0/48.0,};


	elem_t sumX = 0, sumY = 0, sumZ = 0;
	for( int i=0;i<16;++i){
		elem_t p = 0;
		BCC_DELTA( p, x,y,z, domain[i][0], domain[i][1], domain[i][2] );
		sumX += p * filterX[i];
		sumY += p * filterY[i];
		sumZ += p * filterZ[i];
	}

	N[X] = sumX / m_halfScale[X];
	N[Y] = sumY / m_halfScale[Y];
	N[Z] = sumZ / m_halfScale[Z];

}

elem_t CBCCData::ComputeNormalXCentralDifference( int x,int y, int z) const
{
	return 0.25 * m_invHalfScale[X] * (Fetch(x+1,y,z) - Fetch(x-1,y,z));
}

elem_t CBCCData::ComputeNormalYCentralDifference( int x,int y, int z) const
{
	return 0.25 * m_invHalfScale[Y] * (Fetch(x,y+1,z) - Fetch(x,y-1,z));
}

elem_t CBCCData::ComputeNormalZCentralDifference( int x,int y, int z) const
{
	return 0.25 * m_invHalfScale[Z] * (Fetch(x,y,z+2) - Fetch(x,y,z-2));
}

void CBCCData::PrecompNormalCentralDifference(vector3_t N, int x, int y, int z)
{
	elem_t q[6] = {0};

	if( z+2 < m_zSize ) q[0] = Fetch( x,y,z+2);
	if( z-2 >= 0 ) q[1] = Fetch( x,y,z-2);
	if( x+1 < m_xSize ) q[2] = Fetch( x+1,y,z);
	if( x-1 >= 0 ) q[3] = Fetch( x-1,y,z);
	if( y+1 < m_ySize ) q[4] = Fetch( x,y+1,z);
	if( y-1 >= 0 ) q[5] = Fetch( x,y-1,z);


	N[X] = (q[2] - q[3]) / (4.0f * m_halfScale[X]);
	N[Y] = (q[4] - q[5]) / (4.0f * m_halfScale[Y]);
	N[Z] = (q[0] - q[1]) / (4.0f * m_halfScale[Z]);
}

void CBCCData::CalculateNormals( void (*fnNormal)(const vector3_t p, vector3_t N) )
{

 	if( m_bOnTheFlyNormals || m_bEpsilonNormalMode ){
		printf("CBCCData: Skipping precomputation of gradients\n");
		m_sepNormals[X] = NULL;
		m_sepNormals[Y] = NULL;
		m_sepNormals[Z] = NULL;
		return;
	}
	
	printf("Precalculating gradients\n");


	FILE *fp = NULL;

	m_sepNormals[X] = new elem_t ** [m_zSize+2*GUARD_BAND];
	m_sepNormals[Y] = new elem_t ** [m_zSize+2*GUARD_BAND];
	m_sepNormals[Z] = new elem_t ** [m_zSize+2*GUARD_BAND];

	for(int z=0;z<(m_zSize+2*GUARD_BAND);++z){

		m_sepNormals[X][z] = new elem_t * [ m_ySize+2*GUARD_BAND ];
		m_sepNormals[Y][z] = new elem_t * [ m_ySize+2*GUARD_BAND ];
		m_sepNormals[Z][z] = new elem_t * [ m_ySize+2*GUARD_BAND ];

		for(int y=0;y<(m_ySize+2*GUARD_BAND);++y){

			m_sepNormals[X][z][y] = new elem_t [ m_xSize+2*GUARD_BAND];
			m_sepNormals[Y][z][y] = new elem_t [ m_xSize+2*GUARD_BAND ];
			m_sepNormals[Z][z][y] = new elem_t [ m_xSize+2*GUARD_BAND];

			memset( m_sepNormals[X][z][y], 0, sizeof( elem_t) * (m_xSize +2*GUARD_BAND) );
			memset( m_sepNormals[Y][z][y], 0, sizeof( elem_t) * (m_xSize +2*GUARD_BAND) );
			memset( m_sepNormals[Z][z][y], 0, sizeof( elem_t) * (m_xSize +2*GUARD_BAND) );
		}
	}


	if( fnNormal ){

		if( m_gridNormalOutputFile[0] )
			fp = fopen( m_gridNormalOutputFile, "w" );
		else
			fp = fopen( "BCCData_Normal_Data.txt","w");
	}



	for(int x=0;x<m_xSize;x++){

		for(int y=0;y<m_ySize;++y){


			for(int z=0;z<m_zSize;++z){


				vector3_t &v = m_pppNormals[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND];

				switch (m_normalEst){
					case BCC_NORMAL_SOCD:
						PrecompNormalCentralDifference( v, x,y,z );
						break;
					case BCC_NORMAL_ICD:
						PrecompNormalInterpolatedCentralDifference(v, x,y,z );
						break;
					case BCC_NORMAL_BCD:
						PrecompNormalBoxDifference(v,x,y,z);
						break;
					case BCC_NORMAL_SIZE_OPTIMAL16:
						PrecompNormalSizeOptimalFilter( v, x,y,z );
						break;
					case BCC_NORMAL_OPTIMAL26:
						PrecompNormalOptimalFilter(v,x,y,z);
						break;
					case BCC_NORMAL_SEMI_OPTIMAL18:
						PrecompNormalSemiOptimalFilter(v,x,y,z);
						break;
				};

				if( fnNormal ){

					int index[3] = {x,y,z};

					int bccX,bccY,bccZ;

					IndexToBCC( bccX,bccY,bccZ, index );

					vector2_t screenSpace;
					vector3_t worldSpace;
					vector3_t bcc = { bccX, bccY, bccZ };

					BCCToWorld( worldSpace, bcc );

					camGetScreenSpacePointFromWorldSpace( screenSpace, worldSpace );

					vector3_t N;

					fnNormal( worldSpace, N );

					fprintf(fp,"%d %d %d; %d %d; %f %f %f; %f %f %f;\n",x,y,z, (int)(screenSpace[X]+0.5),(int)(screenSpace[Y]+0.5), N[X], N[Y], N[Z], v[X],v[Y],v[Z]);

				}
			
				m_sepNormals[X][z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND] = v[X];
				m_sepNormals[Y][z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND] = v[Y];
				m_sepNormals[Z][z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND] = v[Z];


			}
		}
	}
	

	if( fnNormal ){
		fclose( fp );
	}

}
/*
template<typename T, typename V, typename SV>
		V vuVolumeInterpolator<T,V,SV>::getValueBCCBoxsplineCubicTurbo(const vuVector& pos) const {

			vuLattice::MultiIndex id;

			float x = pos[0], y = pos[1], z = pos[2];
			register float alpha = .5f*(x+y), beta = .5f*(x+z), gamma = .5f*(y+z);
			int alpha0 = (int)floor(alpha); alpha = alpha - alpha0;
			int beta0 = (int)floor(beta); beta = beta - beta0;
			int gamma0 = (int)floor(gamma); gamma = gamma - gamma0;
			id[0] = (alpha0+beta0-gamma0);id[1]=(alpha0-beta0+gamma0);id[2]=(-alpha0+beta0+gamma0);
  //  int x0 = (alpha0+beta0-gamma0),y0=(alpha0-beta0+gamma0),z0=(-alpha0+beta0+gamma0);

			float p1[8],p2[8],p3[8],p4[8];
			register float mymax,mymid,mymin;


			if(ODD(id[2]))
			{
				--id[0]>>=1; --id[1]>>=1;
				if(alpha >= beta)
				{
					if(beta >= gamma)
					{
						mymax = alpha; mymid = beta; mymin = gamma;
						RegionOne_FillData_Odd;
					}
					else if(alpha >= gamma)
					{
						mymax = alpha; mymid = gamma; mymin = beta;
						RegionTwo_FillData_Odd;
					}
					else
					{
						mymax = gamma; mymid = alpha; mymin = beta;
						RegionThree_FillData_Odd;
					}
				}
				else
				{
					if(alpha >= gamma)
					{
						mymax = beta; mymid = alpha; mymin = gamma;
						RegionFour_FillData_Odd;
					}
					else if(beta >= gamma)
					{
						mymax = beta; mymid = gamma; mymin = alpha;
						RegionFive_FillData_Odd;
					}
					else
					{
						mymax = gamma; mymid = beta; mymin = alpha;
						RegionSix_FillData_Odd;
					}
				}
			}
			else
			{
				id[0]>>=1; id[1]>>=1;
				if(alpha >= beta)
				{
					if(beta >= gamma)
					{
						mymax = alpha; mymid = beta; mymin = gamma;
						RegionOne_FillData_Even;
					}
					else if(alpha >= gamma)
					{
						mymax = alpha; mymid = gamma; mymin = beta;
						RegionTwo_FillData_Even;
					}
					else
					{
						mymax = gamma; mymid = alpha; mymin = beta;
						RegionThree_FillData_Even;
					}
				}
				else
				{
					if(alpha >= gamma)
					{
						mymax = beta; mymid = alpha; mymin = gamma;
						RegionFour_FillData_Even;
					}
					else if(beta >= gamma)
					{
						mymax = beta; mymid = gamma; mymin = alpha;
						RegionFive_FillData_Even;
					}
					else
					{
						mymax = gamma; mymid = beta; mymin = alpha;
						RegionSix_FillData_Even;
					}
				}
			}

			register float result = 0.0f;
			alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f;
			CONV_PPIPED(result,p1,alpha,beta,gamma);
			alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f;
			CONV_PPIPED(result,p2,alpha,beta,gamma);
			alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax);
			CONV_PPIPED(result,p3,alpha,beta,gamma);
			alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f);
			CONV_PPIPED(result,p4,alpha,beta,gamma);


			return result;
		}
		*/

	int argTable[8][4][6] = 
		{
			//0
			{			
			  {-1,-1,-1,3,2,1},
			  {1,1,-1,1,2,3},
			  {-1,1,1,1,2,3},
			  {1,-1,1,3,2,1},
			},
			
			{
				//1,
			},
			
			//2
			{
				{-1,-1,-1,2,3,1},
				{1,1,-1,2,1,3},
				{1,-1,1,2,1,3},
				{-1,1,1,2,3,1},
			},
			
			//3
			{
				{-1,-1,-1,1,3,2},
				{-1,1,1,3,1,2},
				{1,-1,1,3,1,2},
				{1,1,-1,1,3,2},
			},
			
			//4
			{
				{-1,-1,-1,3,1,2},
				{1,-1,1,1,3,2},
				{-1,1,1,1,3,2},
				{1,1,-1,3,1,2},
			},
			
			//5
			{
				{-1,-1,-1,2,1,3},
			    {1,-1,1,2,3,1},
			    {1,1,-1,2,3,1},
			    {-1,1,1,2,1,3},
			},
			
			{
				//6
			},
			
			//7
			{
				{-1,-1,-1,1,2,3},			  
				{-1,1,1,3,2,1},
				{ 1,1,-1,3,2,1},
				{ 1,-1,1,1,2,3},
			},
						
		};
		
		int offsets[8][3][3] = {
			{
				//0
				{1,1,1},
				{-2,0,0},
				{1,-1,1},
			},
			
			{
				//1
			},
			
			{
				//2
				{1,1,1},
				{0,-2,0},
				{-1,1,1},
			},
			
			{
				//3
				{1,1,1},
				{0,-2,0},
				{1,1,-1},				
			},
			
			{
				//4
				{1,1,1},
				{-2,0,0},
				{1,1,-1},
			},
			
			{
				//5
				{1,1,1},
				{0,0,-2},
				{-1,1,1},
			},
			
			{
				//6
			},
			
			{
				//7
				{1,1,1},
				{0,0,-2},
				{1,-1,1},
			},
		};
		
	int minmaxTable[8][3] = {
		{2,1,0 }, //0
		{}, // 1
		{1, 2, 0}, //2
		{1,0,2}, //3
		{2,0,1}, //4
		{0,2,1}, //5
		{}, //6
		{0,1,2} //7
	};


double CBCCData::InterpolateQuinticBoxSplineSteroid( const vector3_t bcc) const
{
	
	int index3D[3];
	
	register double abg[3] = {.5f*(bcc[X]+bcc[Y]), .5f*(bcc[X]+bcc[Z]), .5f*(bcc[Y]+bcc[Z])};
	
	int alpha0 = (int)floor(abg[0]); abg[0] = abg[0] - alpha0;
	int beta0 = (int)floor(abg[1]); abg[1] = abg[1]- beta0;
	int gamma0 = (int)floor(abg[2]); abg[2] = abg[2] - gamma0;
	int x0 =
	(alpha0+beta0-gamma0),y0=(alpha0-beta0+gamma0),z0=(-alpha0+beta0+gamma0);
		
	double p[4][8];
	register double result = 0.0f;
				
	register double mymax,mymid,mymin;
	
	int alpha_GE_beta  = abg[0] >= abg[1];
	int beta_GE_gamma  = abg[1] >= abg[2];
	int alpha_GE_gamma = abg[0] >= abg[2];			

	
	int i = (alpha_GE_beta << 2) + (beta_GE_gamma << 1) + alpha_GE_gamma;
	
	int *pArgTable    = &argTable[i][0][0];
	int *pMinMaxTable = &minmaxTable[i][0];
	int *pOffsetTable = &offsets[i][0][0];
	double *pp = &p[0][0];
	
	mymax = abg[*pMinMaxTable++];
	mymid = abg[*pMinMaxTable++],
	mymin = abg[*pMinMaxTable++];	
	
	for(int l=0;l<3;++l){	
		
// 		FILL_PPIPED(p[l],argTable[i][l][0],argTable[i][l][1],argTable[i][l][2],
// 				   argTable[i][l][3],argTable[i][l][4],argTable[i][l][5]);
				   
		FILL_PPIPED(pp,pArgTable[0],pArgTable[1],pArgTable[2],
 				   pArgTable[3],pArgTable[4],pArgTable[5]);
		
		pArgTable += 6;		   
		pp += 8;
	    
// 		x0 += offsets[i][l][0];
// 		y0 += offsets[i][l][1];
// 		z0 += offsets[i][l][2];

	    x0 += *pOffsetTable++;
		y0 += *pOffsetTable++;
		z0 += *pOffsetTable++;
	}
	
	FILL_PPIPED(pp,pArgTable[0],pArgTable[1],pArgTable[2],
 				   pArgTable[3],pArgTable[4],pArgTable[5]);
	

	pp = &p[0][0];
	double alpha,beta,gamma;
	
	
	alpha=mymax-1.0;beta=mymid-1.0;gamma=mymin-1.0;
	CONV_PPIPED(result,pp,alpha,beta,gamma);
	pp+=8;
	alpha=-mymin;beta=mymax-mymin-1.0;gamma=mymid-mymin-1.0;
	CONV_PPIPED(result,pp,alpha,beta,gamma);
	pp+=8;
	alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax);	
	CONV_PPIPED(result,pp,alpha,beta,gamma);
	pp+=8;
	alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0);
	CONV_PPIPED(result,pp,alpha,beta,gamma);	
		
	return (result);
}

// #define INCADDVECTOR(v,v0,v1){ \
// 	v[0] += v0[0] + v1[0]; \
// 	v[1] += v0[1] + v1[1]; \
// 	v[2] += v0[2] + v1[2]; \
// }
//    
// #define VECTOR(v,x,y,z) {v[0]=(x);v[1]=(y);v[2]=(z);}   
// 
// double CBCCData::InterpolateQuinticBoxSplineSteroid( const vector3_t pos) const
// {
// 	
// 	int permuteMatrix[8][8] = {
// 		{0,3,2,1,4,5,6,7},
// 		{0,1,2,3,4,5,6,7}, // This case should never be used
// 		{0,2,3,1,5,4,6,7},
// 		{0,1,3,2,6,4,5,7},
// 		{0,3,1,2,4,6,5,7},
// 		{0,2,1,3,5,6,4,7},
// 		{0,1,2,3,4,5,6,7}, // This case should never be used
// 		{0,1,2,3,6,5,4,7},
// 	};
// 
// 	int p3Offsets[8][3] ={
// 		{-1,1,1},  //0
// 		{0,0,0},   //1 // this should not be used
// 		{1,-1,1},  //2 
// 		{1,-1,1},  //3
// 		{-1,1,1},  //4
// 		{1,1,-1},  //5
// 		{0,0,0},   //6 // this should not be used
// 		{1,1,-1},  //7
// 	};
// 	
// 	int p4Offsets[8][3] = {
// 		{0,0,2},  //0
// 		{0,0,0},   //1 // this should not be used
// 		{0,0,2},  //2 
// 		{2,0,0},  //3
// 		{0,2,0},  //4
// 		{0,2,0},  //5
// 		{0,0,0},   //6 // this should not be used
// 		{2,0,0},  //7
// 	};
// 	
// 	int q2Offsets[8][3] = {
// 		{1,1,-1},  //0
// 		{0,0,0},    //1 // this should not be used
// 		{1,1,-1},  //2 
// 		{-1,1,1},  //3
// 		{1,-1,1},  //4
// 		{1,-1,1},  //5
// 		{0,0,0},    //6 // this should not be used
// 		{-1,1,1},  //7
// 	};
// 	
// 	int q4Offsets[8][3] = {
// 		{1,-1,1},  //0
// 		{0,0,0},    //1 // this should not be used
// 		{-1,1,1},  //2 
// 		{1,1,-1},  //3
// 		{1,1,-1},  //4
// 		{-1,1,1},  //5
// 		{0,0,0},    //6 // this should not be used
// 		{1,-1,1},  //7
// 	};
// 	
// 	
// 	int x0,y0,z0;
// 	int x1,y1,z1;
// 	
// 	double x = pos[0], y = pos[1], z = pos[2];
// 	register double alpha = 0.5*(x+y), beta = 0.5*(x+z), gamma = 0.5*(y+z);
// 	int alpha0 = (int)floor(alpha); alpha = alpha - alpha0;
// 	int beta0 = (int)floor(beta); beta = beta - beta0;
// 	int gamma0 = (int)floor(gamma); gamma = gamma - gamma0;
// 	
// 	x0 = x1 = alpha0 + beta0 - gamma0;
// 	y0 = y1 = alpha0-beta0+gamma0;
// 	z0 = z1 = beta0+gamma0-alpha0;
// 	
// 	int index3D[3];
// 	
// 	double p[4][8];
// 	
// 	//int Q2[3],Q3[3],Q4[3];	
// 	
// 	// variables needed for the sorting
// 	int alpha_GE_beta = alpha >= beta;
// 	int beta_GE_gamma = beta >= gamma;
// 	int alpha_GE_gamma = alpha >= gamma;
// 	
// 	double mymax, mymid, mymin;
// 	// sorting
// 	
// 	mymax = max (alpha, max (beta, gamma));
// 	mymin = min (alpha, min (beta, gamma));
// 	mymid = (alpha + beta + gamma) - mymax - mymin;
// 		
// 	
// 	int i = (alpha_GE_beta << 2) + (beta_GE_gamma << 1) + alpha_GE_gamma;
// 	
// // 	int I1[3],I2[3];
// // 	
// // 	I1[0] = (i == 7);I1[1] = (i == 5); I1[2] = (i == 4);
// // 	
// // 	I2[0] = (i == 3); I2[1] = (i == 2); I2[2] = (i == 0);
// // 	
// 	int P2[3] = {1,1,1};	
// 	int P3[3],P4[3];
// 	int Q1[3] = {-1,-1,-1};	
// 	int Q2[3],Q4[3];
// 	
// 	for(int k=0;k<3;++k){
// 		P3[k] = p3Offsets[i][k];
// 		P4[k] = p4Offsets[i][k];
// 		Q2[k] = q2Offsets[i][k];
// 		Q4[k] = q4Offsets[i][k];
// 	}
// 
// 	int Q3[3] = {P3[0],P3[1],P3[2]};
// 	
// // 	int P3[3] = {1,1,-1};
// // 		
// // 	int c0 = I1[2] + I2[2];
// // 	int c1 = I2[0] + I2[1];
// // 	int v0[3] = {c0*-2,0,2*c0};
// // 	int v1[3] = {0,c1*-2,c1*2};
// // 	INCADDVECTOR(P3,v0,v1);
// 	
// // 	int P4[3] = {2,0,0};
// // 	c0 = I1[1] + I1[2];
// // 	c1 = I2[1] + I2[2];
// // 	
// // 	VECTOR(v0,-2*c0,2*c0,0);
// // 	VECTOR(v1,c1*-2,0,c1*2);
// // 	INCADDVECTOR(P4,v0,v1);
// // 	
// // 	
// // 	int Q2[3] = {-1,1,1};	
// // 	VECTOR(v0,2*c0,-2*c0,0);
// // 	VECTOR(v1,c1*2,0,c1*-2);
// // 	INCADDVECTOR(Q2,v0,v1);
// // 	
// // 	int Q4[3] = {1,-1,1};	
// // 	c0 = I1[1] + I2[1];
// // 	c1 = I1[2] + I2[0];	
// // 	VECTOR(v0,-2*c0,2*c0,0);
// // 	VECTOR(v1,0,c1*2,c1*-2);
// // 	INCADDVECTOR(Q4,v0,v1);
// 	
// 	
// 		
// 	FILL_PPIPED_STEROID(p[0],Q1);
// 	
// 	x0 = x1+P2[0]; y0 = y1+P2[1]; z0 = z1+P2[2];
// 	FILL_PPIPED_STEROID(p[1],Q2);
// 	
// 	x0 = x1+P3[0]; y0 = y1+P3[1]; z0 = z1+P3[2];
// 	FILL_PPIPED_STEROID(p[2],Q3);
// 	
// 	x0 = x1+P4[0]; y0 = y1+P4[1]; z0 = z1+P4[2];
// 	FILL_PPIPED_STEROID(p[3],Q4);
// 	
// 	int *permVector = permuteMatrix[i];
// 	
// 	double p_perm[4][8];
// 	for(int k=0;k<4;++k){
// 		Permute( p_perm[k], p[k], permVector );
// 	}
// 	
// 	// do the convolution. Before each CONV_PPIPED call
// 	// mymin, mymid and mymax are used to transform
// 	// each parallelepiped and tetrahedron of focus.
// 	double result = 0.0;
// 	alpha=mymax-1.0; beta=mymid-1.0; gamma=mymin-1.0;
// 	CONV_PPIPED(result, p_perm[0], alpha, beta, gamma);
// 	alpha=-mymin; beta=mymax-mymin-1.0; gamma=mymid-mymin-1.0;
// 	CONV_PPIPED(result, p_perm[1], alpha, beta, gamma);
// 	alpha=(-mymax+mymid); beta=(-mymax+mymin); gamma=(-mymax);
// 	CONV_PPIPED(result, p_perm[2], alpha, beta, gamma);
// 	alpha=(-mymid+mymin); beta=(-mymid); gamma=(mymax-mymid-1.0);
// 	CONV_PPIPED(result, p_perm[3], alpha, beta, gamma);
// 	
// 	if( ((abs(pos[0])- 32) < 0.001) &&  
// 		((abs(pos[1])-32) < 0.01) &&  
// 		((abs(pos[2])-32) < 0.01) )
// 	{
// 		printf("%f\n",result);
// 	}
// 	return result;	
// }

double CBCCData::InterpolateLinearBoxSplineTurbo( const vector3_t pos ) const
{
		
	int id[3];


	double x = pos[0], y = pos[1], z = pos[2];
	register double alpha = .5f*(x+y), beta = .5f*(x+z), gamma = .5f*(y+z);
	int alpha0 = (int)floor(alpha); alpha = alpha - alpha0;
	int beta0 = (int)floor(beta); beta = beta - beta0;
	int gamma0 = (int)floor(gamma); gamma = gamma - gamma0;
	id[0] = (alpha0+beta0-gamma0);id[1]=(alpha0-beta0+gamma0);id[2]=(-alpha0+beta0+gamma0);
//  int x0 = (alpha0+beta0-gamma0),y0=(alpha0-beta0+gamma0),z0=(-alpha0+beta0+gamma0);


	double p1,p2,p3,p4;
	register double mymax,mymid,mymin;


	if(ODD(id[2]))
	{
		--id[0]>>=1; --id[1]>>=1;
		if(alpha >= beta)
		{
			if(beta >= gamma)
			{
				mymax = alpha; mymid = beta; mymin = gamma;
				p1=Fetch( id[0], id[1], id[2] );
				id[0]++;        p4=Fetch( id[0], id[1], id[2] );
				id[1]++;id[2]++;        p2=Fetch( id[0], id[1], id[2] );
				id[2]-=2;       p3=Fetch( id[0], id[1], id[2] );
			}
			else if(alpha >= gamma)
			{
				mymax = alpha; mymid = gamma; mymin = beta;
				p1=Fetch( id[0], id[1], id[2] );
				id[1]++;        p4=Fetch( id[0], id[1], id[2] );
				id[0]++;id[2]++;        p2=Fetch( id[0], id[1], id[2] );
				id[2]-=2;       p3=Fetch( id[0], id[1], id[2] );
			}
			else
			{
				mymax = gamma; mymid = alpha; mymin = beta;
				p1=Fetch( id[0], id[1], id[2] );
				id[1]++;        p4=Fetch( id[0], id[1], id[2] );
				id[2]++;        p3=Fetch( id[0], id[1], id[2] );
				id[0]++;        p2=Fetch( id[0], id[1], id[2] );
			}
		}
		else
		{
			if(alpha >= gamma)
			{
				mymax = beta; mymid = alpha; mymin = gamma;
				p1=Fetch( id[0], id[1], id[2] );
				id[0]++;        p4=Fetch( id[0], id[1], id[2] );
				id[2]++;        p3=Fetch( id[0], id[1], id[2] );
				id[1]++;        p2=Fetch( id[0], id[1], id[2] );
			}
			else if(beta >= gamma)
			{
				mymax = beta; mymid = gamma; mymin = alpha;
				p1=Fetch( id[0], id[1], id[2] );
				id[2]+=2;        p4=Fetch( id[0], id[1], id[2] );
				id[0]++;id[2]--;        p3=Fetch( id[0], id[1], id[2] );
				id[1]++;        p2=Fetch( id[0], id[1], id[2] );
	
			}
			else
			{
				mymax = gamma; mymid = beta; mymin = alpha;
				p1=Fetch( id[0], id[1], id[2] );
				id[2]+=2;        p4=Fetch( id[0], id[1], id[2] );
				id[1]++;id[2]--;        p3=Fetch( id[0], id[1], id[2] );
				id[0]++;        p2=Fetch( id[0], id[1], id[2] );

			}
		}
	}
	else
	{
		id[0]>>=1; id[1]>>=1;
		if(alpha >= beta)
		{
			if(beta >= gamma)
			{
				mymax = alpha; mymid = beta; mymin = gamma;
				p1=Fetch( id[0], id[1], id[2] );
				id[0]++;        p4=Fetch( id[0], id[1], id[2] );
				id[0]--;id[2]++;        p2=Fetch( id[0], id[1], id[2] );
				id[2]-=2;       p3=Fetch( id[0], id[1], id[2] );
			}
			else if(alpha >= gamma)
			{
				mymax = alpha; mymid = gamma; mymin = beta;
				p1=Fetch( id[0], id[1], id[2] );
				id[1]++;        p4=Fetch( id[0], id[1], id[2] );
				id[1]--;id[2]++;        p2=Fetch( id[0], id[1], id[2] );
				id[2]-=2;       p3=Fetch( id[0], id[1], id[2] );
			}
			else
			{
				mymax = gamma; mymid = alpha; mymin = beta;
				p1=Fetch( id[0], id[1], id[2] );
				id[1]++;        p4=Fetch( id[0], id[1], id[2] );
				id[1]--;id[2]++;        p2=Fetch( id[0], id[1], id[2] );
				id[0]--;        p3=Fetch( id[0], id[1], id[2] );
			}
		}
		else
		{
			if(alpha >= gamma)
			{
				mymax = beta; mymid = alpha; mymin = gamma;
				p1=Fetch( id[0], id[1], id[2] );
				id[0]++;        p4=Fetch( id[0], id[1], id[2] );
				id[0]--;id[2]++;        p2=Fetch( id[0], id[1], id[2] );
				id[1]--;        p3=Fetch( id[0], id[1], id[2] );
			}
			else if(beta >= gamma)
			{
				mymax = beta; mymid = gamma; mymin = alpha;
				p1=Fetch( id[0], id[1], id[2] );
				id[2]++;        p2=Fetch( id[0], id[1], id[2] );
				id[1]--;        p3=Fetch( id[0], id[1], id[2] );
				id[1]++;id[2]++;        p4=Fetch( id[0], id[1], id[2] );
			}
			else
			{
				mymax = gamma; mymid = beta; mymin = alpha;
				p1=Fetch( id[0], id[1], id[2] );
				id[2]++;        p2=Fetch( id[0], id[1], id[2] );
				id[0]--;        p3=Fetch( id[0], id[1], id[2] );
				id[0]++;id[2]++;        p4=Fetch( id[0], id[1], id[2] );
			}
		}
	}

	return p1 + mymax * (p3 - p1) + mymin * (p2 - p4) + mymid * (p4 - p3);  
}
			
// double CBCCData::InterpolateQuinticBoxSplineTurbo( const vector3_t pos ) const
// {
// 	int id[3];
// 		
// 	double x = pos[0], y = pos[1], z = pos[2];
// 	double alpha = .5f*(x+y), beta = .5f*(x+z), gamma = .5f*(y+z);
// 	int alpha0 = (int)floor(alpha); alpha = alpha - alpha0;
// 	int beta0 = (int)floor(beta); beta = beta - beta0;
// 	int gamma0 = (int)floor(gamma); gamma = gamma - gamma0;
// 	id[0] = (alpha0+beta0-gamma0);id[1]=(alpha0-beta0+gamma0);id[2]=(-alpha0+beta0+gamma0);
//   
// 	//  int x0 = (alpha0+beta0-gamma0),y0=(alpha0-beta0+gamma0),z0=(-alpha0+beta0+gamma0);
// 
// 	double p1[8],p2[8],p3[8],p4[8];
// 	double mymax,mymid,mymin;
// 
// 
// 	if(ODD(id[2]))
// 	{
// 		--id[0]>>=1; --id[1]>>=1;
// 		if(alpha >= beta)
// 		{
// 			if(beta >= gamma)
// 			{
// 				mymax = alpha; mymid = beta; mymin = gamma;
// 				RegionOne_FillData_Odd;
// 			}
// 			else if(alpha >= gamma)
// 			{
// 				mymax = alpha; mymid = gamma; mymin = beta;
// 				RegionTwo_FillData_Odd;
// 			}
// 			else
// 			{
// 				mymax = gamma; mymid = alpha; mymin = beta;
// 				RegionThree_FillData_Odd;
// 			}
// 		}
// 		else
// 		{
// 			if(alpha >= gamma)
// 			{
// 				mymax = beta; mymid = alpha; mymin = gamma;
// 				RegionFour_FillData_Odd;
// 			}
// 			else if(beta >= gamma)
// 			{
// 				mymax = beta; mymid = gamma; mymin = alpha;
// 				RegionFive_FillData_Odd;
// 			}
// 			else
// 			{
// 				mymax = gamma; mymid = beta; mymin = alpha;
// 				RegionSix_FillData_Odd;
// 			}
// 		}
// 	}
// 	else
// 	{
// 		id[0]>>=1; id[1]>>=1;
// 		if(alpha >= beta)
// 		{
// 			if(beta >= gamma)
// 			{
// 				mymax = alpha; mymid = beta; mymin = gamma;
// 				RegionOne_FillData_Even;
// 			}
// 			else if(alpha >= gamma)
// 			{
// 				mymax = alpha; mymid = gamma; mymin = beta;
// 				RegionTwo_FillData_Even;
// 			}
// 			else
// 			{
// 				mymax = gamma; mymid = alpha; mymin = beta;
// 				RegionThree_FillData_Even;
// 			}
// 		}
// 		else
// 		{
// 			if(alpha >= gamma)
// 			{
// 				mymax = beta; mymid = alpha; mymin = gamma;
// 				RegionFour_FillData_Even;
// 			}
// 			else if(beta >= gamma)
// 			{
// 				mymax = beta; mymid = gamma; mymin = alpha;
// 				RegionFive_FillData_Even;
// 			}
// 			else
// 			{
// 				mymax = gamma; mymid = beta; mymin = alpha;
// 				RegionSix_FillData_Even;
// 			}
// 		}
// 	}
// 
// 	double result = 0.0f;
// 	alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f;
// 	CONV_PPIPED(result,p1,alpha,beta,gamma);
// 	alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f;
// 	CONV_PPIPED(result,p2,alpha,beta,gamma);
// 	alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax);
// 	CONV_PPIPED(result,p3,alpha,beta,gamma);
// 	alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f);
// 	CONV_PPIPED(result,p4,alpha,beta,gamma);
// 
// 
// 	return result;
// 
// }

// void   CBCCData::Fill_1_Even(int id[],double p1[],double p2[], double p3[], double p4[],
// 							 double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const
// {
// 	*mymax = alpha; *mymid = beta; *mymin = gamma;
// 	RegionOne_FillData_Even
// }
// 
// void   CBCCData::Fill_2_Even(int id[],double p1[],double p2[], double p3[], double p4[],
// 							 double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const
// {
// 	*mymax = alpha; *mymid = gamma; *mymin = beta;
// 	RegionTwo_FillData_Even
// }
// 
// void   CBCCData::Fill_3_Even(int id[],double p1[],double p2[], double p3[], double p4[],
// 							 double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const
// {
// 	*mymax = gamma; *mymid = alpha; *mymin = beta;
// 	RegionThree_FillData_Even
// }
// 
// void   CBCCData::Fill_4_Even(int id[],double p1[],double p2[], double p3[], double p4[],
// 							 double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const
// {
// 	*mymax = beta; *mymid = alpha; *mymin = gamma;
// 	RegionFour_FillData_Even
// }
// 
// void   CBCCData::Fill_5_Even(int id[],double p1[],double p2[], double p3[], double p4[],
// 							 double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const
// {
// 	*mymax = beta; *mymid = gamma; *mymin = alpha;
// 	RegionFive_FillData_Even
// }
// 
// void   CBCCData::Fill_6_Even(int id[],double p1[],double p2[], double p3[], double p4[],
// 							 double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const
// {
// 	*mymax = gamma; *mymid = beta; *mymin = alpha;
// 	RegionSix_FillData_Even
// }
// 
// 
// void   CBCCData::Fill_1_Odd(int id[],double p1[],double p2[], double p3[], double p4[],
// 							double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const
// {
// 	*mymax = alpha; *mymid = beta; *mymin = gamma;
// 	RegionOne_FillData_Odd
// }
// 
// void   CBCCData::Fill_2_Odd(int id[],double p1[],double p2[], double p3[], double p4[],
// 							double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const
// {
// 	*mymax = alpha; *mymid = gamma; *mymin = beta;
// 	RegionTwo_FillData_Odd
// }
// 
// void   CBCCData::Fill_3_Odd(int id[],double p1[],double p2[], double p3[], double p4[],
// 							double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const
// {
// 	*mymax = gamma; *mymid = alpha; *mymin = beta;
// 	RegionThree_FillData_Odd
// }
// 
// void   CBCCData::Fill_4_Odd(int id[],double p1[],double p2[], double p3[], double p4[],
// 							double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const
// {
// 	*mymax = beta; *mymid = alpha; *mymin = gamma;
// 	RegionFour_FillData_Odd
// }
// 
// void   CBCCData::Fill_5_Odd(int id[],double p1[],double p2[], double p3[], double p4[],
// 							double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const
// {
// 	*mymax = beta; *mymid = gamma; *mymin = alpha;
// 	RegionFive_FillData_Odd
// }
// void   CBCCData::Fill_6_Odd(int id[],double p1[],double p2[], double p3[], double p4[], 
// 							double alpha, double beta,double gamma, double *mymax,double *mymid, double *mymin) const
// {
// 	*mymax = gamma; *mymid = beta; *mymin = alpha;
// 	RegionSix_FillData_Odd
// }
// 
// CBCCData::FNFILL CBCCData::m_fnFill[16] = {
// 	&CBCCData::Fill_1_Even,
// 	0,
// 	&CBCCData::Fill_2_Even,
// 	&CBCCData::Fill_3_Even,
// 	&CBCCData::Fill_4_Even,
// 	&CBCCData::Fill_5_Even,
// 	0,
// 	&CBCCData::Fill_6_Even,
// 	&CBCCData::Fill_1_Odd,
// 	0,
// 	&CBCCData::Fill_2_Odd,
// 	&CBCCData::Fill_3_Odd,
// 	&CBCCData::Fill_4_Odd,
// 	&CBCCData::Fill_5_Odd,
// 	0,
// 	&CBCCData::Fill_6_Odd
// };

double inline mypow2(double x)
{
	return x*x;
}

double inline mypow3(double x)
{
	return x*x*x;
}

#define CONV_PPIPED_OLD_INPLACE(value, p) \
   alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f; \
   value += (-10.0f*p[0]+4.0f*(p[1]+p[2]+p[3])-2.0f*(p[4]+p[5]+p[6])+p[7]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*p[0]-2.0f*(p[2]+p[3])+p[4]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*p[0]-2.0f*(p[1]+p[3])+p[5]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*p[0]-2.0f*(p[1]+p[2])+p[6]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*p[0]+p[3]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*p[0]+p[2]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*p[0]+p[1]) * rho22(gamma, alpha-1, beta-1) + \
   (p[0]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f; \
   value += (-10.0f*p[8]+4.0f*(p[9]+p[10]+p[11])-2.0f*(p[12]+p[13]+p[14])+p[15]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*p[8]-2.0f*(p[10]+p[11])+p[12]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*p[8]-2.0f*(p[9]+p[11])+p[13]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*p[8]-2.0f*(p[9]+p[10])+p[14]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*p[8]+p[11]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*p[8]+p[10]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*p[8]+p[9]) * rho22(gamma, alpha-1, beta-1) + \
   (p[8]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax); \
   value += (-10.0f*p[16]+4.0f*(p[17]+p[18]+p[19])-2.0f*(p[20]+p[21]+p[22])+p[23]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*p[16]-2.0f*(p[18]+p[19])+p[20]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*p[16]-2.0f*(p[17]+p[19])+p[21]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*p[16]-2.0f*(p[17]+p[18])+p[22]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*p[16]+p[19]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*p[16]+p[18]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*p[16]+p[17]) * rho22(gamma, alpha-1, beta-1) + \
   (p[16]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f); \
   value += (-10.0f*p[24]+4.0f*(p[25]+p[26]+p[27])-2.0f*(p[28]+p[29]+p[30])+p[31]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*p[24]-2.0f*(p[26]+p[27])+p[28]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*p[24]-2.0f*(p[25]+p[27])+p[29]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*p[24]-2.0f*(p[25]+p[26])+p[30]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*p[24]+p[27]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*p[24]+p[26]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*p[24]+p[25]) * rho22(gamma, alpha-1, beta-1) + \
   (p[24]) * rho22(alpha-1, beta-1, gamma-1);\
   
double CBCCData::InterpolateQuinticBoxSplineTurbo( const vector3_t pos ) const
{
	int id[3];
	
	double x = pos[0], y = pos[1], z = pos[2];
	double alpha = .5*(x+y), beta = .5*(x+z), gamma = .5*(y+z);
	//int alpha0 = (int)floor(alpha); alpha = alpha - alpha0;
	//int beta0 = (int)floor(beta); beta = beta - beta0;
	//int gamma0 = (int)floor(gamma); gamma = gamma - gamma0;
	
	int alpha0 = (int)(alpha); alpha = alpha - alpha0;
	int beta0 = (int)(beta); beta = beta - beta0;
	int gamma0 = (int)(gamma); gamma = gamma - gamma0;
	id[0] = (alpha0+beta0-gamma0);id[1]=(alpha0-beta0+gamma0);id[2]=(-alpha0+beta0+gamma0);
  
	//  int x0 = (alpha0+beta0-gamma0),y0=(alpha0-beta0+gamma0),z0=(-alpha0+beta0+gamma0);

	//double p1[8],p2[8],p3[8],p4[8];	
	double p[32];
	double mymax,mymid,mymin;
/*
	int alpha_GE_beta  = alpha >= beta;
	int beta_GE_gamma  = beta >= gamma;
	int alpha_GE_gamma = alpha >= gamma;	
	int odd = ODD(id[2]); 
	
	int i =  (odd << 3) + (alpha_GE_beta << 2) + (beta_GE_gamma << 1) + alpha_GE_gamma;
	
	id[0] -= odd;
	id[1] -= odd;
	
	id[0] >>= 1;
	id[1] >>= 1;
		
	switch( i ){
		case 0:
			mymax = gamma; mymid = beta; mymin = alpha;
				RegionSix_FillData_Even;break;
		case 2:
			mymax = beta; mymid = gamma; mymin = alpha;
				RegionFive_FillData_Even;break;
		case 3:
			mymax = beta; mymid = alpha; mymin = gamma;
				RegionFour_FillData_Even;break;
		case 4:
			mymax = gamma; mymid = alpha; mymin = beta;
				RegionThree_FillData_Even;break;
		case 5:
			mymax = alpha; mymid = gamma; mymin = beta;
				RegionTwo_FillData_Even;break;
		case 7:
			mymax = alpha; mymid = beta; mymin = gamma;
				RegionOne_FillData_Even;break;
		
		case 8:
			mymax = gamma; mymid = beta; mymin = alpha;
				RegionSix_FillData_Odd;break;
		case 10:
			mymax = beta; mymid = gamma; mymin = alpha;
				RegionFive_FillData_Odd;break;
		case 11:
			mymax = beta; mymid = alpha; mymin = gamma;
				RegionFour_FillData_Odd;break;
		case 12:
			mymax = gamma; mymid = alpha; mymin = beta;
				RegionThree_FillData_Odd;break;
		case 13:
			mymax = alpha; mymid = gamma; mymin = beta;
				RegionTwo_FillData_Odd;break;
				
		case 15:
			mymax = alpha; mymid = beta; mymin = gamma;
				RegionOne_FillData_Odd;
	}
*/
	//	id[0] = (alpha0+beta0-gamma0);id[1]=(alpha0-beta0+gamma0);id[2]=(-alpha0+beta0+gamma0);
			
	if(ODD(id[2]))
	{
		--id[0]>>=1; --id[1]>>=1;
		if(alpha >= beta)
		{
			if(beta >= gamma)
			{
				mymax = alpha; mymid = beta; mymin = gamma;
				RegionOne_FillData_Odd;
			}
			else if(alpha >= gamma)
			{
				mymax = alpha; mymid = gamma; mymin = beta;
				RegionTwo_FillData_Odd;
			}
			else
			{
				mymax = gamma; mymid = alpha; mymin = beta;
				RegionThree_FillData_Odd;
			}
		}
		else
		{
			if(alpha >= gamma)
			{
				mymax = beta; mymid = alpha; mymin = gamma;
				RegionFour_FillData_Odd;
			}
			else if(beta >= gamma)
			{
				mymax = beta; mymid = gamma; mymin = alpha;
				RegionFive_FillData_Odd;
			}
			else
			{
				mymax = gamma; mymid = beta; mymin = alpha;
				RegionSix_FillData_Odd;
			}
		}
	}
	else
	{
		id[0]>>=1; id[1]>>=1;
		if(alpha >= beta)
		{
			if(beta >= gamma)
			{
				mymax = alpha; mymid = beta; mymin = gamma;
				RegionOne_FillData_Even;
			}
			else if(alpha >= gamma)
			{
				mymax = alpha; mymid = gamma; mymin = beta;
				RegionTwo_FillData_Even;
			}
			else
			{
				mymax = gamma; mymid = alpha; mymin = beta;
				RegionThree_FillData_Even;
			}
		}
		else
		{
			if(alpha >= gamma)
			{
				mymax = beta; mymid = alpha; mymin = gamma;
				RegionFour_FillData_Even;
			}
			else if(beta >= gamma)
			{
				mymax = beta; mymid = gamma; mymin = alpha;
				RegionFive_FillData_Even;
			}
			else
			{
				mymax = gamma; mymid = beta; mymin = alpha;
				RegionSix_FillData_Even;
			}
		}
	}

// // 	register double p123;
// //     register double p0;
// //     register double alpha2;
// //     register double alpha3;
// //     register double base_rho;
// //     register double base_rho2;
  	register double result = 0.0;
	
// 	alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f;
// 	CONV_PPIPED_OLD(result,p[0],alpha,beta,gamma);
// 	alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f;
// 	CONV_PPIPED_OLD(result,p[1],alpha,beta,gamma);
// 	alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax);
// 	CONV_PPIPED_OLD(result,p[2],alpha,beta,gamma);
// 	alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f);
// 	CONV_PPIPED_OLD(result,p[3],alpha,beta,gamma);

//     double *pdata = &p[0][0];
  	CONV_PPIPED_OLD_INPLACE(result,p);

   
	return result;		

}

// double CBCCData::InterpolateQuinticBoxSplineTurbo( const vector3_t pos ) const
// {
// 	int id[3];
// 	int startX,startY,startZ;
// 	
// 	double x = pos[0], y = pos[1], z = pos[2];
// 	double alpha = .5*(x+y), beta = .5*(x+z), gamma = .5*(y+z);
// 	int alpha0 = (int)floor(alpha); alpha = alpha - alpha0;
// 	int beta0 = (int)floor(beta); beta = beta - beta0;
// 	int gamma0 = (int)floor(gamma); gamma = gamma - gamma0;
// 	id[0] = (alpha0+beta0-gamma0);id[1]=(alpha0-beta0+gamma0);id[2]=(-alpha0+beta0+gamma0);
//   
// 	//  int x0 = (alpha0+beta0-gamma0),y0=(alpha0-beta0+gamma0),z0=(-alpha0+beta0+gamma0);
// 
// 	//double p1[8],p2[8],p3[8],p4[8];	
// 	//double p[4][8];
// 	
// 	register double pdata[32];
// 	
// 	double mymax,mymid,mymin;
// 
// 	int alpha_GE_beta  = alpha >= beta;
// 	int beta_GE_gamma  = beta >= gamma;
// 	int alpha_GE_gamma = alpha >= gamma;	
// 	int odd = ODD(id[2]); 
// 	
// 	int i =  (odd << 3) + (alpha_GE_beta << 2) + (beta_GE_gamma << 1) + alpha_GE_gamma;
// 	
// 	id[0] -= odd;
// 	id[1] -= odd;
// 	
// 	id[0] >>= 1;
// 	id[1] >>= 1;
// 	
// 	startX = id[0];
// 	startY = id[1];
// 	startZ = id[2];
// 	
// 	double result = 0.0;
// 	
// 	switch( i ){
// 		case 0:
// 			mymax = gamma; mymid = beta; mymin = alpha;
// 				//RegionSix_FillData_Even;
// 				CONV_PPIPED_OLD_INPLACE_RegionSix_Even(result);
// 				return result;
// 				break;
// 		case 2:
// 			mymax = beta; mymid = gamma; mymin = alpha;
// 				//RegionFive_FillData_Even;
// 				CONV_PPIPED_OLD_INPLACE_RegionFive_Even(result);
// 				return result;
// 				break;
// 		case 3:
// 			mymax = beta; mymid = alpha; mymin = gamma;
// 				//RegionFour_FillData_Even;
// 				CONV_PPIPED_OLD_INPLACE_RegionFour_Even(result);
// 				return result;
// 				break;
// 		case 4:
// 			mymax = gamma; mymid = alpha; mymin = beta;
// 				//RegionThree_FillData_Even;
// 				CONV_PPIPED_OLD_INPLACE_RegionThree_Even(result);
// 				return result;
// 				break;
// 		case 5:
// 			mymax = alpha; mymid = gamma; mymin = beta;
// 				//RegionTwo_FillData_Even;
// 				CONV_PPIPED_OLD_INPLACE_RegionTwo_Even(result);
// 				return result;
// 				break;
// 		case 7:
// 			mymax = alpha; mymid = beta; mymin = gamma;
// 				//RegionOne_FillData_Even;
// 				CONV_PPIPED_OLD_INPLACE_RegionOne_Even(result);
// 				return result;
// 				break;
// 		
// 		case 8:
// 			mymax = gamma; mymid = beta; mymin = alpha;
// 				//RegionSix_FillData_Odd;
// 				CONV_PPIPED_OLD_INPLACE_RegionSix_Odd(result);
// 				return result;
// 				break;
// 		case 10:
// 			mymax = beta; mymid = gamma; mymin = alpha;
// 				//RegionFive_FillData_Odd;
// 				CONV_PPIPED_OLD_INPLACE_RegionFive_Odd(result);
// 				return result;
// 				break;
// 		case 11:
// 			mymax = beta; mymid = alpha; mymin = gamma;
// 				//RegionFour_FillData_Odd;
// 				CONV_PPIPED_OLD_INPLACE_RegionFour_Odd(result);
// 				return result;
// 				break;
// 		case 12:
// 			mymax = gamma; mymid = alpha; mymin = beta;
// 				//RegionThree_FillData_Odd;
// 				CONV_PPIPED_OLD_INPLACE_RegionThree_Odd(result);
// 				return result;
// 				break;
// 		case 13:
// 			mymax = alpha; mymid = gamma; mymin = beta;
// 				//RegionTwo_FillData_Odd;
// 				CONV_PPIPED_OLD_INPLACE_RegionTwo_Odd(result);
// 				return result;
// 				break;
// 				
// 		case 15:
// 			mymax = alpha; mymid = beta; mymin = gamma;
// 				//RegionOne_FillData_Odd;
// 				CONV_PPIPED_OLD_INPLACE_RegionOne_Odd(result);
// 				return result;
// 	}
// 	
// 	return 0;
// 	
// }

void   CBCCData::InterpolateNormalQuinticBoxSplineTurbo( const vector3_t pos, vector3_t N ) const
{
					
	int id[3];

	double x = pos[0], y = pos[1], z = pos[2];
	register double alpha = .5f*(x+y), beta = .5f*(x+z), gamma = .5f*(y+z);
	//int alpha0 = (int)floor(alpha); alpha = alpha - alpha0;
	//int beta0 = (int)floor(beta); beta = beta - beta0;
	//int gamma0 = (int)floor(gamma); gamma = gamma - gamma0;
	int alpha0 = (int)(alpha); alpha = alpha - alpha0;
	int beta0 = (int)(beta); beta = beta - beta0;
	int gamma0 = (int)(gamma); gamma = gamma - gamma0;
	id[0] = (alpha0+beta0-gamma0);id[1]=(alpha0-beta0+gamma0);id[2]=(-alpha0+beta0+gamma0);
		
//  int x0 = (alpha0+beta0-gamma0),y0=(alpha0-beta0+gamma0),z0=(-alpha0+beta0+gamma0);

	double p1x[8],p2x[8],p3x[8],p4x[8];
	double p1y[8],p2y[8],p3y[8],p4y[8];
	double p1z[8],p2z[8],p3z[8],p4z[8];
	
	register double mymax,mymid,mymin;
/*
	int alpha_GE_beta  = alpha >= beta;
	int beta_GE_gamma  = beta >= gamma;
	int alpha_GE_gamma = alpha >= gamma;	
	int odd = ODD(id[2]); 
	
	int region =  (odd << 3) + (alpha_GE_beta << 2) + (beta_GE_gamma << 1) + alpha_GE_gamma;
	
	
	
	id[0] -= odd;
	id[1] -= odd;
	
	id[0] >>= 1;
	id[1] >>= 1;
		
	switch( region ){
		case 0:			
			mymax = gamma; mymid = beta; mymin = alpha;
				Norm_RegionSix_FillData_Even;break;
		case 2:
			mymax = beta; mymid = gamma; mymin = alpha;
				Norm_RegionFive_FillData_Even;break;
		case 3:
			mymax = beta; mymid = alpha; mymin = gamma;
				Norm_RegionFour_FillData_Even;break;
		case 4:
			mymax = gamma; mymid = alpha; mymin = beta;
				Norm_RegionThree_FillData_Even;break;
		case 5:
			mymax = alpha; mymid = gamma; mymin = beta;
				Norm_RegionTwo_FillData_Even;break;
		case 7:
			mymax = alpha; mymid = beta; mymin = gamma;
				Norm_RegionOne_FillData_Even;break;
		
		case 8:
			mymax = gamma; mymid = beta; mymin = alpha;
				Norm_RegionSix_FillData_Odd;break;
		case 10:
			mymax = beta; mymid = gamma; mymin = alpha;
				Norm_RegionFive_FillData_Odd;break;
		case 11:
			mymax = beta; mymid = alpha; mymin = gamma;
				Norm_RegionFour_FillData_Odd;break;
		case 12:
			mymax = gamma; mymid = alpha; mymin = beta;
				Norm_RegionThree_FillData_Odd;break;
		case 13:
			mymax = alpha; mymid = gamma; mymin = beta;
				Norm_RegionTwo_FillData_Odd;break;
				
		case 15:
			mymax = alpha; mymid = beta; mymin = gamma;
				Norm_RegionOne_FillData_Odd;
	}
		
	id[0] = (alpha0+beta0-gamma0);id[1]=(alpha0-beta0+gamma0);id[2]=(-alpha0+beta0+gamma0);
*/

	if(ODD(id[2]))
	{
		--id[0]>>=1; --id[1]>>=1;
		if(alpha >= beta)
		{
			if(beta >= gamma)
			{
				mymax = alpha; mymid = beta; mymin = gamma;
				Norm_RegionOne_FillData_Odd;
			}
			else if(alpha >= gamma)
			{
				mymax = alpha; mymid = gamma; mymin = beta;
				Norm_RegionTwo_FillData_Odd;
			}
			else
			{
				mymax = gamma; mymid = alpha; mymin = beta;
				Norm_RegionThree_FillData_Odd;
			}
		}
		else
		{
			if(alpha >= gamma)
			{
				mymax = beta; mymid = alpha; mymin = gamma;
				Norm_RegionFour_FillData_Odd;
			}
			else if(beta >= gamma)
			{
				mymax = beta; mymid = gamma; mymin = alpha;
				Norm_RegionFive_FillData_Odd;
			}
			else
			{
				mymax = gamma; mymid = beta; mymin = alpha;
				Norm_RegionSix_FillData_Odd;
			}
		}
	}
	else
	{
		id[0]>>=1; id[1]>>=1;
		if(alpha >= beta)
		{
			if(beta >= gamma)
			{
				mymax = alpha; mymid = beta; mymin = gamma;
				Norm_RegionOne_FillData_Even;
			}
			else if(alpha >= gamma)
			{
				mymax = alpha; mymid = gamma; mymin = beta;
				Norm_RegionTwo_FillData_Even;
			}
			else
			{
				mymax = gamma; mymid = alpha; mymin = beta;
				Norm_RegionThree_FillData_Even;
			}
		}
		else
		{
			if(alpha >= gamma)
			{
				mymax = beta; mymid = alpha; mymin = gamma;
				Norm_RegionFour_FillData_Even;
			}
			else if(beta >= gamma)
			{
				mymax = beta; mymid = gamma; mymin = alpha;
				Norm_RegionFive_FillData_Even;
			}
			else
			{
				mymax = gamma; mymid = beta; mymin = alpha;
				Norm_RegionSix_FillData_Even;
			}
		}
	}
	
	register double resultx = 0.0f;
	register double resulty = 0.0f;
	register double resultz = 0.0f;
	
	double alphax,alphay,alphaz;
	double betax,betay,betaz;
	double gammax,gammay,gammaz;
	
	alphax=mymax-1.0f;betax=mymid-1.0f;gammax=mymin-1.0f;
	alphaz=alphay=alphax;
	betaz=betay=betax;
	gammaz=gammay=gammax;
	CONV_PPIPED_OLD(resultx,p1x,alphax,betax,gammax);
	CONV_PPIPED_OLD(resulty,p1y,alphay,betay,gammay);
	CONV_PPIPED_OLD(resultz,p1z,alphaz,betaz,gammaz);
	
	alphax=-mymin;betax=mymax-mymin-1.0f;gammax=mymid-mymin-1.0f;
	alphaz=alphay=alphax;
	betaz=betay=betax;
	gammaz=gammay=gammax;
	CONV_PPIPED_OLD(resultx,p2x,alphax,betax,gammax);
	CONV_PPIPED_OLD(resulty,p2y,alphay,betay,gammay);
	CONV_PPIPED_OLD(resultz,p2z,alphaz,betaz,gammaz);
	
	
	alphax=(-mymax+mymid);betax=(-mymax+mymin);gammax=(-mymax);
	alphaz=alphay=alphax;
	betaz=betay=betax;
	gammaz=gammay=gammax;	
	CONV_PPIPED_OLD(resultx,p3x,alphax,betax,gammax);
	CONV_PPIPED_OLD(resulty,p3y,alphay,betay,gammay);
	CONV_PPIPED_OLD(resultz,p3z,alphaz,betaz,gammaz);
	
	
	alphax=(-mymid+mymin);betax=(-mymid);gammax=(mymax-mymid-1.0f);
	alphaz=alphay=alphax;
	betaz=betay=betax;
	gammaz=gammay=gammax;
	CONV_PPIPED_OLD(resultx,p4x,alphax,betax,gammax);
	CONV_PPIPED_OLD(resulty,p4y,alphay,betay,gammay);
	CONV_PPIPED_OLD(resultz,p4z,alphaz,betaz,gammaz);
	
	
	
	
	N[0] = resultx;
	N[1] = resulty;
	N[2] = resultz;
			

}


// void   CBCCData::InterpolateNormalQuinticBoxSplineTurbo( const vector3_t pos, vector3_t N ) const
// {
// 		
// 	for(int ordinate=0;ordinate<3;ordinate++){		
// 	
// 		int id[3];
// 
// 		double x = pos[0], y = pos[1], z = pos[2];
// 		register double alpha = .5f*(x+y), beta = .5f*(x+z), gamma = .5f*(y+z);
// 		int alpha0 = (int)floor(alpha); alpha = alpha - alpha0;
// 		int beta0 = (int)floor(beta); beta = beta - beta0;
// 		int gamma0 = (int)floor(gamma); gamma = gamma - gamma0;
// 	
// 	  
// 	//  int x0 = (alpha0+beta0-gamma0),y0=(alpha0-beta0+gamma0),z0=(-alpha0+beta0+gamma0);
// 
// 		double p1[8],p2[8],p3[8],p4[8];
// 		register double mymax,mymid,mymin;
// 
// 	
// 		
// 		id[0] = (alpha0+beta0-gamma0);id[1]=(alpha0-beta0+gamma0);id[2]=(-alpha0+beta0+gamma0);
// 			
// 		if(ODD(id[2]))
// 		{
// 			--id[0]>>=1; --id[1]>>=1;
// 			if(alpha >= beta)
// 			{
// 				if(beta >= gamma)
// 				{
// 					mymax = alpha; mymid = beta; mymin = gamma;
// 					Norm_RegionOne_FillData_Odd(ordinate);
// 				}
// 				else if(alpha >= gamma)
// 				{
// 					mymax = alpha; mymid = gamma; mymin = beta;
// 					Norm_RegionTwo_FillData_Odd(ordinate);
// 				}
// 				else
// 				{
// 					mymax = gamma; mymid = alpha; mymin = beta;
// 					Norm_RegionThree_FillData_Odd(ordinate);
// 				}
// 			}
// 			else
// 			{
// 				if(alpha >= gamma)
// 				{
// 					mymax = beta; mymid = alpha; mymin = gamma;
// 					Norm_RegionFour_FillData_Odd(ordinate);
// 				}
// 				else if(beta >= gamma)
// 				{
// 					mymax = beta; mymid = gamma; mymin = alpha;
// 					Norm_RegionFive_FillData_Odd(ordinate);
// 				}
// 				else
// 				{
// 					mymax = gamma; mymid = beta; mymin = alpha;
// 					Norm_RegionSix_FillData_Odd(ordinate);
// 				}
// 			}
// 		}
// 		else
// 		{
// 			id[0]>>=1; id[1]>>=1;
// 			if(alpha >= beta)
// 			{
// 				if(beta >= gamma)
// 				{
// 					mymax = alpha; mymid = beta; mymin = gamma;
// 					Norm_RegionOne_FillData_Even(ordinate);
// 				}
// 				else if(alpha >= gamma)
// 				{
// 					mymax = alpha; mymid = gamma; mymin = beta;
// 					Norm_RegionTwo_FillData_Even(ordinate);
// 				}
// 				else
// 				{
// 					mymax = gamma; mymid = alpha; mymin = beta;
// 					Norm_RegionThree_FillData_Even(ordinate);
// 				}
// 			}
// 			else
// 			{
// 				if(alpha >= gamma)
// 				{
// 					mymax = beta; mymid = alpha; mymin = gamma;
// 					Norm_RegionFour_FillData_Even(ordinate);
// 				}
// 				else if(beta >= gamma)
// 				{
// 					mymax = beta; mymid = gamma; mymin = alpha;
// 					Norm_RegionFive_FillData_Even(ordinate);
// 				}
// 				else
// 				{
// 					mymax = gamma; mymid = beta; mymin = alpha;
// 					Norm_RegionSix_FillData_Even(ordinate);
// 				}
// 			}
// 		}
// 	
// 		register double result = 0.0f;
// 		alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f;
// 		CONV_PPIPED(result,p1,alpha,beta,gamma);
// 		alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f;
// 		CONV_PPIPED(result,p2,alpha,beta,gamma);
// 		alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax);
// 		CONV_PPIPED(result,p3,alpha,beta,gamma);
// 		alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f);
// 		CONV_PPIPED(result,p4,alpha,beta,gamma);
// 	
// 	
// 		N[ordinate] = result;
// 	}
// }


void   CBCCData::InterpolateNormalLinearBoxSplineTurbo( const vector3_t pos, vector3_t N ) const
{	
	for(int ordinate = 0; ordinate<3;ordinate++){
		int id[3];
	
		double x = pos[0], y = pos[1], z = pos[2];
		register double alpha = .5f*(x+y), beta = .5f*(x+z), gamma = .5f*(y+z);
		int alpha0 = (int)floor(alpha); alpha = alpha - alpha0;
		int beta0 = (int)floor(beta); beta = beta - beta0;
		int gamma0 = (int)floor(gamma); gamma = gamma - gamma0;
		id[0] = (alpha0+beta0-gamma0);id[1]=(alpha0-beta0+gamma0);id[2]=(-alpha0+beta0+gamma0);
	//  int x0 = (alpha0+beta0-gamma0),y0=(alpha0-beta0+gamma0),z0=(-alpha0+beta0+gamma0);
	
	
		double p1,p2,p3,p4;
		register double mymax,mymid,mymin;
	
	
		if(ODD(id[2]))
		{
			--id[0]>>=1; --id[1]>>=1;
			if(alpha >= beta)
			{
				if(beta >= gamma)
				{
					mymax = alpha; mymid = beta; mymin = gamma;
					p1=GetValueNorm(id,ordinate);
					id[0]++;        p4=GetValueNorm(id,ordinate);
					id[1]++;id[2]++;        p2=GetValueNorm(id,ordinate);
					id[2]-=2;       p3=GetValueNorm(id,ordinate);
				}
				else if(alpha >= gamma)
				{
					mymax = alpha; mymid = gamma; mymin = beta;
					p1=GetValueNorm(id,ordinate);
					id[1]++;        p4=GetValueNorm(id,ordinate);
					id[0]++;id[2]++;        p2=GetValueNorm(id,ordinate);
					id[2]-=2;       p3=GetValueNorm(id,ordinate);
				}
				else
				{
					mymax = gamma; mymid = alpha; mymin = beta;
					p1=GetValueNorm(id,ordinate);
					id[1]++;        p4=GetValueNorm(id,ordinate);
					id[2]++;        p3=GetValueNorm(id,ordinate);
					id[0]++;        p2=GetValueNorm(id,ordinate);
				}
			}
			else
			{
				if(alpha >= gamma)
				{
					mymax = beta; mymid = alpha; mymin = gamma;
					p1=GetValueNorm(id,ordinate);
					id[0]++;        p4=GetValueNorm(id,ordinate);
					id[2]++;        p3=GetValueNorm(id,ordinate);
					id[1]++;        p2=GetValueNorm(id,ordinate);
				}
				else if(beta >= gamma)
				{
					mymax = beta; mymid = gamma; mymin = alpha;
					p1=GetValueNorm(id,ordinate);
					id[2]+=2;        p4=GetValueNorm(id,ordinate);
					id[0]++;id[2]--;        p3=GetValueNorm(id,ordinate);
					id[1]++;        p2=GetValueNorm(id,ordinate);
		
				}
				else
				{
					mymax = gamma; mymid = beta; mymin = alpha;
					p1=GetValueNorm(id,ordinate);
					id[2]+=2;        p4=GetValueNorm(id,ordinate);
					id[1]++;id[2]--;        p3=GetValueNorm(id,ordinate);
					id[0]++;        p2=GetValueNorm(id,ordinate);
	
				}
			}
		}
		else
		{
			id[0]>>=1; id[1]>>=1;
			if(alpha >= beta)
			{
				if(beta >= gamma)
				{
					mymax = alpha; mymid = beta; mymin = gamma;
					p1=GetValueNorm(id,ordinate);
					id[0]++;        p4=GetValueNorm(id,ordinate);
					id[0]--;id[2]++;        p2=GetValueNorm(id,ordinate);
					id[2]-=2;       p3=GetValueNorm(id,ordinate);
				}
				else if(alpha >= gamma)
				{
					mymax = alpha; mymid = gamma; mymin = beta;
					p1=GetValueNorm(id,ordinate);
					id[1]++;        p4=GetValueNorm(id,ordinate);
					id[1]--;id[2]++;        p2=GetValueNorm(id,ordinate);
					id[2]-=2;       p3=GetValueNorm(id,ordinate);
				}
				else
				{
					mymax = gamma; mymid = alpha; mymin = beta;
					p1=GetValueNorm(id,ordinate);
					id[1]++;        p4=GetValueNorm(id,ordinate);
					id[1]--;id[2]++;        p2=GetValueNorm(id,ordinate);
					id[0]--;        p3=GetValueNorm(id,ordinate);
				}
			}
			else
			{
				if(alpha >= gamma)
				{
					mymax = beta; mymid = alpha; mymin = gamma;
					p1=GetValueNorm(id,ordinate);
					id[0]++;        p4=GetValueNorm(id,ordinate);
					id[0]--;id[2]++;        p2=GetValueNorm(id,ordinate);
					id[1]--;        p3=GetValueNorm(id,ordinate);
				}
				else if(beta >= gamma)
				{
					mymax = beta; mymid = gamma; mymin = alpha;
					p1=GetValueNorm(id,ordinate);
					id[2]++;        p2=GetValueNorm(id,ordinate);
					id[1]--;        p3=GetValueNorm(id,ordinate);
					id[1]++;id[2]++;        p4=GetValueNorm(id,ordinate);
				}
				else
				{
					mymax = gamma; mymid = beta; mymin = alpha;
					p1=GetValueNorm(id,ordinate);
					id[2]++;        p2=GetValueNorm(id,ordinate);
					id[0]--;        p3=GetValueNorm(id,ordinate);
					id[0]++;id[2]++;        p4=GetValueNorm(id,ordinate);
				}
			}
		}
	
		N[ordinate] = p1 + mymax * (p3 - p1) + mymin * (p2 - p4) + mymid * (p4 - p3);
	}
}
	
double CBCCData::Interpolate( const vector3_t bcc ) const
{
	int index3D[3];

	//int gridId = 0;

	if( m_reconFilter == BCC_LINEAR_BOX_SPLINE ){
// #ifdef LINEAR

		// linear box spline

		double alpha = .5*(bcc[X]+bcc[Y]), beta = .5*(bcc[X]+bcc[Z]), gamma = .5*(bcc[Y]+bcc[Z]);
		int alpha0 = (int)floor(alpha); alpha = alpha - alpha0;
		int beta0 = (int)floor(beta); beta = beta - beta0;
		int gamma0 = (int)floor(gamma); gamma = gamma - gamma0;
		int x0 =
		(alpha0+beta0-gamma0),y0=(alpha0-beta0+gamma0),z0=(-alpha0+beta0+gamma0);

		//double kern = 0.0, accumKernel = 0.0;

		BCCToIndex(index3D,x0,y0,z0);
		double p1 = GetValue( index3D );// * kern; accumKernel += kern;
		BCCToIndex(index3D,x0+1,y0+1,z0+1);

		double p2 = GetValue( index3D );// * kern; accumKernel += kern;

		double p3, p4;

		double mymax,mymid,mymin;

		if(alpha >= beta)
		{
		  if(beta >= gamma)
			{
			  mymax = alpha; mymid = beta; mymin = gamma;
			  BCCToIndex(index3D,x0+1,y0+1,z0-1);
			  p3 = GetValue( index3D );// * kern; accumKernel += kern;
			  BCCToIndex(index3D,x0+2,y0,z0);
			  p4 = GetValue( index3D );// * kern; accumKernel += kern;
			}
		  else if(alpha >= gamma)
			{
			  mymax = alpha; mymid = gamma; mymin = beta;
			  BCCToIndex(index3D,x0+1,y0+1,z0-1);
			  p3 = GetValue( index3D );// * kern; accumKernel += kern;
			  BCCToIndex(index3D,x0,y0+2,z0);
			  p4 = GetValue( index3D );// * kern; accumKernel += kern;
			}
		  else
			{
			  mymax = gamma; mymid = alpha; mymin = beta;
			  BCCToIndex(index3D,x0-1,y0+1,z0+1);
			  p3 = GetValue( index3D );// * kern; accumKernel += kern;
			  BCCToIndex(index3D,x0,y0+2,z0);
			  p4 = GetValue( index3D );// * kern; accumKernel += kern;
			}
		}
		else
		{
		  if(alpha >= gamma)
			{
			  mymax = beta; mymid = alpha; mymin = gamma;
			  BCCToIndex(index3D,x0+1,y0-1,z0+1);
			  p3 = GetValue( index3D );// * kern; accumKernel += kern;
			  BCCToIndex(index3D,x0+2,y0,z0);
			  p4 = GetValue( index3D );// * kern; accumKernel += kern;
		}
		  else if(beta >= gamma)
			{
			  mymax = beta; mymid = gamma; mymin = alpha;
			  BCCToIndex(index3D,x0+1,y0-1,z0+1);
			  p3 = GetValue( index3D );// * kern; accumKernel += kern;
			  BCCToIndex(index3D,x0,y0,z0+2);
			  p4 = GetValue( index3D );// * kern; accumKernel += kern;
			}
		  else
			{
			  mymax = gamma; mymid = beta; mymin = alpha;
			  BCCToIndex(index3D,x0-1,y0+1,z0+1);
			  p3 = GetValue( index3D );// * kern; accumKernel += kern;
			  BCCToIndex(index3D,x0,y0,z0+2);
			  p4 = GetValue( index3D );// * kern; accumKernel += kern;
			}
		}

		return (double)(p1 + mymax * (p3 - p1) + mymin * (p2 - p4) + mymid * (p4 - p3));
	}

	else // QUINTIC_BOX_SPLINE
	{
 //#else
 
		register double alpha = .5f*(bcc[X]+bcc[Y]), beta = .5f*(bcc[X]+bcc[Z]), gamma = .5f*(bcc[Y]+bcc[Z]);

		int alpha0 = (int)floor(alpha); alpha = alpha - alpha0;
		int beta0 = (int)floor(beta); beta = beta - beta0;
		int gamma0 = (int)floor(gamma); gamma = gamma - gamma0;
		int x0 =
		(alpha0+beta0-gamma0),y0=(alpha0-beta0+gamma0),z0=(-alpha0+beta0+gamma0);

		double p1[8],p2[8],p3[8],p4[8];
		register double result = 0.0f;


		register double mymax,mymid,mymin;

		if(alpha >= beta)
		{
		  if(beta >= gamma)
			{
			  mymax = alpha; mymid = beta; mymin = gamma;
			  FILL_PPIPED(p1,-1,-1,-1,1,2,3);
			  x0++; y0++;z0++;
			  FILL_PPIPED(p2,-1,1,1,3,2,1);
			  z0 -= 2;
			  FILL_PPIPED(p3,1,1,-1,3,2,1);
			  x0++;y0--;z0++;
			  FILL_PPIPED(p4,1,-1,1,1,2,3);

			}
		  else if(alpha >= gamma)
			{
			  mymax = alpha; mymid = gamma; mymin = beta;
			  FILL_PPIPED(p1,-1,-1,-1,2,1,3);
			  x0++; y0++;z0++;
			  FILL_PPIPED(p2,1,-1,1,2,3,1);
			  z0 -= 2;
			  FILL_PPIPED(p3,1,1,-1,2,3,1);
			  x0--;y0++;z0++;
			  FILL_PPIPED(p4,-1,1,1,2,1,3);
			}
		  else
			{
			  mymax = gamma; mymid = alpha; mymin = beta;
			  FILL_PPIPED(p1,-1,-1,-1,3,1,2);
			  x0++; y0++;z0++;
			  FILL_PPIPED(p2,1,-1,1,1,3,2);
			  x0 -= 2;
			  FILL_PPIPED(p3,-1,1,1,1,3,2);
			  x0++;y0++;z0--;
			  FILL_PPIPED(p4,1,1,-1,3,1,2);
			}
		}
		else
		{
		  if(alpha >= gamma)
			{
			  mymax = beta; mymid = alpha; mymin = gamma;
			  FILL_PPIPED(p1,-1,-1,-1,1,3,2);
			  x0++; y0++; z0++;
			  FILL_PPIPED(p2,-1,1,1,3,1,2);
			  y0 -= 2;
			  FILL_PPIPED(p3,1,-1,1,3,1,2);
			  x0++; y0++; z0--;
			  FILL_PPIPED(p4,1,1,-1,1,3,2);
			}
		  else if(beta >= gamma)
			{
			  mymax = beta; mymid = gamma; mymin = alpha;
			  FILL_PPIPED(p1,-1,-1,-1,2,3,1);
			  x0++; y0++; z0++;
			  FILL_PPIPED(p2,1,1,-1,2,1,3);
			  y0 -= 2;
			  FILL_PPIPED(p3,1,-1,1,2,1,3);
			  x0--;y0++;z0++;
			  FILL_PPIPED(p4,-1,1,1,2,3,1);
			}
		  else
			{
			  mymax = gamma; mymid = beta; mymin = alpha;
			  FILL_PPIPED(p1,-1,-1,-1,3,2,1);
			  x0++; y0++; z0++;
			  FILL_PPIPED(p2,1,1,-1,1,2,3);
			  x0 -= 2;
			  FILL_PPIPED(p3,-1,1,1,1,2,3);
			  x0++;y0--;z0++;
			  FILL_PPIPED(p4,1,-1,1,3,2,1);
			}
		}


		alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f;
		CONV_PPIPED(result,p1,alpha,beta,gamma);
		alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f;
		CONV_PPIPED(result,p2,alpha,beta,gamma);
		alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax);
		CONV_PPIPED(result,p3,alpha,beta,gamma);
		alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f);
		CONV_PPIPED(result,p4,alpha,beta,gamma);

		return (double)(result);
 }

//#endif

}

//double CBCCData::Get( const vector3_t p ) const


void CBCCData::Evaluate( const vector3_t min, const vector3_t max, const vector3_t step,
		elem_t (*fn)(const vector3_t p), void (*fnNormal)(const vector3_t p, vector3_t N) )
{

	DestroyData();

	m_xSize = (int)((max[X] - min[X]) / step[X])+1;
	m_ySize = (int)((max[Y] - min[Y]) / step[Y])+1;
	m_zSize = (int)((max[Z] - min[Z]) / step[Z])+1;

	m_halfScale[X] = step[X] / 2.0f;
	m_halfScale[Y] = step[Y] / 2.0f;
	m_halfScale[Z] = step[Z] / 1.0f;

	m_invHalfScale[X] = 1.0/m_halfScale[X];
	m_invHalfScale[Y] = 1.0/m_halfScale[Y];
	m_invHalfScale[Z] = 1.0/m_halfScale[Z];


	m_pppData = new elem_t** [ m_zSize + 2*GUARD_BAND ];
	for( int z = 0; z < (m_zSize+2*GUARD_BAND); z++ ){
		m_pppData[z] = new  elem_t *  [ m_ySize + 2*GUARD_BAND ];
		for( int y = 0; y< (m_ySize+2*GUARD_BAND); y++ ){
			m_pppData[z][y] = new elem_t [ m_xSize + 2*GUARD_BAND ];
			memset( m_pppData[z][y], 0, sizeof( elem_t) * (m_xSize + 2*GUARD_BAND) );
		}
	}



	//m_pVolume = new CVolume( m_xSize, m_ySize, m_zSize, step[X], step[Y], step[Z] );

	// Nasty hack :(
	m_pVolume = new CVolume( (m_xSize - 0.5) * step[X], (m_ySize - 0.5) * step[Y], (m_zSize - 1) * step[Z], step[X], step[Y], step[Z] );

	m_translate[X] = m_pVolume->m_HXSize / m_halfScale[X];
	m_translate[Y] = m_pVolume->m_HYSize / m_halfScale[Y];
	m_translate[Z] = m_pVolume->m_HZSize / m_halfScale[Z];


	for(int z=0;z<m_zSize;++z){
		for(int y=0;y<m_ySize;++y){
			for(int x=0;x<m_xSize;++x){

				vector3_t p;

				p[X] = x*step[X] + min[X];
				p[Y] = y*step[Y] + min[Y];

				if( z % 2 == 1){
					p[X] += step[X]/2.0f;
					p[Y] += step[Y]/2.0f;
				}

				p[Z] = z*step[Z] + min[Z];

				m_pppData[z+GUARD_BAND][y+GUARD_BAND][x+GUARD_BAND] = fn( p );



			}
		}
	}

	// Creating normals memory

	m_pppNormals = new vector3_t** [ m_zSize + 2*GUARD_BAND];

	for(int z = 0; z< (m_zSize+2*GUARD_BAND);++z){

		m_pppNormals[z] = new vector3_t* [ m_ySize + 2*GUARD_BAND];

		for(int y=0;y<(m_ySize+2*GUARD_BAND);++y){
			m_pppNormals[z][y] = new vector3_t [ m_xSize + 2*GUARD_BAND];
			memset(m_pppNormals[z][y], 0, sizeof(vector3_t) * (m_xSize+2*GUARD_BAND) );
		}
	}

	CalculateNormals( fnNormal );

	
}
/*
void  CBCCData::GetNormal( const vector3_t p, vector3_t N ) const
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
*/

		/*		
{
	if( m_bEpsilonNormalMode ){
		CBaseData::GetNormal( p, N );
		return;
	}

  int index3D[3];
  vector3_t bcc;
  WorldToBCC(bcc,p);

  for(int ordinate=0;ordinate<3;ordinate++){
		//#ifdef LINEAR
		if( m_normalReconFilter == BCC_LINEAR_BOX_SPLINE ){
			// linear box spline

		  double alpha = .5*(bcc[X]+bcc[Y]), beta = .5*(bcc[X]+bcc[Z]), gamma = .5*(bcc[Y]+bcc[Z]);
		  int alpha0 = (int)floor(alpha); alpha = alpha - alpha0;
		  int beta0 = (int)floor(beta); beta = beta - beta0;
		  int gamma0 = (int)floor(gamma); gamma = gamma - gamma0;
		  int x0 =
			(alpha0+beta0-gamma0),y0=(alpha0-beta0+gamma0),z0=(-alpha0+beta0+gamma0);

		  //double kern = 0.0, accumKernel = 0.0;

		  BCCToIndex(index3D,x0,y0,z0);
		  double p1 = GetValueNorm(index3D, ordinate);// * kern; accumKernel += kern;
		  BCCToIndex(index3D,x0+1,y0+1,z0+1);

		  double p2 = GetValueNorm(index3D, ordinate);// * kern; accumKernel += kern;

		  double p3, p4;

		  double mymax,mymid,mymin;

		  if(alpha >= beta)
			{
			  if(beta >= gamma)
				{
				  mymax = alpha; mymid = beta; mymin = gamma;
				  BCCToIndex(index3D,x0+1,y0+1,z0-1);
				  p3 = GetValueNorm(index3D, ordinate);// * kern; accumKernel += kern;
				  BCCToIndex(index3D,x0+2,y0,z0);
				  p4 = GetValueNorm(index3D, ordinate);// * kern; accumKernel += kern;
				}
			  else if(alpha >= gamma)
				{
				  mymax = alpha; mymid = gamma; mymin = beta;
				  BCCToIndex(index3D,x0+1,y0+1,z0-1);
				  p3 = GetValueNorm(index3D, ordinate);// * kern; accumKernel += kern;
				  BCCToIndex(index3D,x0,y0+2,z0);
				  p4 = GetValueNorm(index3D, ordinate);// * kern; accumKernel += kern;
				}
			  else
				{
				  mymax = gamma; mymid = alpha; mymin = beta;
				  BCCToIndex(index3D,x0-1,y0+1,z0+1);
				  p3 = GetValueNorm(index3D, ordinate);// * kern; accumKernel += kern;
				  BCCToIndex(index3D,x0,y0+2,z0);
				  p4 = GetValueNorm(index3D, ordinate);// * kern; accumKernel += kern;
				}
			}
		  else
			{
			  if(alpha >= gamma)
				{
				  mymax = beta; mymid = alpha; mymin = gamma;
				  BCCToIndex(index3D,x0+1,y0-1,z0+1);
				  p3 = GetValueNorm(index3D, ordinate);// * kern; accumKernel += kern;
				  BCCToIndex(index3D,x0+2,y0,z0);
				  p4 = GetValueNorm(index3D, ordinate);// * kern; accumKernel += kern;
			}
			  else if(beta >= gamma)
				{
				  mymax = beta; mymid = gamma; mymin = alpha;
				  BCCToIndex(index3D,x0+1,y0-1,z0+1);
				  p3 = GetValueNorm(index3D, ordinate);// * kern; accumKernel += kern;
				  BCCToIndex(index3D,x0,y0,z0+2);
				  p4 = GetValueNorm(index3D, ordinate);// * kern; accumKernel += kern;
				}
			  else
				{
				  mymax = gamma; mymid = beta; mymin = alpha;
				  BCCToIndex(index3D,x0-1,y0+1,z0+1);
				  p3 = GetValueNorm(index3D, ordinate);// * kern; accumKernel += kern;
				  BCCToIndex(index3D,x0,y0,z0+2);
				  p4 = GetValueNorm(index3D, ordinate);// * kern; accumKernel += kern;
				}
			}

		   N[ordinate] = elem_t(p1 + mymax * (p3 - p1) + mymin * (p2 - p4) + mymid * (p4 - p3));


		}
		else{ // Quintic Box Spline
		//#else

		  register double alpha = .5f*(bcc[X]+bcc[Y]), beta = .5f*(bcc[X]+bcc[Z]), gamma = .5f*(bcc[Y]+bcc[Z]);

		  int alpha0 = (int)floor(alpha); alpha = alpha - alpha0;
		  int beta0 = (int)floor(beta); beta = beta - beta0;
		  int gamma0 = (int)floor(gamma); gamma = gamma - gamma0;
		  int x0 =
			(alpha0+beta0-gamma0),y0=(alpha0-beta0+gamma0),z0=(-alpha0+beta0+gamma0);

		  double p1[8],p2[8],p3[8],p4[8];
		  register double result = 0.0f;


		  register double mymax,mymid,mymin;

		  if(alpha >= beta)
			{
			  if(beta >= gamma)
				{
				  mymax = alpha; mymid = beta; mymin = gamma;
				  FILL_PPIPED_NORM(p1,-1,-1,-1,1,2,3,ordinate);
				  x0++; y0++;z0++;
				  FILL_PPIPED_NORM(p2,-1,1,1,3,2,1,ordinate);
				  z0 -= 2;
				  FILL_PPIPED_NORM(p3,1,1,-1,3,2,1,ordinate);
				  x0++;y0--;z0++;
				  FILL_PPIPED_NORM(p4,1,-1,1,1,2,3,ordinate);

				}
			  else if(alpha >= gamma)
				{
				  mymax = alpha; mymid = gamma; mymin = beta;
				  FILL_PPIPED_NORM(p1,-1,-1,-1,2,1,3,ordinate);
				  x0++; y0++;z0++;
				  FILL_PPIPED_NORM(p2,1,-1,1,2,3,1,ordinate);
				  z0 -= 2;
				  FILL_PPIPED_NORM(p3,1,1,-1,2,3,1,ordinate);
				  x0--;y0++;z0++;
				  FILL_PPIPED_NORM(p4,-1,1,1,2,1,3,ordinate);
				}
			  else
				{
				  mymax = gamma; mymid = alpha; mymin = beta;
				  FILL_PPIPED_NORM(p1,-1,-1,-1,3,1,2,ordinate);
				  x0++; y0++;z0++;
				  FILL_PPIPED_NORM(p2,1,-1,1,1,3,2,ordinate);
				  x0 -= 2;
				  FILL_PPIPED_NORM(p3,-1,1,1,1,3,2,ordinate);
				  x0++;y0++;z0--;
				  FILL_PPIPED_NORM(p4,1,1,-1,3,1,2,ordinate);
				}
			}
		  else
			{
			  if(alpha >= gamma)
				{
				  mymax = beta; mymid = alpha; mymin = gamma;
				  FILL_PPIPED_NORM(p1,-1,-1,-1,1,3,2,ordinate);
				  x0++; y0++; z0++;
				  FILL_PPIPED_NORM(p2,-1,1,1,3,1,2,ordinate);
				  y0 -= 2;
				  FILL_PPIPED_NORM(p3,1,-1,1,3,1,2,ordinate);
				  x0++; y0++; z0--;
				  FILL_PPIPED_NORM(p4,1,1,-1,1,3,2,ordinate);
				}
			  else if(beta >= gamma)
				{
				  mymax = beta; mymid = gamma; mymin = alpha;
				  FILL_PPIPED_NORM(p1,-1,-1,-1,2,3,1,ordinate);
				  x0++; y0++; z0++;
				  FILL_PPIPED_NORM(p2,1,1,-1,2,1,3,ordinate);
				  y0 -= 2;
				  FILL_PPIPED_NORM(p3,1,-1,1,2,1,3,ordinate);
				  x0--;y0++;z0++;
				  FILL_PPIPED_NORM(p4,-1,1,1,2,3,1,ordinate);
				}
			  else
				{
				  mymax = gamma; mymid = beta; mymin = alpha;
				  FILL_PPIPED_NORM(p1,-1,-1,-1,3,2,1,ordinate);
				  x0++; y0++; z0++;
				  FILL_PPIPED_NORM(p2,1,1,-1,1,2,3,ordinate);
				  x0 -= 2;
				  FILL_PPIPED_NORM(p3,-1,1,1,1,2,3,ordinate);
				  x0++;y0--;z0++;
				  FILL_PPIPED_NORM(p4,1,-1,1,3,2,1,ordinate);
				}
			}


		  alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f;
		  CONV_PPIPED(result,p1,alpha,beta,gamma);
		  alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f;
		  CONV_PPIPED(result,p2,alpha,beta,gamma);
		  alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax);
		  CONV_PPIPED(result,p3,alpha,beta,gamma);
		  alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f);
		  CONV_PPIPED(result,p4,alpha,beta,gamma);

		  N[ordinate] = (elem_t)result;

		} // Quintic Box Spline
		//#endif
	}

}
		*/		
		
