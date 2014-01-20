// GLM - a math library for OpenGL
// Programmed by Markus Ilmola
// markus.ilmola@pp.inet.fi
// http://markus_ilmola.tripod.com
// http://personal.inet.fi/koti/markus.ilmola

// Primitiv generation stuff.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "glm.h"

void glmGeneratePlanef(float l, float w, DWORD lSeg, DWORD wSeg,
                       float *vertex,    DWORD vStride,
                       float *texCoord,  DWORD tStride,
                       float *normal,    DWORD nStride,
                       DWORD *indicates,  DWORD *vCount, DWORD *iCount)
{
  int vertices, triangles;
  BYTE *v, *t, *p;
  DWORD x, z, i, j;
  float *vert, *tex, *norm;
  if (lSeg<1) lSeg=1;
  if (wSeg<1) wSeg=1;
  if (vStride==0) vStride=3*sizeof(float);
  if (tStride==0) tStride=2*sizeof(float);
  if (nStride==0) nStride=3*sizeof(float);
  vertices=(lSeg+1)*(wSeg+1);
  triangles=2*lSeg*wSeg;
  if (vCount) (*vCount)=vertices;
  if (iCount) (*iCount)=3*triangles;

  v=(BYTE *)vertex;
  t=(BYTE *)texCoord;
  for (z=0; z<wSeg+1; z++)
  {
    for (x=0; x<lSeg+1; x++)
    {
      if (vertex)
      {
        vert=(float *)v;
        vert[0]=-(l/2)+x*(l/lSeg);
        vert[2]= 0.0f;
        vert[1]=-(w/2)+z*(w/wSeg);
        v+=vStride;
	  }
	  if (texCoord)
      {
        tex=(float *)t;
        tex[0]=x*(1.0/lSeg);
        tex[1]=z*(1.0/wSeg);
        t+=tStride;
      }
    }
  }

  // Generate normals.
  if (normal)
  {
    for (i=0, p=(BYTE *)normal; i<vertices; i++, p+=nStride)
    {
      norm=(float *)p;
	  norm[0]=0.0f;
      norm[1]=0.0f;
      norm[2]=1.0f;
    }
  }

  // Generate indicates.
  if (indicates)
  {
    for (z=0; z<wSeg; z++)
    {
      for (x=0; x<lSeg; x++)
      {
        i=6*(z*lSeg+x);
        j=z*(lSeg+1)+x;
        indicates[i+2]=j;
        indicates[i+1]=j+(lSeg+1);
        indicates[i]=j+1;
        indicates[i+5]=j+1;
        indicates[i+4]=j+(lSeg+1);
        indicates[i+3]=j+(lSeg+2);
      }
    }
  }
}

void glmGeneratePlaned(double l, double w, DWORD lSeg, DWORD wSeg,
                       double *vertex,    DWORD vStride,
                       double *texCoord,  DWORD tStride,
                       double *normal,    DWORD nStride,
                       DWORD *indicates,  DWORD *vCount, DWORD *iCount)
{
  int vertices, triangles;
  BYTE *v, *t, *p;
  DWORD x, z, i, j;
  double *vert, *tex, *norm;

  if (lSeg<1) lSeg=1;
  if (wSeg<1) wSeg=1;

  if (vStride==0) vStride=3*sizeof(double);
  if (tStride==0) tStride=2*sizeof(double);
  if (nStride==0) nStride=3*sizeof(double);

  vertices=(lSeg+1)*(wSeg+1);
  triangles=2*lSeg*wSeg;
  if (vCount) (*vCount)=vertices;
  if (iCount) (*iCount)=3*triangles;

  v=(BYTE *)vertex;
  t=(BYTE *)texCoord;
  for (z=0; z<wSeg+1; z++)
  {
    for (x=0; x<lSeg+1; x++)
    {
      if (vertex)
      {
        vert=(double *)v;
        vert[0]=-(l/2)+x*(l/lSeg);
        vert[2]= 0.0f;
        vert[1]=-(w/2)+z*(w/wSeg);
        v+=vStride;
	  }
	  if (texCoord)
      {
        tex=(double *)t;
        tex[0]=x*(1.0/lSeg);
        tex[1]=z*(1.0/wSeg);
        t+=tStride;
      }
    }
  }

  // Generate normals.
  if (normal)
  {
    for (i=0, p=(BYTE *)normal; i<vertices; i++, p+=nStride)
    {
      norm=(double *)p;
      norm[0]=0.0f;
      norm[1]=0.0f;
      norm[2]=1.0f;
    }
  }

  // Generate indicates.
  if (indicates)
  {
    for (z=0; z<wSeg; z++)
    {
      for (x=0; x<lSeg; x++)
      {
        i=6*(z*lSeg+x);
        j=z*(lSeg+1)+x;
        indicates[i+2]=j;
        indicates[i+1]=j+(lSeg+1);
        indicates[i]=j+1;
        indicates[i+5]=j+1;
        indicates[i+4]=j+(lSeg+1);
        indicates[i+3]=j+(lSeg+2);
      }
    }
  }
}

void glmGenerateSpheref(float r, DWORD stacks, DWORD slices,
                        float *vertex,    DWORD vStride,
                        float *texCoord,  DWORD tStride,
                        float *normal,    DWORD nStride,
                        DWORD *indicates,  DWORD *vCount, DWORD *iCount)
{
  int vertices, triangles;
  BYTE *v, *t, *n;
  DWORD i, j, x;
  float sine1, cosine1, sine2, cosine2;
  float *vert, *norm, *tex;

  if (stacks<2) stacks=2;
  if (slices<3) slices=3;

  if (vStride==0) vStride=3*sizeof(float);
  if (tStride==0) tStride=2*sizeof(float);
  if (nStride==0) nStride=3*sizeof(float);

  vertices=(stacks+1)*(slices+1);
  triangles=2*stacks*slices;
  if (vCount) (*vCount)=vertices;
  if (iCount) (*iCount)=3*triangles;

  v=(BYTE *)vertex;
  t=(BYTE *)texCoord;
  n=(BYTE *)normal;
  for (i=0; i<stacks+1; i++)
  {
    for (j=0; j<slices+1; j++)
    {
      sine1=sin(GLM_PI*((float)i/stacks));
      cosine1=cos(GLM_PI*((float)i/stacks));
      sine2=sin(GLM_2PI*(float)j/slices);
      cosine2=cos(GLM_2PI*(float)j/slices);

      if (vertex)
      {
        vert=(float *)v;
        vert[0]=r*sine1*sine2;
        vert[1]=r*sine1*cosine2;
        vert[2]=r*cosine1;
        v+=vStride;
      }

      if (normal)
      {
        norm=(float *)n;
        norm[0]=sine1*sine2;
        norm[1]=sine1*cosine2;
        norm[2]=cosine1;
        n+=nStride;
      }

      if (texCoord)
      {
        tex=(float *)t;
        tex[0]=j*(1.0/slices);
        tex[1]=i*(1.0/stacks);
        t+=tStride;
      }
    }
  }

  if (indicates)
  {
    for (i=0; i<stacks; i++)
    {
      for (j=0; j<slices; j++)
      {
        x=6*(i*slices+j);
        indicates[x  ]=i*(slices+1)+j;
        indicates[x+1]=i*(slices+1)+j+1;
        indicates[x+2]=(i+1)*(slices+1)+j;
        indicates[x+5]=(i+1)*(slices+1)+j;
        indicates[x+3]=i*(slices+1)+j+1;
        indicates[x+4]=(i+1)*(slices+1)+j+1;
      }
    }
  }

}

void glmGenerateSphered(double r, DWORD stacks, DWORD slices,
                        double *vertex,    DWORD vStride,
                        double *texCoord,  DWORD tStride,
                        double *normal,    DWORD nStride,
                        DWORD *indicates,  DWORD *vCount, DWORD *iCount)
{
  int vertices, triangles;
  BYTE *v, *t, *n;
  DWORD i, j, x;
  double sine1, cosine1, sine2, cosine2;
  double *vert, *norm, *tex;

  if (stacks<2) stacks=2;
  if (slices<3) slices=3;

  if (vStride==0) vStride=3*sizeof(double);
  if (tStride==0) tStride=2*sizeof(double);
  if (nStride==0) nStride=3*sizeof(double);

  vertices=(stacks+1)*(slices+1);
  triangles=2*stacks*slices;
  if (vCount) (*vCount)=vertices;
  if (iCount) (*iCount)=3*triangles;

  v=(BYTE *)vertex;
  t=(BYTE *)texCoord;
  n=(BYTE *)normal;
  for (i=0; i<stacks+1; i++)
  {
    for (j=0; j<slices+1; j++)
    {
      sine1=sin(GLM_PI*((double)i/stacks));
      cosine1=cos(GLM_PI*((double)i/stacks));
      sine2=sin(GLM_2PI*(double)j/slices);
      cosine2=cos(GLM_2PI*(double)j/slices);

      if (vertex)
      {
        vert=(double *)v;
        vert[0]=r*sine1*sine2;
        vert[1]=r*sine1*cosine2;
        vert[2]=r*cosine1;
        v+=vStride;
      }

      if (normal)
      {
        norm=(double *)n;
        norm[0]=sine1*sine2;
        norm[1]=sine1*cosine2;
        norm[2]=cosine1;
        n+=nStride;
      }

      if (texCoord)
      {
        tex=(double *)t;
        tex[0]=j*(1.0/slices);
        tex[1]=i*(1.0/stacks);
        t+=tStride;
      }
    }
  }

  if (indicates)
  {
    for (i=0; i<stacks; i++)
    {
      for (j=0; j<slices; j++)
      {
        x=6*(i*slices+j);
        indicates[x  ]=i*(slices+1)+j;
        indicates[x+1]=i*(slices+1)+j+1;
        indicates[x+2]=(i+1)*(slices+1)+j;
        indicates[x+5]=(i+1)*(slices+1)+j;
		indicates[x+3]=i*(slices+1)+j+1;
        indicates[x+4]=(i+1)*(slices+1)+j+1;
      }
    }
  }
}


void glmGenerateCylinderf(float r, float h, DWORD stacks, DWORD slices,
                          float *vertex,    DWORD vStride,
                          float *texCoord,  DWORD tStride,
                          float *normal,    DWORD nStride,
                          DWORD *indicates,  DWORD *vCount, DWORD *iCount)
{
  int vertices, triangles;
  BYTE *v, *t, *n;
  DWORD i, x, y;
  float sine, cosine;
  float *vert, *norm, *tex;

  if (stacks<1) stacks=1;
  if (slices<3) slices=3;

  if (vStride==0) vStride=3*sizeof(float);
  if (tStride==0) tStride=2*sizeof(float);
  if (nStride==0) nStride=3*sizeof(float);

  vertices=(stacks+1)*(slices+1);
  triangles=2*stacks*slices;
  if (vCount) (*vCount)=vertices;
  if (iCount) (*iCount)=3*triangles;

  v=(BYTE *)vertex;
  t=(BYTE *)texCoord;
  n=(BYTE *)normal;
  for (y=0; y<stacks+1; y++)
  {
    for (x=0; x<slices+1; x++)
    {
      sine=sin(GLM_2PI*(float)x/slices);
      cosine=cos(GLM_2PI*(float)x/slices);

      if (vertex)
      {
        vert=(float *)v;
        vert[1]=r*sine;
        vert[0]=r*cosine;
        vert[2]=-(h/2)+y*(h/stacks);
        v+=vStride;
	  }

      if (normal)
      {
        norm=(float *)n;
        norm[1]=sine;
        norm[0]=cosine;
        norm[2]=0;
        n+=nStride;
      }

      if (texCoord)
      {
        tex=(float *)t;
        tex[0]=x*(1.0/slices);
        tex[1]=y*(1.0/stacks);
        t+=tStride;
	  }
    }
  }

  if (indicates)
  {
    for (y=0; y<stacks; y++)
    {
      for (x=0; x<slices; x++)
      {
        i=6*(y*slices+x);
        indicates[i  ]=y*(slices+1)+x;
        indicates[i+1]=y*(slices+1)+x+1;
        indicates[i+2]=(y+1)*(slices+1)+x;
        indicates[i+5]=(y+1)*(slices+1)+x;
		indicates[i+3]=y*(slices+1)+x+1;
        indicates[i+4]=(y+1)*(slices+1)+x+1;
      }
    }
  }
}

void glmGenerateCylinderd(double r, double h, DWORD stacks, DWORD slices,
                          double *vertex,    DWORD vStride,
                          double *texCoord,  DWORD tStride,
                          double *normal,    DWORD nStride,
                          DWORD *indicates,  DWORD *vCount, DWORD *iCount)
{
  int vertices, triangles;
  BYTE *v, *t, *n;
  DWORD i, x, y;
  double sine, cosine;
  double *vert, *norm, *tex;

  if (stacks<1) stacks=1;
  if (slices<3) slices=3;

  if (vStride==0) vStride=3*sizeof(double);
  if (tStride==0) tStride=2*sizeof(double);
  if (nStride==0) nStride=3*sizeof(double);

  vertices=(stacks+1)*(slices+1);
  triangles=2*stacks*slices;
  if (vCount) (*vCount)=vertices;
  if (iCount) (*iCount)=3*triangles;

  v=(BYTE *)vertex;
  t=(BYTE *)texCoord;
  n=(BYTE *)normal;
  for (y=0; y<stacks+1; y++)
  {
    for (x=0; x<slices+1; x++)
    {
      sine=sin(GLM_2PI*(double)x/slices);
      cosine=cos(GLM_2PI*(double)x/slices);

      if (vertex)
      {
        vert=(double *)v;
        vert[1]=r*sine;
        vert[0]=r*cosine;
        vert[2]=-(h/2)+y*(h/stacks);
        v+=vStride;
	  }

      if (normal)
      {
        norm=(double *)n;
        norm[1]=sine;
        norm[0]=cosine;
        norm[2]=0;
        n+=nStride;
      }

      if (texCoord)
      {
        tex=(double *)t;
        tex[0]=x*(1.0/slices);
        tex[1]=y*(1.0/stacks);
        t+=tStride;
	  }
    }
  }

  if (indicates)
  {
    for (y=0; y<stacks; y++)
    {
      for (x=0; x<slices; x++)
      {
        i=6*(y*slices+x);
        indicates[i  ]=y*(slices+1)+x;
        indicates[i+1]=y*(slices+1)+x+1;
        indicates[i+2]=(y+1)*(slices+1)+x;
        indicates[i+5]=(y+1)*(slices+1)+x;
		indicates[i+3]=y*(slices+1)+x+1;
        indicates[i+4]=(y+1)*(slices+1)+x+1;
      }
    }
  }
}


void glmGenerateDiscf(float r, DWORD slices,
                      float *vertex,    DWORD vStride,
                      float *texCoord,  DWORD tStride,
                      float *normal,    DWORD nStride,
                      DWORD *indicates,  DWORD *vCount, DWORD *iCount)
{
  int vertices, triangles;
  BYTE *v, *t, *p;
  float sine, cosine;
  float *vert, *tex, *norm;
  DWORD i, j;

  if (slices<3) slices=3;

  if (vStride==0) vStride=3*sizeof(float);
  if (tStride==0) tStride=2*sizeof(float);
  if (nStride==0) nStride=3*sizeof(float);

  vertices=slices+1;
  triangles=slices;
  if (vCount) (*vCount)=vertices;
  if (iCount) (*iCount)=3*triangles;

  v=(BYTE *)vertex;
  t=(BYTE *)texCoord;
  for (i=0; i<slices; i++)
  {
    sine=sin(GLM_2PI*(float)i/slices);
    cosine=cos(GLM_2PI*(float)i/slices);

    if (vertex)
    {
      vert=(float *)v;
      vert[0]=r*sine;
      vert[1]=r*cosine;
      vert[2]=0;
      v+=vStride;
    }

    if (texCoord)
    {
      tex=(float *)t;
      tex[0]=0.5*sine+0.5;
      tex[1]=0.5*cosine+0.5;
      t+=tStride;
    }
  }

  if (vertex)
  {
    vert=(float *)v;
    vert[0]=0;
    vert[1]=0;
    vert[2]=0;
  }

  if (texCoord)
  {
    tex=(float *)t;
    tex[0]=0.5;
    tex[1]=0.5;
  }

  // Generate normals.
  if (normal)
  {
    for (i=0, p=(BYTE *)normal; i<vertices+1; i++, p+=nStride)
    {
      norm=(float *)p;
      norm[0]=0.0f;
      norm[1]=0.0f;
      norm[2]=1.0f;
    }
  }


  if (indicates)
  {
    for (i=0; i<slices; i++)
    {
      j=3*i;
      indicates[j]=i+1;
      if (i==slices-1) indicates[j]=0;
	  indicates[j+1]=i;
      indicates[j+2]=vertices-1;
    }
  }

}

void glmGenerateDiscd(double r, DWORD slices,
                      double *vertex,    DWORD vStride,
                      double *texCoord,  DWORD tStride,
                      double *normal,    DWORD nStride,
                      DWORD *indicates,  DWORD *vCount, DWORD *iCount)
{
  int vertices, triangles;
  BYTE *v, *t, *p;
  double sine, cosine;
  double *vert, *tex, *norm;
  DWORD i, j;

  if (slices<3) slices=3;

  if (vStride==0) vStride=3*sizeof(double);
  if (tStride==0) tStride=2*sizeof(double);
  if (nStride==0) nStride=3*sizeof(double);

  vertices=slices+1;
  triangles=slices;
  if (vCount) (*vCount)=vertices;
  if (iCount) (*iCount)=3*triangles;

  v=(BYTE *)vertex;
  t=(BYTE *)texCoord;
  for (i=0; i<slices; i++)
  {
    sine=sin(GLM_2PI*(double)i/slices);
    cosine=cos(GLM_2PI*(double)i/slices);

    if (vertex)
    {
      vert=(double *)v;
      vert[0]=r*sine;
      vert[1]=r*cosine;
      vert[2]=0;
      v+=vStride;
    }

    if (texCoord)
    {
      tex=(double *)t;
      tex[0]=0.5*sine+0.5;
      tex[1]=0.5*cosine+0.5;
      t+=tStride;
    }
  }

  if (vertex)
  {
    vert=(double *)v;
    vert[0]=0;
    vert[1]=0;
    vert[2]=0;
  }

  if (texCoord)
  {
    tex=(double *)t;
    tex[0]=0.5;
    tex[1]=0.5;
  }

  // Generate normals.
  if (normal)
  {
    for (i=0, p=(BYTE *)normal; i<vertices+1; i++, p+=nStride)
    {
      norm=(double *)p;
      norm[0]=0.0f;
      norm[1]=0.0f;
      norm[2]=1.0f;
    }
  }


  if (indicates)
  {
    for (i=0; i<slices; i++)
    {
      j=3*i;
      indicates[j]=i+1;
      if (i==slices-1) indicates[j]=0;
	  indicates[j+1]=i;
      indicates[j+2]=vertices-1;
    }
  }
}


void glmGenerateTorusf(float a, float c, DWORD slices, DWORD stacks,
                       float *vertex,    DWORD vStride,
                       float *texCoord,  DWORD tStride,
                       float *normal,    DWORD nStride,
                       DWORD *indicates,  DWORD *vCount, DWORD *iCount)
{
  int vertices, triangles;
  BYTE *v, *t, *n;
  DWORD i, j, x, y;
  float sine1, cosine1, sine2, cosine2;
  float *vert, *tex, *norm;

  if (slices<3) slices=3;
  if (stacks<3) stacks=3;

  if (vStride==0) vStride=3*sizeof(float);
  if (tStride==0) tStride=2*sizeof(float);
  if (nStride==0) nStride=3*sizeof(float);

  vertices=(slices+1)*(stacks+1);
  triangles=2*slices*stacks;
  if (vCount) (*vCount)=vertices;
  if (iCount) (*iCount)=3*triangles;

  v=(BYTE *)vertex;
  t=(BYTE *)texCoord;
  n=(BYTE *)normal;
  for (i=0; i<stacks+1; i++)
  {
    for (j=0; j<slices+1; j++)
    {
      sine1=sin(GLM_2PI*(float)i/stacks);
      cosine1=cos(GLM_2PI*(float)i/stacks);
      sine2=sin(GLM_2PI*(float)j/slices);
      cosine2=cos(GLM_2PI*(float)j/slices);

	  if (vertex)
      {
        vert=(float *)v;
        vert[0]=(a+c*sine2)*cosine1;
        vert[1]=(a+c*sine2)*sine1;
        vert[2]=c*cosine2;
        v+=vStride;
      }

	  if (normal)
      {
        norm=(float *)n;
        norm[0]=sine2*cosine1;
        norm[1]=sine2*sine1;
        norm[2]=cosine2;
        n+=nStride;
      }

      if (texCoord)
      {
        tex=(float *)t;
        tex[0]=j*(1.0/slices);
        tex[1]=i*(1.0/stacks);
        t+=tStride;
      }
    }
  }


  if (indicates)
  {
    for (y=0; y<stacks; y++)
    {
      for (x=0; x<slices; x++)
      {
        i=6*(y*slices+x);
        indicates[i  ]=y*(slices+1)+x;
        indicates[i+1]=y*(slices+1)+x+1;
        indicates[i+2]=(y+1)*(slices+1)+x;
        indicates[i+5]=(y+1)*(slices+1)+x;
		indicates[i+3]=y*(slices+1)+x+1;
        indicates[i+4]=(y+1)*(slices+1)+x+1;
      }
    }
  }
}

void glmGenerateTorusd(double a, double c, DWORD slices, DWORD stacks,
                       double *vertex,    DWORD vStride,
                       double *texCoord,  DWORD tStride,
                       double *normal,    DWORD nStride,
                       DWORD *indicates,  DWORD *vCount, DWORD *iCount)
{
  int vertices, triangles;
  BYTE *v, *t, *n;
  DWORD i, j, x, y;
  double sine1, cosine1, sine2, cosine2;
  double *vert, *tex, *norm;

  if (slices<3) slices=3;
  if (stacks<3) stacks=3;

  if (vStride==0) vStride=3*sizeof(double);
  if (tStride==0) tStride=2*sizeof(double);
  if (nStride==0) nStride=3*sizeof(double);

  vertices=(slices+1)*(stacks+1);
  triangles=2*slices*stacks;
  if (vCount) (*vCount)=vertices;
  if (iCount) (*iCount)=3*triangles;

  v=(BYTE *)vertex;
  t=(BYTE *)texCoord;
  n=(BYTE *)normal;
  for (i=0; i<stacks+1; i++)
  {
    for (j=0; j<slices+1; j++)
    {
      sine1=sin(GLM_2PI*(double)i/stacks);
      cosine1=cos(GLM_2PI*(double)i/stacks);
      sine2=sin(GLM_2PI*(double)j/slices);
      cosine2=cos(GLM_2PI*(double)j/slices);


	  if (vertex)
      {
        vert=(double *)v;
        vert[0]=(a+c*sine2)*cosine1;
        vert[1]=(a+c*sine2)*sine1;
        vert[2]=c*cosine2;
        v+=vStride;
      }

	  if (normal)
      {
        norm=(double *)n;
        norm[0]=sine2*cosine1;
        norm[1]=sine2*sine1;
        norm[2]=cosine2;
        n+=nStride;
      }

      if (texCoord)
      {
        tex=(double *)t;
        tex[0]=j*(1.0/slices);
        tex[1]=i*(1.0/stacks);
        t+=tStride;
      }
    }
  }


  if (indicates)
  {
    for (y=0; y<stacks; y++)
    {
      for (x=0; x<slices; x++)
      {
        i=6*(y*slices+x);
        indicates[i  ]=y*(slices+1)+x;
        indicates[i+1]=y*(slices+1)+x+1;
        indicates[i+2]=(y+1)*(slices+1)+x;
        indicates[i+5]=(y+1)*(slices+1)+x;
		indicates[i+3]=y*(slices+1)+x+1;
        indicates[i+4]=(y+1)*(slices+1)+x+1;
      }
    }
  }
}



void glmGenerateConef(float r, float h, DWORD stacks, DWORD slices,
                      float *vertex,    DWORD vStride,
                      float *texCoord,  DWORD tStride,
                      float *normal,    DWORD nStride,
                      DWORD *indicates,  DWORD *vCount, DWORD *iCount)
{
  int vertices, triangles;
  BYTE *v, *t, *n;
  float nk;
  DWORD x, y, i;
  float sine, cosine;
  float *vert, *tex, *norm;

  if (stacks<1) stacks=1;
  if (slices<3) slices=3;

  if (vStride==0) vStride=3*sizeof(float);
  if (tStride==0) tStride=2*sizeof(float);
  if (nStride==0) nStride=3*sizeof(float);

  vertices=(stacks+1)*(slices+1);
  triangles=2*stacks*slices;
  if (vCount) (*vCount)=vertices;
  if (iCount) (*iCount)=3*triangles;

  v=(BYTE *)vertex;
  t=(BYTE *)texCoord;
  n=(BYTE *)normal;
  nk=h/sqrt(h*h+r*r);
  for (y=0; y<stacks+1; y++)
  {
    for (x=0; x<slices+1; x++)
    {
      sine=sin(GLM_2PI*(float)x/slices);
      cosine=cos(GLM_2PI*(float)x/slices);

      if (vertex)
      {
        vert=(float *)v;
        vert[1]=(1.0-(float)y/stacks)*(r*sine);
        vert[0]=(1.0-(float)y/stacks)*(r*cosine);
        vert[2]=-(h/2)+y*(h/stacks);
        v+=vStride;
	  }

      if (normal)
      {
        norm=(float *)n;
        norm[1]=sine*nk;
        norm[0]=cosine*nk;
        norm[2]=nk;
        n+=nStride;
      }

      if (texCoord)
      {
        tex=(float *)t;
        tex[0]=x*(1.0/slices);
        tex[1]=y*(1.0/stacks);
        t+=tStride;
	  }
    }
  }

  if (indicates)
  {
    for (y=0; y<stacks; y++)
    {
      for (x=0; x<slices; x++)
      {
        i=6*(y*slices+x);
        indicates[i  ]=y*(slices+1)+x;
        indicates[i+1]=y*(slices+1)+x+1;
        indicates[i+2]=(y+1)*(slices+1)+x;
        indicates[i+5]=(y+1)*(slices+1)+x;
		indicates[i+3]=y*(slices+1)+x+1;
        indicates[i+4]=(y+1)*(slices+1)+x+1;
      }
    }
  }
}


void glmGenerateConed(double r, double h, DWORD stacks, DWORD slices,
                      double *vertex,    DWORD vStride,
                      double *texCoord,  DWORD tStride,
                      double *normal,    DWORD nStride,
                      DWORD *indicates,  DWORD *vCount, DWORD *iCount)
{
  int vertices, triangles;
  BYTE *v, *t, *n;
  double nk;
  DWORD x, y, i;
  double sine, cosine;
  double *vert, *tex, *norm;

  if (stacks<1) stacks=1;
  if (slices<3) slices=3;

  if (vStride==0) vStride=3*sizeof(double);
  if (tStride==0) tStride=2*sizeof(double);
  if (nStride==0) nStride=3*sizeof(double);

  vertices=(stacks+1)*(slices+1);
  triangles=2*stacks*slices;
  if (vCount) (*vCount)=vertices;
  if (iCount) (*iCount)=3*triangles;

  v=(BYTE *)vertex;
  t=(BYTE *)texCoord;
  n=(BYTE *)normal;
  nk=h/sqrt(h*h+r*r);
  for (y=0; y<stacks+1; y++)
  {
    for (x=0; x<slices+1; x++)
    {
      sine=sin(GLM_2PI*(double)x/slices);
      cosine=cos(GLM_2PI*(double)x/slices);

      if (vertex)
      {
        vert=(double *)v;
        vert[1]=(1.0-(double)y/stacks)*(r*sine);
        vert[0]=(1.0-(double)y/stacks)*(r*cosine);
        vert[2]=-(h/2)+y*(h/stacks);
        v+=vStride;
	  }

      if (normal)
      {
        norm=(double *)n;
        norm[1]=sine*nk;
        norm[0]=cosine*nk;
        norm[2]=nk;
        n+=nStride;
      }

      if (texCoord)
      {
        tex=(double *)t;
        tex[0]=x*(1.0/slices);
        tex[1]=y*(1.0/stacks);
        t+=tStride;
	  }
    }
  }

  if (indicates)
  {
    for (y=0; y<stacks; y++)
    {
      for (x=0; x<slices; x++)
      {
        i=6*(y*slices+x);
        indicates[i  ]=y*(slices+1)+x;
        indicates[i+1]=y*(slices+1)+x+1;
        indicates[i+2]=(y+1)*(slices+1)+x;
        indicates[i+5]=(y+1)*(slices+1)+x;
		indicates[i+3]=y*(slices+1)+x+1;
        indicates[i+4]=(y+1)*(slices+1)+x+1;
      }
    }
  }
}

