#ifdef WIN32
#include <windows.h>
#endif

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif

#include <stdio.h>


#include "pthread.h"
#include "glpreview.h"

static unsigned char *g_pImageContainer = NULL;
static int   g_width = 0;
static int	 g_height =0;
static int   g_interval;

static int   g_end = 0;
static int	 g_preview = 0;
static int   g_currentRasterPosX;
static int   g_currentRasterPosY;
static int   g_glutWnd = 0;

pthread_t g_previewThread;

void prevRender()
{
	if( g_end ) return ;

	if( g_pImageContainer) {
		
		glViewport(0,0,g_width,g_height);
				
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(0,g_width,0,g_height,-1,1);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
				
		glRasterPos2f(0,g_height);
		glPixelZoom(1,-1);
		glDrawPixels( g_width,g_height, 0x80E0 /*GL_BGR*/, GL_UNSIGNED_BYTE, g_pImageContainer );

		glPixelZoom(1,1);
		
		// glColor3f(0.8,0.8,0.8);
		// 		glBegin(GL_LINES);
		// 			glVertex2f(0,g_height - g_currentRasterPosY);
		// 			glVertex2f(g_currentRasterPosX, g_height - g_currentRasterPosY);			
		// 		glEnd();

		glutSwapBuffers();
		
	}
}


void prevSetSize( int w, int h)
{
	if( g_end ) return;

	// glViewport(0,0,w,h);
	// 
	// glMatrixMode(GL_PROJECTION);
	// glLoadIdentity();
	// glOrtho(0,w,0,h,-1,1);
	// glMatrixMode(GL_MODELVIEW);
	// glLoadIdentity();

	g_width = w;
	g_height = h;
}

void prevSetCurrentRasterPos( int x, int y )
{
	g_currentRasterPosX = x;
	g_currentRasterPosY = y;
}

void prevSetData( unsigned char *pData )
{
	g_pImageContainer = pData;
}

void prevIdle()
{
	if( g_end ){
		if( g_glutWnd ) {
			glutDestroyWindow( g_glutWnd );
			g_glutWnd = 0;
		}
		pthread_exit(NULL);
	}

	if( g_preview ){
		g_preview = 0;
		glutPostRedisplay();
	}
}

void prevInit()
{
	int   argc = 1;
	const char* argv = "zahid";

	glutInit( &argc, (char**)&argv );
	glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB );
	glutInitWindowSize( g_width, g_height );
	//glutInitWindowPosition( 100, 100 );
	g_glutWnd = glutCreateWindow( "Preview" );

	glViewport(0,0,g_width,g_height);
		
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0,g_width,0,g_height,-1,1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	
	glPixelStorei( GL_PACK_ALIGNMENT, 1 );
	glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );

	glutDisplayFunc( prevRender );
	glutIdleFunc( prevIdle );
	prevPreview();
	
}

void prevPreview()
{
	g_preview = 1;
}

void prevOnTimer( int n )
{		
	if( !g_end ){
		prevRender();	
		glutTimerFunc( g_interval, prevOnTimer, n );
	}
	else{
		if( g_glutWnd ) {
			glutDestroyWindow( g_glutWnd );
			g_glutWnd = 0;
		}
		pthread_exit(NULL);
	}
}

void *prevMain(void *arg)
{
	prevInit();
	
	//glutTimerFunc( g_interval, prevOnTimer, g_interval );
	g_currentRasterPosX = g_currentRasterPosY = 0;
	glutMainLoop();
	return NULL;
}

void prevBegin( unsigned int millis )
{
	g_end = 0;	
	g_interval = millis;
	int ret = pthread_create(&g_previewThread,NULL,prevMain,NULL);
	if( ret ) {
		printf("Preview thread cannot be created\n");
	}
}

void prevEnd()
{
	g_end = 1;
	void *status;
	pthread_join( g_previewThread, &status );

}
