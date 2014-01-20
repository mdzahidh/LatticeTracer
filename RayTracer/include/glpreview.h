#ifndef __GL_PREVIEW__
#define __GL_PREVIEW__


typedef void (*fnProcess)(void);



void prevSetSize( int w, int h );
void prevBegin( unsigned int millis );
void prevEnd();
void prevSetData( unsigned char *pData );
void prevPreview();
void prevSetCurrentRasterPos( int x, int y );

#endif
