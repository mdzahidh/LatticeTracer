#ifndef __IMAGE_WRITER_H__
#define __IMAGE_WRITER_H__

#ifdef __cplusplus
	extern "C" {
#endif

int imgSaveBMP( const char* fname, int w, int h, unsigned char* img ); 
void imgClear( unsigned char r, unsigned char g, unsigned char b, int w, int h, unsigned char *img);

#ifdef __cplusplus
	}
#endif

#endif
