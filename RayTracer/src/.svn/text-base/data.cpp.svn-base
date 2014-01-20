#include "data.h"

elem_t *** LoadVUDFileByte( const char *fname, int &xsize, int &ysize, int &zsize )
{
	int xSize,ySize,zSize;
	FILE *fp = fopen( fname, "rb" );

	if( !fp ) {
		printf("File \"%s\" cannot be opened!\n",fname );
		return NULL;
	}

	char line[256];
	char *s;
	size_t readsize;

	for(int i=0;i<6;++i){
		s =fgets(line, 256, fp ); // reads the first commented line
	}

	char dummy[20];
	//line should have dimension information now
	sscanf(line,"%s %d %d %d", dummy, &xSize, &ySize, &zSize );

	//now just read off the next 7 lines
	for(int i=0;i<7;++i){
		s = fgets(line, 256, fp ); // reads the first commented line
	}


	printf("Reading %s (%dx%dx%d).... ", fname, xSize, ySize, zSize );

	int totalElem = xSize * ySize * zSize;
	unsigned char *pMem = new unsigned char [ totalElem ];

	readsize = fread(pMem, totalElem * sizeof(unsigned char),1,fp);

	fclose(fp);

	//CreateData( xSize, ySize, zSize, xScale, yScale, zScale, pMem );


	elem_t ***pppData = new elem_t ** [ zSize + 2*GUARD_BAND];
	for( int z = 0; z < (zSize + 2*GUARD_BAND); z++ ){
		pppData[z] = new  elem_t *  [ ySize + 2*GUARD_BAND ];
		for( int y = 0; y< (ySize+2*GUARD_BAND); y++ ){
			//pppData[z][ySize - y - 1] = new elem_t [ xSize ];

			pppData[z][y] = new elem_t [ xSize + 2*GUARD_BAND ];

			//memcpy( pppData[z][ ySize - y - 1 ], pOffset, sizeof( atomic_t ) * xSize );
            memset( pppData[z][y], 0, sizeof(elem_t) * (xSize + 2*GUARD_BAND) );

		}
	}


    unsigned char *pOffset = pMem;
	for( int z = 0; z < zSize; z++ ){

		for( int y = 0; y< ySize; y++ ){

			for(int n = 0; n < xSize; ++n ){
				//pppData[z][ySize - y - 1][n] = pOffset[n];
				pppData[z+GUARD_BAND][y+GUARD_BAND][n+GUARD_BAND] = pOffset[n];
			}

			pOffset += xSize;
		}
	}


	DELETE_SAFE_ARRAY( pMem );

	printf("Done!\n");

	xsize = xSize;
	ysize = ySize;
	zsize = zSize;

	return pppData;
}

/*
template <class T> elem_t *** LoadVUDFileFloat( const char *fname, int &xsize, int &ysize, int &zsize )
{
	int xSize,ySize,zSize;
	FILE *fp = fopen( fname, "rb" );

	if( !fp ) {
		printf("File \"%s\" cannot be opened!\n",fname );
		return NULL;
	}

	char line[256];
	for(int i=0;i<6;++i){
		fgets(line, 256, fp ); // reads the first commented line
	}

	char dummy[20];
	//line should have dimension information now
	sscanf(line,"%s %d %d %d", dummy, &xSize, &ySize, &zSize );

	//now just read off the next 7 lines
	for(int i=0;i<7;++i){
		fgets(line, 256, fp ); // reads the first commented line
	}


	printf("Reading %s (%dx%dx%d).... ", fname, xSize, ySize, zSize );

	int totalElem = xSize * ySize * zSize;
	T *pMem = new T [ totalElem ];

	fread(pMem, totalElem * sizeof(T),1,fp);

	fclose(fp);

	float *pOffset = pMem;
	elem_t ***pppData = new elem_t ** [ zSize ];
	for( int z = 0; z < zSize; z++ ){
		pppData[z] = new  elem_t *  [ ySize ];
		for( int y = 0; y< ySize; y++ ){

			//pppData[z][ySize - y - 1] = new elem_t [ xSize ];

			pppData[z][y] = new elem_t [ xSize ];

			//memcpy( pppData[z][ ySize - y - 1 ], pOffset, sizeof( atomic_t ) * xSize );

			/// TODO: Optimization possible in this loop
			for(int n = 0; n < xSize; ++n ){
				pppData[z][y][n] = pOffset[n];

			}

			pOffset += xSize;
		}
	}


	DELETE_SAFE_ARRAY( pMem );

	printf("Done!\n");

	xsize = xSize;
	ysize = ySize;
	zsize = zSize;

	return pppData;
}
*/

elem_t *** LoadVUDFileFloat( const char *fname, int &xsize, int &ysize, int &zsize )
{
	int xSize,ySize,zSize;
	FILE *fp = fopen( fname, "rb" );

	if( !fp ) {
		printf("File \"%s\" cannot be opened!\n",fname );
		return NULL;
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

	int totalElem = xSize * ySize * zSize;
	float *pMem = new float [ totalElem ];

	readsize = fread(pMem, totalElem * sizeof(float),1,fp);

	fclose(fp);


	elem_t ***pppData = new elem_t ** [ zSize + 2*GUARD_BAND ];
	for( int z = 0; z < (zSize +2*GUARD_BAND); z++ ){
		pppData[z] = new  elem_t *  [ ySize + 2*GUARD_BAND ];
		for( int y = 0; y< (ySize+2*GUARD_BAND); y++ ){

			//pppData[z][ySize - y - 1] = new elem_t [ xSize ];

			pppData[z][y] = new elem_t [ xSize + 2*GUARD_BAND ];

            memset( pppData[z][y], 0, sizeof(elem_t) * (xSize + 2*GUARD_BAND) );
			//memcpy( pppData[z][ ySize - y - 1 ], pOffset, sizeof( atomic_t ) * xSize );

//			/// TODO: Optimization possible in this loop
//			for(int n = 0; n < xSize; ++n ){
//				pppData[z+GUARD_BAND][y+GUARD_BAND][n+GUARD_BAND] = pOffset[n];
//
//			}


		}
	}

    float *pOffset = pMem;
    for( int z = 0; z < zSize; z++ ){
		for( int y = 0; y< ySize; y++ ){
			for(int n = 0; n < xSize; ++n ){
				pppData[z+GUARD_BAND][y+GUARD_BAND][n+GUARD_BAND] = pOffset[n];

			}
			pOffset += xSize;
		}
	}


	DELETE_SAFE_ARRAY( pMem );

	printf("Done!\n");

	xsize = xSize;
	ysize = ySize;
	zsize = zSize;

	return pppData;
}
