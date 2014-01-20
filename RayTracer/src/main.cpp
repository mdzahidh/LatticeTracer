#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include <time.h>


#include <stdio.h>
//#include <GL/glut.h>

#include "ccdata.h"
#include "bccdata.h"
#include "hybriddata.h"
#include "getpot.h"
#include "potscript.h"
#include "camera.h"
#include "tracer.h"
#include "transferfunc.h"
#include "glpreview.h"
#include "imagewriter.h"
#include "function.h"

#include <mpi.h>
#define ROOT_PROC_ID 0

FILE *g_DataFile;
bool  g_bPreview = false;
CBaseData	*g_FunctionForSampling = NULL;
vector4_t    g_backgroundColor = {1,1,1,1};

int    g_numProcs;
int	   g_procID;

#ifdef WIN32
unsigned long rtGetTime()
{
	return GetTickCount();
}
#else
clock_t rtGetTime()
{
// 	struct timeval tv;
// 	gettimeofday(&tv,NULL);
// 	return tv.tv_sec*1000 + tv.tv_usec/1000;
	return clock();
}
#endif


void DataCallback( int x, int y, const vector3_t p, int n, double *pValues, vector3_t *pNormals )
{
	fprintf(g_DataFile,"%d %d; %lf %lf %lf;", x,y, p[X], p[Y], p[Z]);
	
	// Remember that 0 is the mother data;
	for(int i=0;i<n;i++){
		fprintf(g_DataFile," %lf; %lf %lf %lf;", pValues[i], pNormals[i][X], pNormals[i][Y], pNormals[i][Z]);
	}
	
	fprintf(g_DataFile,"\n");
}

void ClampColor( color_t color )
{
	if( color[0] < 0 ) color[0] = 0;
	else if( color[0] > 1 ) color[0] = 1.0f;
	
	if( color[1] < 0 ) color[1] = 0;
	else if( color[1] > 1 ) color[1] = 1.0f;
	
	if( color[2] < 0 ) color[2] = 0;
	else if( color[2] > 1 ) color[2] = 1.0f;		
}

void TakeSnapSuperSample(unsigned char *imgOut, int w, int h)
{
	unsigned char *pOffset = imgOut;
	Ray_t ray;	
	vector4_t sampleColor[5];
	vector4_t color;
	vector2_t offs[5] = { {0.0,0.0}, {0.25,0.25}, { -0.25, 0.25 }, { -0.25, -0.25 }, { 0.25, -0.25 } };

#ifdef PRINT_PROGRESS	
	double percent = 0;
	int n = 0;
	int oldN = 0;
	int progressLength = 40;
#endif	

	rtBegin();
	
	unsigned char* buffer = (unsigned char*)malloc( sizeof(unsigned char)*3*w*h);
	memset(buffer,0,sizeof(unsigned char)*h*w*3);
	
	int jtrack = 0;
	jtrack = g_procID;
	for(int j = g_procID; j< h; j += g_numProcs ){
		

#ifdef PRINT_PROGRESS
		if( g_procID == ROOT_PROC_ID){
			percent = (double) j / (double)(h - 1 - g_procID);
			n = (int)(percent * progressLength);

			if( n > oldN ){

				if( n == progressLength/4 + 1 ) printf("25%%");			
				if( n == progressLength/2 + 1 ) printf("50%%");
				if( n == (progressLength*3)/4 + 1 ) printf("75%%");

				printf(">");	
				fflush(stdout);		
				oldN = n;
			}
		}
#endif		
		pOffset = buffer + (j * (w * 3 * sizeof(unsigned char)));
		
		for(int i=0;i<w;++i){
			
			
			int hitCount = 0;
			
			color[R] = 0;
			color[G] = 0;
			color[B] = 0;
			color[A] = 0;
			
			for(int s=0;s<5;++s){
				
				int hit = 0;
				
				camGetWorldSpaceRayFromScreenSpace( i + offs[s][X], j + offs[s][Y], &ray );
				
				hit = rtTraceRay( i, j, sampleColor[s], &ray, g_backgroundColor );
				
				if( hit ){
					color[R] += sampleColor[s][R];
					color[G] += sampleColor[s][G];
					color[B] += sampleColor[s][B];
					color[A] += sampleColor[s][A];
					
					++hitCount;
				}
				  
				
			}
						
			color[R] /= hitCount;
			color[G] /= hitCount;
			color[B] /= hitCount;
			color[A] /= hitCount;						 
			
			if( hitCount )
			{
				ClampColor( color );
						
				*pOffset ++ = (int)(color[Z] * 255.0 + 0.5);
				*pOffset ++ = (int)(color[Y] * 255.0 + 0.5);
				*pOffset ++ = (int)(color[X] * 255.0 + 0.5);
			}
			else{
				*pOffset ++ = (int)(g_backgroundColor[Z] * 255 + 0.5);
				*pOffset ++ = (int)(g_backgroundColor[Y] * 255 + 0.5);
				*pOffset ++ = (int)(g_backgroundColor[X] * 255 + 0.5);
			}

			//prevPreview();
		}

		
		
		
		if(g_bPreview){
			MPI_Reduce(buffer,imgOut,sizeof(unsigned char)*3*w*h,MPI_UNSIGNED_CHAR,MPI_SUM,ROOT_PROC_ID,MPI_COMM_WORLD);
			if(g_procID == ROOT_PROC_ID){
				float ff = (jtrack * 10.0 / (float)h);
				if(  ff > 1 ){
					prevSetData(imgOut);
					prevSetCurrentRasterPos( w, j+1 );
					prevPreview();
					jtrack -= floor( ff ) * h/10.0;
				}
			}
		}
	
		jtrack += g_numProcs;
	}
		
	MPI_Reduce(buffer,imgOut,sizeof(unsigned char)*3*w*h,MPI_UNSIGNED_CHAR,MPI_SUM,ROOT_PROC_ID,MPI_COMM_WORLD);
	
	free(buffer);
	
	rtEnd();

#ifdef PRINT_PROGRESS	
	if( g_procID == ROOT_PROC_ID) 
		printf("100%%\n");
#endif	

}


void TakeSnap(unsigned char *imgOut, int w, int h)
{
	unsigned char *pOffset;
	Ray_t ray;	
	vector4_t color;

#ifdef PRINT_PROGRESS	
	double percent = 0;
	int n = 0;
	int oldN = 0;
	int progressLength = 40;
#endif	


	rtBegin();
	
	
	unsigned char* buffer = (unsigned char*)malloc( sizeof(unsigned char)*3*w*h);
	
	memset(buffer,0,sizeof(unsigned char)*h*w*3);
   	
	int start_j = (g_procID * h) / g_numProcs;
	int end_j   = ((g_procID+1) * h) / g_numProcs;
	
	int jtrack = g_procID;
	for(int j = g_procID; j< h; j += g_numProcs ){
		
		pOffset = buffer + (j * (w * 3 * sizeof(unsigned char)));		
		
#ifdef PRINT_PROGRESS
		if( g_procID == ROOT_PROC_ID ){
			percent = (double) j / (double) ( h-1-g_procID);
			n = (int)(percent * progressLength);

			if( n > oldN ){

				if( n == progressLength/4 + 1 ) printf("25%%");			
				if( n == progressLength/2 + 1 ) printf("50%%");
				if( n == (progressLength*3)/4 + 1 ) printf("75%%");

				printf(">");	
				fflush(stdout);		
				oldN = n;
			}
		}
#endif		
		for(int i=0;i<w;++i){
			
			
			
			camGetWorldSpaceRayFromScreenSpace( (double)i, (double)j, &ray );
			
	
			int hitCount = rtTraceRay( i, j, color, &ray, g_backgroundColor ); 
			
			if( hitCount )
			{
				ClampColor( color );
						
				*pOffset ++ = (int)(color[Z] * 255.0 + 0.5);
				*pOffset ++ = (int)(color[Y] * 255.0 + 0.5);
				*pOffset ++ = (int)(color[X] * 255.0 + 0.5);
			}
			else{
				*pOffset ++ = (int)(g_backgroundColor[Z] * 255 + 0.5);
				*pOffset ++ = (int)(g_backgroundColor[Y] * 255 + 0.5);
				*pOffset ++ = (int)(g_backgroundColor[X] * 255 + 0.5);
			}


			//prevPreview();
		}

		if(g_bPreview){
			MPI_Reduce(buffer,imgOut,sizeof(unsigned char)*3*w*h,MPI_UNSIGNED_CHAR,MPI_SUM,ROOT_PROC_ID,MPI_COMM_WORLD);
			if(g_procID == ROOT_PROC_ID){
				float ff = (jtrack * 10.0 / (float)h);
				if(  ff > 1 ){
					prevSetData(imgOut);
					prevSetCurrentRasterPos( w, j+1 );
					prevPreview();
					jtrack -= floor( ff ) * h/10.0;
				}
			}
		}
	
		jtrack += g_numProcs;
			
	}
		
	
	rtEnd();

#ifdef PRINT_PROGRESS	
	if( g_procID == ROOT_PROC_ID) printf("100%%\n");
#endif
	
	MPI_Reduce(buffer,imgOut,sizeof(unsigned char)*3*w*h,MPI_UNSIGNED_CHAR,MPI_SUM,ROOT_PROC_ID,MPI_COMM_WORLD);
	
	free(buffer);
}

elem_t EvaluateCallBack( const vector3_t p )
{
	if( g_FunctionForSampling ){
		return g_FunctionForSampling->Get(p);
	}

	return 0;
}

void EvaluateNormalCallBack( const vector3_t p, vector3_t N )
{
	if( g_FunctionForSampling ){
		g_FunctionForSampling->GetNormal( p, N );
		return;
	}
	else{
		N[X] = N[Y] = N[Z] = 0;
	}
}

CBaseData* CreateGrid( PotGrid *pGridScript )
{
	CBaseData *pGrid = NULL;
	CBaseData::DATA_FORMAT fmt;

	if( pGridScript->m_dataFiles.m_dataFormat == POTGRID_BYTE){
		fmt = CBaseData::DATA_FORMAT_BYTE;
	}
	else{
		fmt = CBaseData::DATA_FORMAT_FLOAT;
	}

	switch( pGridScript->m_type ){
		case POTGRID_BCC:
			
			if( pGridScript->m_dataType == POTGRID_FILE && pGridScript->m_normalType != POTGRID_FILE ){
				
				pGrid = new CBCCData ( pGridScript->m_name.c_str(),
									   pGridScript->m_dataFiles.m_scalarFile.c_str(),
									   fmt,
					                   pGridScript->m_xScale,
									   pGridScript->m_yScale,
									   pGridScript->m_zScale,
									   pGridScript->m_normalEstimatorFilter.c_str(),
									   pGridScript->m_bNormalEstimatorOnTheFly,
									   pGridScript->m_gridPointDumpFile.c_str() );
				
				pGrid->SetReconstructionFilterByName( pGridScript->m_reconstructionFilter.c_str() ); 
				pGrid->SetNormalReconstructionFilterByName( pGridScript->m_normalReconstructionFilter.c_str());

				pGrid->SetEpsilon( pGridScript->m_epsilon );

			}
			else if( pGridScript->m_dataType == POTGRID_FILE && pGridScript->m_normalType == POTGRID_FILE )
			{
				pGrid = new CBCCData( pGridScript->m_name.c_str(),
					pGridScript->m_dataFiles.m_scalarFile.c_str(),
					fmt,
					pGridScript->m_dataFiles.m_normalXFile.c_str(),
					pGridScript->m_dataFiles.m_normalYFile.c_str(),
					pGridScript->m_dataFiles.m_normalZFile.c_str(),
					pGridScript->m_xScale,
					pGridScript->m_yScale,
					pGridScript->m_zScale,					
					pGridScript->m_gridPointDumpFile.c_str() );

			   
				pGrid->SetReconstructionFilterByName( pGridScript->m_reconstructionFilter.c_str() ); 
				pGrid->SetNormalReconstructionFilterByName( pGridScript->m_normalReconstructionFilter.c_str());
				pGrid->SetEpsilon( pGridScript->m_epsilon );
			}
			else if( pGridScript->m_dataType == POTGRID_FUNCTION) {
				pGrid = CFunctionSystem::GetFunctionByName( pGridScript->m_dataFunctionName );
			}
			else if( pGridScript->m_dataType == POTGRID_SAMPLED ){

				g_FunctionForSampling = CFunctionSystem::GetFunctionByName( pGridScript->m_dataFunctionName );
				
				if( !g_FunctionForSampling ){
					return NULL;
				}
			
				pGrid = new CBCCData( pGridScript->m_name.c_str(), 
						pGridScript->m_normalEstimatorFilter.c_str(),
						pGridScript->m_bNormalEstimatorOnTheFly,
						pGridScript->m_gridPointDumpFile.c_str() );

				pGrid->SetReconstructionFilterByName( pGridScript->m_reconstructionFilter.c_str() ); 
				pGrid->SetNormalReconstructionFilterByName( pGridScript->m_normalReconstructionFilter.c_str());
				pGrid->SetEpsilon( pGridScript->m_epsilon );

				vector3_t steps;
				
				steps[X] = (pGridScript->m_dataFunctionMax[X] - pGridScript->m_dataFunctionMin[X]) / (pGridScript->m_dataFunctionGrid[X]- 0.5);
				steps[Y] = (pGridScript->m_dataFunctionMax[Y] - pGridScript->m_dataFunctionMin[Y]) / (pGridScript->m_dataFunctionGrid[Y]- 0.5);
				steps[Z] = (pGridScript->m_dataFunctionMax[Z] - pGridScript->m_dataFunctionMin[Z]) / (pGridScript->m_dataFunctionGrid[Z]- 1.0);

				printf("Scales: (XScale=%f, YScale=%f, ZScale=%f)\n",steps[X], steps[Y], steps[Z]);

				pGrid->Evaluate( pGridScript->m_dataFunctionMin, pGridScript->m_dataFunctionMax, steps, EvaluateCallBack, EvaluateNormalCallBack );

				// Resetting things to original state
				g_FunctionForSampling = NULL;

				
			}

			break;
		
		case POTGRID_CC:

			if( pGridScript->m_dataType == POTGRID_FILE && pGridScript->m_normalType != POTGRID_FILE ){
				
				pGrid = new CCCData  ( pGridScript->m_name.c_str(),
									   pGridScript->m_dataFiles.m_scalarFile.c_str(),
									   fmt,
					                   pGridScript->m_xScale,
									   pGridScript->m_yScale,
									   pGridScript->m_zScale,
									   pGridScript->m_normalEstimatorFilter.c_str(),
									   pGridScript->m_bNormalEstimatorOnTheFly,
									   pGridScript->m_gridPointDumpFile.c_str() );
				
				pGrid->SetReconstructionFilterByName( pGridScript->m_reconstructionFilter.c_str() ); 
				pGrid->SetNormalReconstructionFilterByName( pGridScript->m_normalReconstructionFilter.c_str());

				pGrid->SetEpsilon( pGridScript->m_epsilon );

			}
			else if( pGridScript->m_dataType == POTGRID_FILE && pGridScript->m_normalType == POTGRID_FILE )
			{
				pGrid = new CCCData( pGridScript->m_name.c_str(),
					pGridScript->m_dataFiles.m_scalarFile.c_str(),
					fmt,
					pGridScript->m_dataFiles.m_normalXFile.c_str(),
					pGridScript->m_dataFiles.m_normalYFile.c_str(),
					pGridScript->m_dataFiles.m_normalZFile.c_str(),
					pGridScript->m_xScale,
					pGridScript->m_yScale,
					pGridScript->m_zScale,					
					pGridScript->m_gridPointDumpFile.c_str() );

			   
				pGrid->SetReconstructionFilterByName( pGridScript->m_reconstructionFilter.c_str() ); 
				pGrid->SetNormalReconstructionFilterByName( pGridScript->m_normalReconstructionFilter.c_str());
				pGrid->SetEpsilon( pGridScript->m_epsilon );
			}
			else if( pGridScript->m_dataType == POTGRID_FUNCTION) {
				pGrid = CFunctionSystem::GetFunctionByName( pGridScript->m_dataFunctionName );
			}
			else if( pGridScript->m_dataType == POTGRID_SAMPLED ){

				g_FunctionForSampling = CFunctionSystem::GetFunctionByName( pGridScript->m_dataFunctionName );
				
				if( !g_FunctionForSampling ){
					return NULL;
				}
			
				pGrid = new CCCData( pGridScript->m_name.c_str(), 
					                  pGridScript->m_normalEstimatorFilter.c_str(), 
							  pGridScript->m_bNormalEstimatorOnTheFly,
									  pGridScript->m_gridPointDumpFile.c_str() );

				pGrid->SetReconstructionFilterByName( pGridScript->m_reconstructionFilter.c_str() ); 
				pGrid->SetNormalReconstructionFilterByName( pGridScript->m_normalReconstructionFilter.c_str());
				pGrid->SetEpsilon( pGridScript->m_epsilon );

				vector3_t steps;
				
				steps[X] = (pGridScript->m_dataFunctionMax[X] - pGridScript->m_dataFunctionMin[X]) / (pGridScript->m_dataFunctionGrid[X]-1);
				steps[Y] = (pGridScript->m_dataFunctionMax[Y] - pGridScript->m_dataFunctionMin[Y]) / (pGridScript->m_dataFunctionGrid[Y]-1);
				steps[Z] = (pGridScript->m_dataFunctionMax[Z] - pGridScript->m_dataFunctionMin[Z]) / (pGridScript->m_dataFunctionGrid[Z]-1);

				pGrid->Evaluate( pGridScript->m_dataFunctionMin, pGridScript->m_dataFunctionMax, steps, EvaluateCallBack, EvaluateNormalCallBack );

				// Resetting things to original state
				g_FunctionForSampling = NULL;

				
			}
		

			break;
		case POTGRID_HYBRID:
			printf("I shouldn't reach here !!!\n");
			break;
	}

	//Set the boundary condition for this grid
	pGrid->SetBoundaryCondition( pGridScript->m_boundaryCond.c_str() );
	pGrid->HookGradientComponentEstimators();
	return pGrid;
}

void ReportPotScript( CPotScript *pScript )
{
	printf("Experiment Type =%d\n",pScript->m_experimentType );
	printf("Seed: %ld\n",pScript->m_seed);
	printf("RandomSamples: %ld\n",pScript->m_randomSamples);

	printf("Random Bounding Box: Min: [%f %f %f] ; Max: [ %f %f %f]\n", 
		pScript->m_randomBoundingMin[X],
		pScript->m_randomBoundingMin[Y],
		pScript->m_randomBoundingMin[Z],
		pScript->m_randomBoundingMax[X],
		pScript->m_randomBoundingMax[Y],
		pScript->m_randomBoundingMax[Z]);

	printf("Number of Grids = %d\n", pScript->m_numGrids );
	printf("Main Grid Idx =  %d\n", pScript->m_motherGridIdx );

	for(int i=0;i<pScript->m_numGrids;++i){

		printf("[Grid%d]\n",i);
		printf("\tName: %s\n",pScript->m_grids[i].m_name.c_str());
		printf("\tType: %d\n",pScript->m_grids[i].m_type);
		printf("\tDataType: %d\n", pScript->m_grids[i].m_dataType);
		printf("\tNormalType: %d\n",pScript->m_grids[i].m_normalType);
		printf("\tReconstructionFilter: %s\n",pScript->m_grids[i].m_reconstructionFilter.c_str());
		printf("\tNormalReconstructionFilter: %s\n",pScript->m_grids[i].m_normalReconstructionFilter.c_str());
		printf("\tNormalEstimatorFilter: %s\n",pScript->m_grids[i].m_normalEstimatorFilter.c_str());
		printf("\tOnTheFlyNormals: %d\n",pScript->m_grids[i].m_bNormalEstimatorOnTheFly);
		printf("\tDataFunction: %s\n",pScript->m_grids[i].m_dataFunctionName.c_str() );
		printf("\tFunctionEvaluateMin: (%f %f %f)\n", pScript->m_grids[i].m_dataFunctionMin[X],pScript->m_grids[i].m_dataFunctionMin[Y],pScript->m_grids[i].m_dataFunctionMin[Z] );
		printf("\tFunctionEvaluateMax: (%f %f %f)\n", pScript->m_grids[i].m_dataFunctionMax[X],pScript->m_grids[i].m_dataFunctionMax[Y],pScript->m_grids[i].m_dataFunctionMax[Z] );
		printf("\tFunctionEvaluateGrid: (%d %d %d)\n", pScript->m_grids[i].m_dataFunctionGrid[X],pScript->m_grids[i].m_dataFunctionGrid[Y],pScript->m_grids[i].m_dataFunctionGrid[Z] );
		printf("\tEpsilon: %f\n",pScript->m_grids[i].m_epsilon );
		printf("\tXScale: %f\n",pScript->m_grids[i].m_xScale );
		printf("\tYScale: %f\n",pScript->m_grids[i].m_yScale );
		printf("\tZScale: %f\n",pScript->m_grids[i].m_zScale );
		printf("\tBoundaryCondition: %s\n",pScript->m_grids[i].m_boundaryCond.c_str());

		printf("\t\tDataFormat:  %d\n", pScript->m_grids[i].m_dataFiles.m_dataFormat );
		printf("\t\tScalarFile:  %s\n", pScript->m_grids[i].m_dataFiles.m_scalarFile.c_str() );
		printf("\t\tNormalXFile: %s\n", pScript->m_grids[i].m_dataFiles.m_normalXFile.c_str() );
		printf("\t\tNormalYFile: %s\n", pScript->m_grids[i].m_dataFiles.m_normalYFile.c_str() );
		printf("\t\tNormalZFile: %s\n", pScript->m_grids[i].m_dataFiles.m_normalZFile.c_str() );

		printf("\tGridPointDumpFile: %s\n",pScript->m_grids[i].m_gridPointDumpFile.c_str());
	}

	printf("[Output]\n");
	printf("\tTracedDumpFile: %s\n",pScript->m_output.m_traceDumpFile.c_str());
	printf("\tImage: %s\n",pScript->m_output.m_image.c_str());

	printf("[Camera]\n");
	printf("\tType: %d\n",pScript->m_camera.m_type );
	printf("\tFOV: %f\n",pScript->m_camera.m_FOV );
	printf("\tNearPlane: %f\n",pScript->m_camera.m_nearPlane);
	printf("\tFarPlane: %f\n",pScript->m_camera.m_farPlane);
	printf("\tPosition: (%f %f %f)\n",pScript->m_camera.m_position[X],pScript->m_camera.m_position[Y],pScript->m_camera.m_position[Z]);
	printf("\tUp: (%f %f %f)\n",pScript->m_camera.m_up[X],pScript->m_camera.m_up[Y],pScript->m_camera.m_up[Z]);
	printf("\tLookAt: (%f %f %f)\n",pScript->m_camera.m_lookAt[X],pScript->m_camera.m_lookAt[Y],pScript->m_camera.m_lookAt[Z]);

	printf("[Tracer]\n");
	printf("\tMode: %d\n",pScript->m_tracer.m_tracerMode );
	printf("\tIsoSurface: %f\n",pScript->m_tracer.m_isoSurface);
	printf("\tIsoSurfaceThreshold: %f\n",pScript->m_tracer.m_isoSurfaceThreshold);
	printf("\tResolution: (%d %d)\n",pScript->m_tracer.m_imgWidth, pScript->m_tracer.m_imgHeight );
	printf("\tPreview: %d\n",pScript->m_tracer.m_bPreview);
	printf("\tShadowed: %d\n",pScript->m_tracer.m_bShadowed);
	printf("\tTwoSidedLighting: %d\n",pScript->m_tracer.m_bTwoSidedLighting);
	printf("\tLighting: %d\n",pScript->m_tracer.m_bLighting);
	printf("\tSuperSampling: %d\n", pScript->m_tracer.m_bSuperSampling);
	printf("\tNumLights: %d\n",pScript->m_tracer.m_numLights);
	printf("\tSamplingRate: %f\n",pScript->m_tracer.m_samplingRate);
	printf("\tShadowOffset: %f\n",pScript->m_tracer.m_shadowOffset);
	
	printf("\tShadowAmbient: (%f %f %f %f)\n",pScript->m_tracer.m_shadowAmbient[R],
	       pScript->m_tracer.m_shadowAmbient[G],
	       pScript->m_tracer.m_shadowAmbient[B],
	       pScript->m_tracer.m_shadowAmbient[A]);
	
	printf("\tAmbient: (%f %f %f %f)\n",pScript->m_tracer.m_ambient[R],
	       pScript->m_tracer.m_ambient[G],
	       pScript->m_tracer.m_ambient[B],
	       pScript->m_tracer.m_ambient[A]);
	
	printf("\tAlphaTransfer: %s\n",pScript->m_tracer.m_alphaTransferName.c_str());
	printf("\tColorTransfer: %s\n",pScript->m_tracer.m_colorTransferName.c_str());

	for(int i=0;i<pScript->m_tracer.m_numLights;++i){
		printf("[Light%d]\n",i);
		printf("\tPosition: (%f %f %f %f)\n", pScript->m_lights[i].m_position[X],
			pScript->m_lights[i].m_position[Y],
			pScript->m_lights[i].m_position[Z],
			pScript->m_lights[i].m_position[W] );
		
		printf("\tcolor: (%f %f %f %f)\n", pScript->m_lights[i].m_color[X],
			pScript->m_lights[i].m_color[Y],
			pScript->m_lights[i].m_color[Z],
			pScript->m_lights[i].m_color[W] );
		
	}
}

void RunTracedExperiment( CPotScript *pScript )
{
		
	CBaseData **ppGrids;

	//Set up camera parameters
	camSetPerspectiveParameters( pScript->m_camera.m_FOV, 
		                        (double)pScript->m_tracer.m_imgWidth/(double)pScript->m_tracer.m_imgHeight,
								pScript->m_camera.m_nearPlane);
	

	camSetImageParameters( pScript->m_tracer.m_imgWidth, pScript->m_tracer.m_imgHeight );		
	
	camSetOrientation( pScript->m_camera.m_position, pScript->m_camera.m_lookAt,  pScript->m_camera.m_up);
	
	if( pScript->m_camera.m_type == POTCAMERA_PERSPECTIVE )
		camInitPerspectiveView();
	else
		camInitOrthogonalView();
	
	ppGrids = new CBaseData* [ pScript->m_numGrids ];
	
	for(int i=0;i<pScript->m_numGrids;++i){

		if( pScript->m_grids[i].m_type == POTGRID_HYBRID) {
			// We will do a hack here to load a Hybrid grid
			// We are making a tempory pGrid that has 
			// the same type as subtype and make the system
			// load the sub grid that way.
			printf("Loading a hybrid grid ...\n");
			PotGrid pGrid = pScript->m_grids[i];
			pGrid.m_type = pGrid.m_subType;
			
			printf("\tLoading the sub grid ... ");
			CBaseData *pSubGrid = CreateGrid( &pGrid );
			printf("DONE!");
			
			pGrid.m_dataType = POTGRID_FUNCTION;
			printf("\tLoading the function grid ... ");
			CBaseData *pFunctionGrid = CreateGrid( &pGrid );
			printf("DONE!");
			
			CBaseData *pHybridGrid = 
					new CHybridData ( pGrid.m_name.c_str(),
					                  pFunctionGrid,
							  pSubGrid );
			
			ppGrids[i] = pHybridGrid;
			printf("DONE!\n");		
		}
		else {
			ppGrids[i] = CreateGrid( &pScript->m_grids[i] );
		}		
	}	

	rtSetSamplingRate( pScript->m_tracer.m_samplingRate );
	
	CTransferFunction::SetAlphaTransferByName( pScript->m_tracer.m_alphaTransferName );
	CTransferFunction::SetColorTransferByName( pScript->m_tracer.m_colorTransferName );
	
	rtSetMode( pScript->m_tracer.m_tracerMode == POTTRACERMODE_DVR ? RT_DVR : RT_ISR );
	rtSetIsoSurfaceThreshold( pScript->m_tracer.m_isoSurfaceThreshold );
	rtSetIsoSurface( pScript->m_tracer.m_isoSurface );
	rtSet(RT_SHADOW, pScript->m_tracer.m_bShadowed );
	rtSet(RT_LIGHTING, RT_TRUE);	
	rtSet( RT_TWO_SIDED_LIGHTING, pScript->m_tracer.m_bTwoSidedLighting );
	rtSet(RT_LIGHTING, pScript->m_tracer.m_bLighting);	
	rtSetShadowAmbient( pScript->m_tracer.m_shadowAmbient);
	rtSetAmbient( pScript->m_tracer.m_ambient);
		
	rtSetShadowOffsetParam( pScript->m_tracer.m_shadowOffset );
	
	
	for(int i=0; i < pScript->m_tracer.m_numLights; ++i ){

		rtAddLight( pScript->m_lights[i].m_position, pScript->m_lights[i].m_color );
	}

	// Adding the mother grid
	rtAddData( ppGrids[ pScript->m_motherGridIdx ] );

	// Adding rest of the grids
	for(int i=0; i < pScript->m_numGrids; ++i ){

		if( i != pScript->m_motherGridIdx ){
			rtAddData( ppGrids[ i ] );
		}

	}

	g_DataFile = NULL;
	
	if( pScript->m_output.m_traceDumpFile.length() ){
		rtSetDataCallback( DataCallback );
		
	
		g_DataFile = fopen( pScript->m_output.m_traceDumpFile.c_str(),"w");
		fprintf(g_DataFile,"width=%d\n",pScript->m_tracer.m_imgWidth);
		fprintf(g_DataFile,"height=%d\n",pScript->m_tracer.m_imgHeight);
		fprintf(g_DataFile,"grids=%d\n",pScript->m_numGrids);
	
		fprintf(g_DataFile,"gridname=%s\n",ppGrids[ pScript->m_motherGridIdx]->GetName() );
	
		for(int i=0; i < pScript->m_numGrids; ++i ){
	
			if( i != pScript->m_motherGridIdx ){
				fprintf(g_DataFile,"gridname=%s\n",ppGrids[ i ]->GetName());
			}
		}
		
	}
	else{
		printf("No Trace Dump File was specified!\n");
	}

	///////

	unsigned long stTime,enTime;
	unsigned char *pImg = new unsigned char [ pScript->m_tracer.m_imgWidth * pScript->m_tracer.m_imgHeight * 3];
	imgClear(0,0,0,pScript->m_tracer.m_imgWidth,pScript->m_tracer.m_imgHeight,pImg );
	
	// TODO setup preview here
	g_bPreview = pScript->m_tracer.m_bPreview;
	
	if((g_procID == ROOT_PROC_ID) && g_bPreview){
		prevSetSize( pScript->m_tracer.m_imgWidth, pScript->m_tracer.m_imgHeight );
		prevSetData( pImg );
		prevBegin( 1000 );
	};

	 

	stTime = rtGetTime();
	
	if( pScript->m_tracer.m_bSuperSampling )
		TakeSnapSuperSample( pImg, pScript->m_tracer.m_imgWidth, pScript->m_tracer.m_imgHeight);
	else 
		TakeSnap( pImg, pScript->m_tracer.m_imgWidth, pScript->m_tracer.m_imgHeight);
	
	enTime = rtGetTime();
			
	//sprintf(fName,pScript->m_output.m_image);
	if( g_procID == ROOT_PROC_ID){
		imgSaveBMP( pScript->m_output.m_image.c_str(), pScript->m_tracer.m_imgWidth, pScript->m_tracer.m_imgHeight, pImg );
		printf("Image saved: %s\n",pScript->m_output.m_image.c_str());
	}
	
	if(g_procID == ROOT_PROC_ID){
		#ifdef WIN32
			printf("Time taken: %f seconds\n",(enTime - stTime) / 1000.0f );
		#else
			printf("Time taken: %f seconds\n",((double)(enTime - stTime)) / CLOCKS_PER_SEC );
		#endif
	}
	
	if( g_DataFile) 
		fclose(g_DataFile);

	
	for(int i=0;i<pScript->m_numGrids;++i){
		
		// Dont delete the function grids 
		if( pScript->m_grids[i].m_dataType != POTGRID_FUNCTION)  
			DELETE_SAFE( ppGrids[i] );

	}

	DELETE_SAFE( ppGrids );
	
	DELETE_SAFE_ARRAY( pImg );

	if( (g_procID == ROOT_PROC_ID) && g_bPreview ){
		prevEnd();
	}
}

void Seed( long s )
{
#ifdef WIN32
	srand(s);
#else
	srand48( s );
#endif

}

double Random()
{
#ifdef WIN32
	return (double)rand() / (double)RAND_MAX;
#else
	return drand48();
#endif
}


void RandomPosition( const vector3_t Min, const vector3_t Max, vector3_t P)
{
	P[X] = Random() * (Max[X] - Min[X]) + Min[X];
	P[Y] = Random() * (Max[Y] - Min[Y]) + Min[Y];
	P[Z] = Random() * (Max[Z] - Min[Z]) + Min[Z];	
} 

void RunRandomExperimentNormal( CPotScript *pScript )
{
	printf("Running random experiment for normal interpolation\n");
	
	if( pScript->m_seed == -1)
		Seed( time(NULL) ) ;		
	else 
		Seed( pScript->m_seed );

	CBaseData **ppGrids;

	ppGrids = new CBaseData* [ pScript->m_numGrids ];
	
	for(int i=0;i<pScript->m_numGrids;++i){

		ppGrids[i] = CreateGrid( &pScript->m_grids[i] );		
	}
	
	
	clock_t stTime = rtGetTime();
	
	for(long i=0;i<pScript->m_randomSamples;++i){
		
		vector3_t P;
		
		RandomPosition( pScript->m_randomBoundingMin, pScript->m_randomBoundingMax, P );
		
 		vector3_t N;
		ppGrids[ pScript->m_motherGridIdx ]->GetNormal( P, N );
	}
	
	clock_t endTime = rtGetTime();
	printf("Time Taken: %f\n",((double)(endTime - stTime)) / CLOCKS_PER_SEC );
	       
	DELETE_SAFE_ARRAY( ppGrids );

}

void RunRandomExperimentScalar( CPotScript *pScript )
{
	printf("Running random experiment for scalar interpolation\n");
	
	if( pScript->m_seed == -1)
		Seed( time(NULL) ) ;		
	else 
		Seed( pScript->m_seed );

	CBaseData **ppGrids;

	ppGrids = new CBaseData* [ pScript->m_numGrids ];
	
	for(int i=0;i<pScript->m_numGrids;++i){

		ppGrids[i] = CreateGrid( &pScript->m_grids[i] );		
	}
	
	
	clock_t stTime = rtGetTime();
	
	for(long i=0;i<pScript->m_randomSamples;++i){
		
		vector3_t P;
		
		RandomPosition( pScript->m_randomBoundingMin, pScript->m_randomBoundingMax, P );
		
		ppGrids[ pScript->m_motherGridIdx ]->Get( P );
	}
	
	clock_t endTime = rtGetTime();
	printf("Time Taken: %f\n",((double)(endTime - stTime)) / CLOCKS_PER_SEC );
	       
	DELETE_SAFE_ARRAY( ppGrids );

}

void RunRandomExperiment( CPotScript *pScript)
{
	if( pScript->m_seed == -1)
		Seed( time(NULL) ) ;		
	else 
		Seed( pScript->m_seed );

	CBaseData **ppGrids;

	ppGrids = new CBaseData* [ pScript->m_numGrids ];
	
	for(int i=0;i<pScript->m_numGrids;++i){

		ppGrids[i] = CreateGrid( &pScript->m_grids[i] );		
	}

	elem_t *pL2ErrorData;
	vector3_t *pL2ErrorNormal;
	elem_t *pL2ErrorNormalMag;
	elem_t *pL2ErrorAngle;

	pL2ErrorData = new elem_t [ pScript->m_numGrids - 1 ];
	pL2ErrorNormal = new vector3_t [ pScript->m_numGrids - 1];
	pL2ErrorAngle  = new elem_t [ pScript->m_numGrids - 1 ];
	pL2ErrorNormalMag = new elem_t [ pScript->m_numGrids - 1 ]; 

	memset( pL2ErrorData, 0, sizeof(elem_t) * (pScript->m_numGrids - 1) );
	memset( pL2ErrorNormal, 0, sizeof(vector3_t) * (pScript->m_numGrids - 1) );
	memset( pL2ErrorNormalMag, 0, sizeof(elem_t) * (pScript->m_numGrids - 1) );

	for(long i=0;i<pScript->m_randomSamples;++i){
		
		vector3_t P;
		
		RandomPosition( pScript->m_randomBoundingMin, pScript->m_randomBoundingMax, P );

		elem_t 		motherValue;
		vector3_t 	motherNormal;

		motherValue = ppGrids[ pScript->m_motherGridIdx ]->Get( P );
		ppGrids[ pScript->m_motherGridIdx ]->GetNormal( P, motherNormal );
 		
		int n = 0;
		for(int g=0;g<pScript->m_numGrids;++g){
			
			if( g != pScript->m_motherGridIdx ){
				elem_t value;
				vector3_t normal;

				value = ppGrids[g]->Get(  P  );
				ppGrids[g]->GetNormal( P,  normal );

				elem_t diffValue = value - motherValue;
				vector3_t diffNormal; 
					
				diffNormal[X] = normal[X] - motherNormal[X];
				diffNormal[Y] = normal[Y] - motherNormal[Y];
				diffNormal[Z] = normal[Z] - motherNormal[Z];
	
				pL2ErrorData [n] += diffValue * diffValue;
				
				pL2ErrorNormal[n][X] += diffNormal[X] * diffNormal[X];
				pL2ErrorNormal[n][Y] += diffNormal[Y] * diffNormal[Y];
				pL2ErrorNormal[n][Z] += diffNormal[Z] * diffNormal[Z];
				
				elem_t angle = (180.0 / PI) * acos( 
						glmDotProduct3d( normal, motherNormal ) / 
						(glmComputeVectorLength3d( normal) * glmComputeVectorLength3d( motherNormal ) )
					       );
				
				pL2ErrorAngle[n] +=  angle * angle;
				
				++n;		
			}
		}
	}
	
	for(int i=0;i<pScript->m_numGrids-1;i++){
		pL2ErrorNormalMag[i] = sqrt((pL2ErrorNormal[i][X] + pL2ErrorNormal[i][Y] + pL2ErrorNormal[i][Z]) / pScript->m_randomSamples );
		pL2ErrorNormal[i][X] = sqrt(pL2ErrorNormal[i][X] / pScript->m_randomSamples);
		pL2ErrorNormal[i][Y] = sqrt(pL2ErrorNormal[i][Y] / pScript->m_randomSamples);
		pL2ErrorNormal[i][Z] = sqrt(pL2ErrorNormal[i][Z] / pScript->m_randomSamples);
		pL2ErrorData[i] = sqrt(pL2ErrorData[i] / pScript->m_randomSamples);
		pL2ErrorAngle[i] = sqrt ( pL2ErrorAngle[i] / pScript->m_randomSamples );
	}
	
	FILE *fp = NULL;
	if( pScript->m_output.m_traceDumpFile.length() ){
		fp = fopen( pScript->m_output.m_traceDumpFile.c_str(), "w" );
		if( !fp ){
			printf("Cannot open output file %s \n", pScript->m_output.m_traceDumpFile.c_str() );	
		}

		if( fp ){				
			fprintf(fp,"Mother Grid: %s\n", pScript->m_grids[ pScript->m_motherGridIdx ].m_name.c_str() );		
		}
	}
	else{
		printf("Didn't specify any output text file to dump data!\n");
	}

	int n = 0;
	for(int i=0;i< pScript->m_numGrids;++i){
		if( i == pScript->m_motherGridIdx ) continue;

		printf("Grid: %s\n",pScript->m_grids[n].m_name.c_str());
		printf("\tL2 norm of data error: %lf\n",pL2ErrorData[n] );
		printf("\tL2 norm of Normal error X: %lf\n",pL2ErrorNormal[n][X] );
		printf("\tL2 norm of Normal error Y: %lf\n",pL2ErrorNormal[n][Y] );
		printf("\tL2 norm of Normal error Z: %lf\n",pL2ErrorNormal[n][Z] );
		printf("\tL2 norm of Normal Magnitude: %lf\n",pL2ErrorNormalMag[n]);
		printf("\tL2 norm of the Angle: %lf\n",pL2ErrorAngle[n]);
	
		if( fp ){
			fprintf(fp,"Grid: %s\n",pScript->m_grids[n].m_name.c_str());
			fprintf(fp,"\tL2 norm of data error: %lf\n",pL2ErrorData[n] );
			fprintf(fp,"\tL2 norm of Normal error X: %lf\n",pL2ErrorNormal[n][X] );
			fprintf(fp,"\tL2 norm of Normal error Y: %lf\n",pL2ErrorNormal[n][Y] );
			fprintf(fp,"\tL2 norm of Normal error Z: %lf\n",pL2ErrorNormal[n][Z] );
			fprintf(fp,"\tL2 norm of Normal Magnitude: %lf\n",pL2ErrorNormalMag[n]);
			fprintf(fp,"\tL2 norm of the Angle: %lf\n",pL2ErrorAngle[n]);
		}
		n++;
	}

	if( fp ) {
		printf("File written.. %s\n", pScript->m_output.m_traceDumpFile.c_str());
		fclose(fp);
	}

	for(int i=0;i<pScript->m_numGrids;++i){
		
		// Dont delete the function grids 
		if( pScript->m_grids[i].m_dataType != POTGRID_FUNCTION)  
			DELETE_SAFE( ppGrids[i] );

	}
	
	DELETE_SAFE_ARRAY( ppGrids );
	DELETE_SAFE_ARRAY( pL2ErrorData );
	DELETE_SAFE_ARRAY( pL2ErrorNormal );
	DELETE_SAFE_ARRAY( pL2ErrorNormalMag );
	DELETE_SAFE_ARRAY( pL2ErrorAngle );
	
}



void RunExperiment( std::string fname )
{
	CPotScript::InitPotScriptSystem();
	CTransferFunction::InitTransferFunctionSystem();
	CFunctionSystem::InitFunctionSystem();

	CPotScript *pScript = CPotScript::ParsePotScript( fname );
	
	if( !pScript ){
		printf("Error: failed to run pot script\n");
		return;
	}
	
	if( g_procID == ROOT_PROC_ID)
		ReportPotScript( pScript );

	if( pScript->m_experimentType == POTEXPERIMENT_TRACED){
		RunTracedExperiment( pScript );
	}
	else{
		//RunRandomExperiment( pScript );
		RunRandomExperimentScalar( pScript );
		
		//RunRandomExperimentNormal( pScript );
	}

	delete pScript;
}

void RunTestExperiment()
{
	CPotScript::InitPotScriptSystem();
	CTransferFunction::InitTransferFunctionSystem();
	CFunctionSystem::InitFunctionSystem();
	
	g_FunctionForSampling = CFunctionSystem::GetFunctionByName( "ML" );
	if( g_FunctionForSampling == NULL ){
		printf("No function found\n");
		return;
	}

	CBCCData *pBCCData = new CBCCData( "TestSampled_ML", "BCC_NORMAL_OPTIMAL26" );
	pBCCData->SetReconstructionFilterByName( "BCC_LINEAR_BOX_SPLINE" );
	pBCCData->SetNormalReconstructionFilterByName( "BCC_LINEAR_BOX_SPLINE" );

	int grid[3] = {41,41,82};
	vector3_t min   = {-1,-1,-1};
	vector3_t max   = { 1, 1, 1};
	vector3_t step = { (max[X] - min[X]) / double(grid[X] - 0.5),(max[Y] - min[Y]) / double(grid[Y] - 0.5), (max[Z] - min[Z]) / double(grid[Z] - 1) };   
	pBCCData->Evaluate( min, max, step, EvaluateCallBack, EvaluateNormalCallBack );

	int index[3] = {16,32,32};
 	
	int bccIndex[3];
	vector3_t bccCoord;
	
	vector3_t worldP;

	

	//printf("Different %f\n", actualData - approxData );	
	
	//printf("Difference Normal  %f %f %f\n", actualNormal[X] - approxNormal[X], actualNormal[Y] - approxNormal[Y], actualNormal[Z] - approxNormal[Z] );

	elem_t dataErrorSum = 0;
	for( int x = 0; x < grid[X]; ++x ){
		for( int y=0;y<grid[Y];++y){
	
			for(int z=0;z<grid[Z];++z){
	
				index[X] = x;
				index[Y] = y;
				index[Z] = z;

				pBCCData->IndexToBCC( bccIndex[X], bccIndex[Y], bccIndex[Z], index );
			
				bccCoord[X] = bccIndex[X];
				bccCoord[Y] = bccIndex[Y];
				bccCoord[Z] = bccIndex[Z];
			
			
				pBCCData->BCCToWorld( worldP, bccCoord );
			
				elem_t actualData = g_FunctionForSampling->Get( worldP );
				elem_t approxData = pBCCData->Get( worldP );
				
				vector3_t actualNormal;
				vector3_t approxNormal;
			
				g_FunctionForSampling->GetNormal( worldP, actualNormal );
				pBCCData->GetNormal( worldP, approxNormal );

				dataErrorSum += fabs(actualData - approxData);			
			}
		}
	}

	printf("Total Error: %f\n", dataErrorSum );
 
	DELETE_SAFE( pBCCData ); 
}

int main(int argc, char **argv )
{	
	GetPot cl(argc,argv);
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&g_numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD,&g_procID);
		
	if( !cl.search("-f")) {
		printf("Usage:\n");
		printf("rayRelease -f [GetPot script file name]\n");
		return -1;
	}

	
	std::string potFile = cl.next( (const char*)0 );
	if( potFile.c_str() ){
		RunExperiment( potFile );
	}	

	MPI_Finalize();
	
	return 0;
}
