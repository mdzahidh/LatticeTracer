#include "potscript.h"
#include "getpot.h"

std::map <std::string,POTGRIDTYPE>       CPotScript::m_potGridTypeMap;
std::map <std::string,POTGRIDDATATYPE>   CPotScript::m_potGridDataTypeMap;
std::map <std::string,POTCAMERATYPE>     CPotScript::m_potCameraTypeMap;
std::map <std::string,POTEXPERIMENTTYPE> CPotScript::m_potExperimentTypeMap;
std::map <std::string, bool>             CPotScript::m_potOnOffMap;  
std::map <std::string,POTGRIDDATAFORMAT> CPotScript::m_potGridDataFormatMap;
std::map <std::string,POTTRACERMODE>     CPotScript::m_potTracerModeMap;

void CPotScript::InitPotScriptSystem()
{
	m_potGridTypeMap["CC"]				= POTGRID_CC;
	m_potGridTypeMap["BCC"]				= POTGRID_BCC;
	m_potGridTypeMap["HYBRID"]			= POTGRID_HYBRID;

	m_potGridDataTypeMap["FILE"]		= POTGRID_FILE;
	m_potGridDataTypeMap["FUNCTION"]	= POTGRID_FUNCTION;
	m_potGridDataTypeMap["SAMPLED"]		= POTGRID_SAMPLED;

	m_potCameraTypeMap["PERSPECTIVE"]	= POTCAMERA_PERSPECTIVE;
	m_potCameraTypeMap["ORTHO"]			= POTCAMERA_ORTHOGRAPHIC;

	m_potExperimentTypeMap["RANDOM"] = POTEXPERIMENT_RANDOM;
	m_potExperimentTypeMap["TRACED"] = POTEXPERIMENT_TRACED;

	m_potOnOffMap["OFF"] = false;
	m_potOnOffMap["ON"] = true;

	m_potGridDataFormatMap["BYTE"]  = POTGRID_BYTE;
	m_potGridDataFormatMap["FLOAT"]  = POTGRID_FLOAT;

	m_potTracerModeMap["DVR"] = POTTRACERMODE_DVR;
	m_potTracerModeMap["ISR"] = POTTRACERMODE_ISR;
}

CPotScript* CPotScript::ParsePotScript( std::string fname )
{

	CPotScript *pScript = new CPotScript();

	GetPot potScript( fname.c_str() );
	std::string token;


	token = potScript("ExperimentType","TRACED");
	
	pScript->m_experimentType = m_potExperimentTypeMap[token];
	pScript->m_seed = potScript("Seed",-1);
	pScript->m_randomSamples = potScript("RandomSamples",0);
	
	pScript->m_randomBoundingMin[X] = potScript("RandomBoundingBoxMin",0.0,X);
	pScript->m_randomBoundingMin[Y] = potScript("RandomBoundingBoxMin",0.0,Y);
	pScript->m_randomBoundingMin[Z] = potScript("RandomBoundingBoxMin",0.0,Z);
	
	pScript->m_randomBoundingMax[X] = potScript("RandomBoundingBoxMax",0.0,X);
	pScript->m_randomBoundingMax[Y] = potScript("RandomBoundingBoxMax",0.0,Y);
	pScript->m_randomBoundingMax[Z] = potScript("RandomBoundingBoxMax",0.0,Z);
	
	pScript->m_numGrids = potScript("NumberOfGrids",0);
	pScript->m_motherGridIdx =  potScript("MainGridIdx",0);
	
	if( pScript->m_numGrids > 8 ){
		printf("Error: System can only support maximum of %d grids\n",MAX_GRIDS);
		DELETE_SAFE( pScript );
		return NULL;
	};

	for(int i=0;i<pScript->m_numGrids;i++){
		std::string secName;
		char buffer[20];
		
		sprintf(buffer,"%d/",i);
		secName = std::string("Grid") + std::string(buffer);
		
		potScript.set_prefix(secName.c_str());

		pScript->m_grids[i].m_name              = potScript("Name","");
		pScript->m_grids[i].m_type              = m_potGridTypeMap[potScript("Type","")];
		pScript->m_grids[i].m_subType           = m_potGridTypeMap[potScript("SubType","")];
		pScript->m_grids[i].m_dataType          = m_potGridDataTypeMap[potScript("DataType","")];
		pScript->m_grids[i].m_normalType        = m_potGridDataTypeMap[potScript("NormalType","")];
		
		//TODO
		pScript->m_grids[i].m_reconstructionFilter    = potScript("ReconstructionFilter","");
		pScript->m_grids[i].m_normalReconstructionFilter = potScript("NormalReconstructionFilter","");
		pScript->m_grids[i].m_normalEstimatorFilter   = potScript("NormalEstimatorFilter","");
		pScript->m_grids[i].m_bNormalEstimatorOnTheFly = m_potOnOffMap[ potScript("OnTheFlyNormals","OFF") ];
				
		//if the normals are read from the file then dont
		//do any on the fly computation of normals
		if( pScript->m_grids[i].m_normalType == POTGRID_FILE ){
			pScript->m_grids[i].m_bNormalEstimatorOnTheFly = false;
		}
		
		
		pScript->m_grids[i].m_dataFunctionName  = potScript("DataFunction","");
		
		pScript->m_grids[i].m_dataFunctionMin[X]   = potScript("FunctionEvaluateMin",0.0,X);
		pScript->m_grids[i].m_dataFunctionMin[Y]   = potScript("FunctionEvaluateMin",0.0,Y);
		pScript->m_grids[i].m_dataFunctionMin[Z]   = potScript("FunctionEvaluateMin",0.0,Z);

		pScript->m_grids[i].m_dataFunctionMax[X]   = potScript("FunctionEvaluateMax",0.0,X);
		pScript->m_grids[i].m_dataFunctionMax[Y]   = potScript("FunctionEvaluateMax",0.0,Y);
		pScript->m_grids[i].m_dataFunctionMax[Z]   = potScript("FunctionEvaluateMax",0.0,Z);

		pScript->m_grids[i].m_dataFunctionGrid[X] =  potScript("FunctionEvaluateGrid",0,X);
		pScript->m_grids[i].m_dataFunctionGrid[Y] =  potScript("FunctionEvaluateGrid",0,Y);
		pScript->m_grids[i].m_dataFunctionGrid[Z] =  potScript("FunctionEvaluateGrid",0,Z);

		pScript->m_grids[i].m_epsilon           = potScript("Epsilon",0.0);
		pScript->m_grids[i].m_xScale            = potScript("XScale", 0.0);
		pScript->m_grids[i].m_yScale            = potScript("YScale", 0.0);
		pScript->m_grids[i].m_zScale            = potScript("ZScale", 0.0);
		pScript->m_grids[i].m_boundaryCond	= potScript("BoundaryCondition","");
		pScript->m_grids[i].m_gridPointDumpFile = potScript("GridPointDumpFile", "");


		secName = secName + std::string("Data/");

		potScript.set_prefix(secName.c_str());

		pScript->m_grids[i].m_dataFiles.m_scalarFile  = potScript("ScalarFile", "");
		pScript->m_grids[i].m_dataFiles.m_normalXFile = potScript("NormalXFile", "");
		pScript->m_grids[i].m_dataFiles.m_normalYFile = potScript("NormalYFile", "");
		pScript->m_grids[i].m_dataFiles.m_normalZFile = potScript("NormalZFile", "");
		pScript->m_grids[i].m_dataFiles.m_dataFormat  = m_potGridDataFormatMap[ potScript("DataFormat","BYTE") ];

	}

	potScript.set_prefix("Output/");
	pScript->m_output.m_traceDumpFile = potScript("TracedDumpFile","");
	pScript->m_output.m_image = potScript("Image","output.bmp");

	potScript.set_prefix("Camera/");
	pScript->m_camera.m_type = m_potCameraTypeMap[ potScript("Type","") ];
	pScript->m_camera.m_FOV  = potScript("FOV",0.0);
	pScript->m_camera.m_nearPlane = potScript("NearPlane",0.0 );
	pScript->m_camera.m_farPlane = potScript("FarPlane",0.0 );
	
	pScript->m_camera.m_position[X] = potScript("Position",0.0,X );
	pScript->m_camera.m_position[Y] = potScript("Position",0.0,Y );
	pScript->m_camera.m_position[Z] = potScript("Position",0.0,Z );

	pScript->m_camera.m_up[X] = potScript("Up",0.0,X );
	pScript->m_camera.m_up[Y] = potScript("Up",0.0,Y );
	pScript->m_camera.m_up[Z] = potScript("Up",0.0,Z );
	
	pScript->m_camera.m_lookAt[X] = potScript("LookAt",0.0,X );
	pScript->m_camera.m_lookAt[Y] = potScript("LookAt",0.0,Y );
	pScript->m_camera.m_lookAt[Z] = potScript("LookAt",0.0,Z );
	
	potScript.set_prefix("Tracer/");
	pScript->m_tracer.m_tracerMode	        = m_potTracerModeMap[ potScript("Mode","") ];
	pScript->m_tracer.m_isoSurface	        = potScript("IsoSurface", 0.0);
	pScript->m_tracer.m_isoSurfaceThreshold = potScript("IsoSurfaceThreshold",EPSILON); 
	pScript->m_tracer.m_imgWidth            = potScript("Resolution",0,0 );
	pScript->m_tracer.m_imgHeight           = potScript("Resolution",0,1 );
	pScript->m_tracer.m_bPreview            = m_potOnOffMap[ potScript("Preview","OFF") ];
	pScript->m_tracer.m_bShadowed           = m_potOnOffMap[ potScript("Shadowed","ON") ];
	pScript->m_tracer.m_bTwoSidedLighting   = m_potOnOffMap[ potScript("TwoSidedLighting","ON") ];
	pScript->m_tracer.m_bLighting           = m_potOnOffMap[ potScript("Lighting","ON") ];
	printf("Two-Sided %d\n",pScript->m_tracer.m_bTwoSidedLighting);
	pScript->m_tracer.m_bSuperSampling    = m_potOnOffMap[ potScript("SuperSampling","OFF") ];
	
	pScript->m_tracer.m_shadowAmbient[X] = potScript("ShadowAmbient",0.3,X);
	pScript->m_tracer.m_shadowAmbient[Y] = potScript("ShadowAmbient",0.3,Y);
	pScript->m_tracer.m_shadowAmbient[Z] = potScript("ShadowAmbient",0.3,Z);
	pScript->m_tracer.m_shadowAmbient[W] = potScript("ShadowAmbient",0.0,W);
	
	pScript->m_tracer.m_ambient[X] = potScript("Ambient",0.0,X);
	pScript->m_tracer.m_ambient[Y] = potScript("Ambient",0.0,Y);
	pScript->m_tracer.m_ambient[Z] = potScript("Ambient",0.0,Z);
	pScript->m_tracer.m_ambient[W] = potScript("Ambient",0.0,W);
	
	pScript->m_tracer.m_numLights         = potScript("NumLights",0);
	pScript->m_tracer.m_samplingRate      = potScript("SamplingRate",0.0);
	pScript->m_tracer.m_shadowOffset      = potScript("ShadowOffset",0.0);
	pScript->m_tracer.m_alphaTransferName = potScript("AlphaTransfer","");
	pScript->m_tracer.m_colorTransferName = potScript("ColorTransfer","");
	

	if( pScript->m_tracer.m_numLights > 8 ) {
		printf("Number of lights cannot be more than %d\n", MAX_LIGHTS );
		DELETE_SAFE( pScript );
		return NULL;
	}

	for(int j=0;j<pScript->m_tracer.m_numLights; ++j ){
		
		std::string secName;
		char buffer[20];
		
		sprintf(buffer,"%d/",j);
		secName = std::string("Light") + std::string(buffer);

		potScript.set_prefix(secName.c_str());

		pScript->m_lights[j].m_position[X] = potScript("Position",0.0,X);
		pScript->m_lights[j].m_position[Y] = potScript("Position",0.0,Y);
		pScript->m_lights[j].m_position[Z] = potScript("Position",0.0,Z);
		pScript->m_lights[j].m_position[W] = potScript("Position",0.0,W);

		pScript->m_lights[j].m_color[R] = potScript("Color",0.0,R);
		pScript->m_lights[j].m_color[G] = potScript("Color",0.0,G);
		pScript->m_lights[j].m_color[B] = potScript("Color",0.0,B);
		pScript->m_lights[j].m_color[A] = potScript("Color",0.0,A);

	}

	return pScript;
}
