#ifndef __POT_SCRIPT_H__

#include "utils.h"
#include "tracer.h"

#include <string>
#include <map>

#define MAX_GRIDS 8

// For all the enumeration,
// the first element will always 
// be the default. So if an
// option is not given properly
// or spelling mistakes are made
// then the first element will be
// taken as default.

typedef enum PotExperimentType {
	POTEXPERIMENT_TRACED,
	POTEXPERIMENT_RANDOM,
}POTEXPERIMENTTYPE;

typedef enum PotGridType_s{	

	POTGRID_BCC,
	POTGRID_CC,
	POTGRID_HYBRID,
}POTGRIDTYPE;

typedef enum PotTracerMode_s{
	POTTRACERMODE_DVR,
	POTTRACERMODE_ISR,
}POTTRACERMODE;

typedef enum PotGridDataFormat_s{
	POTGRID_BYTE,
	POTGRID_FLOAT
}POTGRIDDATAFORMAT;

typedef enum PotGridDataType_s{
	
	POTGRID_FUNCTION,
	POTGRID_FILE,
	POTGRID_SAMPLED,

} POTGRIDDATATYPE;

typedef enum PotCameraType_s {

	POTCAMERA_PERSPECTIVE,
	POTCAMERA_ORTHOGRAPHIC

} POTCAMERATYPE;

typedef struct PotGridDataFiles_s{
	
	std::string			m_scalarFile;	
	std::string			m_normalXFile;
	std::string			m_normalYFile;
	std::string			m_normalZFile;
	POTGRIDDATAFORMAT	m_dataFormat;
}PotGridDataFiles;


typedef struct PotGrid_s{
	
	std::string			m_name;
	POTGRIDTYPE			m_type;
	POTGRIDTYPE			m_subType; //only useful for hybrid
	POTGRIDDATATYPE		m_dataType;
	POTGRIDDATATYPE		m_normalType;
	std::string			m_reconstructionFilter;
	std::string			m_normalReconstructionFilter;
	std::string     	m_normalEstimatorFilter;
	bool			m_bNormalEstimatorOnTheFly;
	
	std::string			m_dataFunctionName;
	vector3_t			m_dataFunctionMin;
	vector3_t			m_dataFunctionMax;
	int					m_dataFunctionGrid[3];

	elem_t				m_epsilon;
	elem_t				m_xScale;
	elem_t				m_yScale;
	elem_t				m_zScale;
	std::string			m_boundaryCond;

	PotGridDataFiles	m_dataFiles;

	std::string			m_gridPointDumpFile;

} PotGrid;

typedef struct PotOutput_s{
	std::string			m_traceDumpFile;
	std::string			m_image;
} PotOutput;

typedef struct PotCamera_s {
	POTCAMERATYPE		m_type;
	elem_t				m_FOV;
	elem_t				m_nearPlane;
	elem_t				m_farPlane;
	vector3_t			m_position;
	vector3_t			m_up;
	vector3_t			m_lookAt;
}PotCamera;

typedef struct PotTracer_s{

	POTTRACERMODE			m_tracerMode;
	double				m_isoSurface;
	double				m_isoSurfaceThreshold;
	int					m_imgWidth;
	int					m_imgHeight;

	bool				m_bPreview;
	bool				m_bShadowed;
	bool				m_bTwoSidedLighting;
	bool				m_bLighting;
	bool				m_bSuperSampling;

	int					m_numLights;
	
	elem_t				m_samplingRate;
	elem_t				m_shadowOffset;

	std::string			m_alphaTransferName;
	std::string			m_colorTransferName;
	
	vector4_t			m_shadowAmbient;
	vector4_t			m_ambient;
} PotTracer;

typedef struct PotLight_s {

	vector4_t			m_position;
	color_t				m_color;

} PotLight;

class CPotScript{
public:

	POTEXPERIMENTTYPE	m_experimentType;
	long int		m_seed;
	long int		m_randomSamples;
	
	vector3_t		m_randomBoundingMin;
	vector3_t		m_randomBoundingMax;
	
	int					m_numGrids; 
	int					m_motherGridIdx;
	

	PotGrid				m_grids[ MAX_GRIDS ];
	PotLight			m_lights[ MAX_LIGHTS ];

	PotOutput			m_output;
	PotCamera			m_camera;
	PotTracer			m_tracer;

	static std::map <std::string,POTGRIDTYPE>       m_potGridTypeMap;
	static std::map <std::string,POTGRIDDATATYPE>   m_potGridDataTypeMap;
	static std::map <std::string,POTCAMERATYPE>     m_potCameraTypeMap;
	static std::map <std::string,POTEXPERIMENTTYPE> m_potExperimentTypeMap;
	static std::map <std::string,POTGRIDDATAFORMAT> m_potGridDataFormatMap;
	static std::map <std::string, bool>             m_potOnOffMap;  
	static std::map <std::string,POTTRACERMODE>     m_potTracerModeMap;
public:
	
	static CPotScript* ParsePotScript( std::string fname );
	static void InitPotScriptSystem();
};


#endif
