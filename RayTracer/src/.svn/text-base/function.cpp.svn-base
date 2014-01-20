#include "function.h"
#include "mldata.h"
#include "hamdata.h"

CMLData g_ML("ML",1);
CHAMData g_HAM("HAM",1);

std::map<std::string, CBaseData*> CFunctionSystem::m_functionNameMap;

void CFunctionSystem::InitFunctionSystem()
{
	m_functionNameMap["ML"] = &g_ML;
	m_functionNameMap["HAM"] = &g_HAM;
}

CBaseData * CFunctionSystem::GetFunctionByName( std::string name )
{
	if( m_functionNameMap.find( name ) == m_functionNameMap.end() ){
		printf("The function %s doesn't exist\n", name.c_str());
		return NULL;
	}

	return m_functionNameMap[ name ];
}
