#ifndef __FUNCTION_H__
#define __FUNCTION_H__

#include "data.h"
#include <map>

class CFunctionSystem {
public:
	static CBaseData * GetFunctionByName( std::string name );

	static std::map<std::string,CBaseData*> m_functionNameMap;

	static void InitFunctionSystem();
};

#endif
