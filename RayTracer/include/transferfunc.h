#ifndef __TRANSFER_FUNC__
#define __TRANSFER_FUNC__

#include "tracer.h"

#include <string>
#include <map>

class CTransferFunction{
private:

	static std::map<std::string,fnAlphaTransfer> m_alphaMap;
	static std::map<std::string,fnColorTransfer> m_colorMap;

public :

	static bool SetAlphaTransferByName( std::string name );
	static bool SetColorTransferByName( std::string name );
	static void InitTransferFunctionSystem();
};

#endif
