#ifndef __HYBRID_DATA_H__
#define __HYBRID_DATA_H__

#include "data.h"
#include "utils.h"


class CHybridData : public CBaseData {
private:
	CBaseData *m_pFunction;
	CBaseData *m_pSampled;		
public:
	CHybridData( const char* name,  CBaseData *pFunction, CBaseData *pSampled ) : 
		CBaseData( std::string(name) ), 
		m_pFunction( pFunction ),
		m_pSampled( pSampled )
		{		
			//Hacking
			m_pVolume = new CVolume(1,1,1,1,1,1);
			*m_pVolume = *pFunction->GetVolume();	
		}
	virtual ~CHybridData()
	{
		// We will only destroy the sampled grid
		// because funciton grid is static in the sytem
		// and we only have pointer
		
		DELETE_SAFE( m_pSampled );
	}
		
	// We will get the data value from the function grid
	virtual double Get( const vector3_t p ) const { 		
		return m_pFunction->Get( p ); 
	}
	
	// wheares we will get the normal from the sample grid
	virtual void  GetNormal( const vector3_t p, vector3_t N ) const {
		m_pSampled->GetNormal( p, N );	
	}
};

#endif
