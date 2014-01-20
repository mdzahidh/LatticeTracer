#ifndef __HAMDATA_H__
#define __HAMDATA_H__

#include <math.h>

#include "data.h"

class CHAMData : public CBaseData {
private:
	double m_alpha;
	double m_freq;
			
	double m_2pifm;
	double m_2pifmalpha;
	double m_beta;
public:
	CHAMData( const char *name, float size ) : CBaseData( std::string(name)), m_alpha(0.25), m_freq(6.0), m_beta(2.0)
	{
		m_2pifm = 2 * PI * m_freq;
		m_2pifmalpha = m_2pifm * m_alpha;
		 
		m_pVolume = new CVolume( 3,3,3, size,size,size );
	}
	
	virtual ~CHAMData()
	{
		DELETE_SAFE( m_pVolume );
	}
	
	virtual double Get( const vector3_t x ) const 
	{	
		double rSq = (x[X] * x[X]) + (x[Y] * x[Y]) + (x[Z]*x[Z]);
		
		if (rSq < EPSILON){			
			return 0;
		}
				
		double r = sqrt(rSq);
					   
		return m_beta * r - m_alpha * cos( (m_2pifm * x[Z]) / r ); 
	}
	
	virtual void GetNormal( const vector3_t x, vector3_t N ) const
	{
		double rSq = (x[X] * x[X]) + (x[Y] * x[Y]) + (x[Z]*x[Z]);
		
		if (rSq < EPSILON){
			N[X] = N[Y] = N[Z] = 0;
			return;
		}
				
		double r = sqrt(rSq);		
		double oneOverR = 1.0/r;
		
		double twicePiFmZ = m_2pifm * x[Z];
		double twicePiFmZOverR = twicePiFmZ / r;
		double zOverRCube =x[Z] / (rSq * r);
		
		double sinPart = sin( twicePiFmZOverR );
		double scale = (oneOverR * m_beta - m_2pifmalpha * sinPart * zOverRCube);
		
		N[X] = x[X] * scale;
		N[Y] = x[Y] * scale;
		N[Z] = x[Z] * oneOverR * m_beta + m_2pifmalpha * sinPart * ( oneOverR - zOverRCube * x[Z]);
	}	
};

#endif
