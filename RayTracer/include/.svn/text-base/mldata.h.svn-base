#ifndef __MLDATA_H__
#define __MLDATA_H__

#include <math.h>
#include "data.h"

/***
This is the analytic form of Marchner-Lobb Dataset
***/

class CMLData : public CBaseData {
private:
	double	m_alpha;
	double	m_freq;
public:

	CMLData( const char *name, float size ) : CBaseData(std::string(name)), m_alpha(0.25), m_freq(6)
	{
		m_pVolume = new CVolume(3,3,3,size,size,size);
	}

	virtual ~CMLData()
	{
		DELETE_SAFE( m_pVolume );
	}

	virtual double Get( const vector3_t x) const
	{
		vector3_t p;
		
		p[X] = x[X] / m_pVolume->m_HXSize;
		p[Y] = x[Y] / m_pVolume->m_HYSize;
		p[Z] = x[Z] / m_pVolume->m_HZSize;

		double r = sqrt(p[X] * p[X] + p[Y] * p[Y]);
		double lambda = cos( PI_OVER_2 * r );
		
		double ret = (1 - sin( PI_OVER_2 * p[Z] ) + m_alpha * ( 1 + cos( PI_2 * m_freq * lambda ) ) ) / (2*(1+m_alpha));

		return ret;
	}

	virtual void GetNormal( const vector3_t x, vector3_t N ) const 
	{
		vector3_t p;
		
		p[X] = x[X] / m_pVolume->m_HXSize;
		p[Y] = x[Y] / m_pVolume->m_HYSize;
		p[Z] = x[Z] / m_pVolume->m_HZSize;

		
		double r = sqrt(p[X] * p[X] + p[Y] * p[Y]);
		double coeff;

		if( (r > EPSILON) ){
			coeff = (m_alpha * sin(PI_2 * m_freq * cos(PI_OVER_2 * r)) * PI_SQUARE * m_freq * sin( PI_OVER_2 * r)) / (2*r*(1+m_alpha));
		}
		else{
			coeff = (PI_CUBE * m_alpha * m_freq * sin( PI * m_freq ) * cos( PI * m_freq )) / ( 2 * (1 + m_alpha) );
		}

		N[X] = coeff * p[X];
		N[Y] = coeff * p[Y];
		N[Z] = - (cos( PI_OVER_2 * p[Z] ) * PI) / (4*(1+m_alpha));

	}

};

#endif
