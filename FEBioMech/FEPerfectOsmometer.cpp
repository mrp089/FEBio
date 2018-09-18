/*
 *  FEPerfectOsmometer.cpp
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 7/14/10.
 *
 */
#include "stdafx.h"
#include "FEPerfectOsmometer.h"
#include "FECore/FEModel.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEPerfectOsmometer, FEElasticMaterial)
	ADD_PARAMETER2(m_phiwr, FE_PARAM_DOUBLE, FE_RANGE_CLOSED(0.0, 1.0), "phiw0");
	ADD_PARAMETER2(m_iosm , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "iosm");
	ADD_PARAMETER2(m_bosm , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "bosm");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// FEPerfectOsmometer
//-----------------------------------------------------------------------------

bool FEPerfectOsmometer::Init()
{
	if (FEElasticMaterial::Init() == false) return false;
	
	m_Rgas = GetFEModel()->GetGlobalConstant("R");
	m_Tabs = GetFEModel()->GetGlobalConstant("T");
	
	if (m_Rgas <= 0) return MaterialError("A positive universal gas constant R must be defined in Globals section");
	if (m_Tabs <= 0) return MaterialError("A positive absolute temperature T must be defined in Globals section");

	return true;
}

//-----------------------------------------------------------------------------
mat3ds FEPerfectOsmometer::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// jacobian
	double J = pt.m_J;
	
	// calculate internal concentration in current configuration
	double iosm = m_iosm*m_phiwr/(J-1+m_phiwr);
	
	// calculate osmotic pressure
	double p = m_Rgas*m_Tabs*(iosm - m_bosm);
	
	// calculate T = -p*I
	mat3dd I(1.0);	// identity tensor
	mat3ds s = -p*I;
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEPerfectOsmometer::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// jacobian
	double J = pt.m_J;

	// calculate internal osmolarity in current configuration
	double iosm = m_iosm*m_phiwr/(J-1+m_phiwr);
	
	// calculate osmotic pressure
	double p = m_Rgas*m_Tabs*(iosm - m_bosm);
	
	// calculate derivative of osmotic pressure w.r.t. J
	double dp = -m_Rgas*m_Tabs*iosm/(J-1+m_phiwr);
	
	mat3dd I(1.0);	// Identity
	
	tens4ds I1 = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	
	// calculate tangent osmotic modulus
	tens4ds c = -J*dp*I1 + p*(2.0*I4 - I1);
	return c;
}
