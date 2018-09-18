#pragma once
#include "FEElasticMaterial.h"
#include "FERemodelingElasticMaterial.h"

//-----------------------------------------------------------------------------
//! This is a neo-Hookean material whose Young's modulus is evaluated from the density
//! according to the power-law relation proposed by Carter and Hayes for trabecular bone

class FECarterHayesOld : public FEElasticMaterial, public FERemodelingInterface
{
public:
	FECarterHayesOld(FEModel* pfem) : FEElasticMaterial(pfem) {}
	
public:
	double	m_c;	//!< c coefficient for calculation of Young's modulus
	double	m_g;	//!< gamma exponent for calculation of Young's modulus
	double	m_v;	//!< prescribed Poisson's ratio

public:
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;
	
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) override;
	
public: // --- remodeling interface ---

	//! calculate strain energy density at material point
	double StrainEnergy(FEMaterialPoint& pt) override;
	
	//! calculate tangent of strain energy density with mass density
	double Tangent_SE_Density(FEMaterialPoint& pt) override;
	
	//! calculate tangent of stress with mass density
	mat3ds Tangent_Stress_Density(FEMaterialPoint& pt) override;

	//! return Young's modulus
	double YoungModulus(double rhor) { return m_c*pow(rhor, m_g);}
	
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};