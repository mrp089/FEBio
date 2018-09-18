/*
 *  FESphericalFiberDistribution.h
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 12/26/11.
 *  Copyright 2011 Columbia University. All rights reserved.
 *
 */

#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material class for the spherical fiber distribution
//!

class FESphericalFiberDistribution : public FEElasticMaterial
{
public:
	FESphericalFiberDistribution(FEModel* pfem);
	
	//! Cauchy stress
	virtual mat3ds Stress(FEMaterialPoint& mp) override;
	
	// Spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp) override;
	
	// Strain energy density
	virtual double StrainEnergyDensity(FEMaterialPoint& mp) override;
	
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
	
public:
	double	m_beta;		// power in power-law relation
	double	m_ksi;		// coefficient in power-law relation
	double	m_alpha;	// coefficient of exponential argument
};