#pragma once
#include "FEBiphasicSolute.h"

//-----------------------------------------------------------------------------
// This class implements a material that has a constant orthotropic diffusivity

class FEDiffConstOrtho :	public FESoluteDiffusivity
{
public:
	//! constructor
	FEDiffConstOrtho(FEModel* pfem);
	
	//! free diffusivity
	double Free_Diffusivity(FEMaterialPoint& pt) override;

	//! Tangent of free diffusivity with respect to concentration
	double Tangent_Free_Diffusivity_Concentration(FEMaterialPoint& mp, const int isol) override;
	
	//! diffusivity
	mat3ds Diffusivity(FEMaterialPoint& pt) override;
	
	//! Tangent of diffusivity with respect to strain
	tens4ds Tangent_Diffusivity_Strain(FEMaterialPoint& mp) override;
	
	//! Tangent of diffusivity with respect to concentration
	mat3ds Tangent_Diffusivity_Concentration(FEMaterialPoint& mp, const int isol) override;
	
	//! data checking
	bool Validate() override;
	
public:
	double	m_free_diff;	//!< free diffusivity
	double	m_diff[3];		//!< principal diffusivities
	
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};