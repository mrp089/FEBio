//
//  FEDamageMaterial.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/18/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEDamageMaterial.h"
#include "FEDamageCriterion.h"
#include "FEDamageCDF.h"
#include "FEUncoupledMaterial.h"
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageMaterial::FEDamageMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	// set material properties
	AddProperty(&m_pBase, "elastic"  );
	AddProperty(&m_pDamg, "damage"   );
	AddProperty(&m_pCrit, "criterion");
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FEDamageMaterial::Init()
{
    FEUncoupledMaterial* m_pMat = dynamic_cast<FEUncoupledMaterial*>((FEElasticMaterial*)m_pBase);
    if (m_pMat != nullptr)
        return MaterialError("Elastic material should not be of type uncoupled");
    
	return FEElasticMaterial::Init();
}

//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FEDamageMaterial::Stress(FEMaterialPoint& pt)
{
    // evaluate the damage
    double d = Damage(pt);
    
    // evaluate the stress
    mat3ds s = m_pBase->Stress(pt);
    
    // return damaged stress
    return s*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate tangent stiffness at material point
tens4dss FEDamageMaterial::Tangent(FEMaterialPoint& pt)
{
    // evaluate the damage
    double d = Damage(pt);
    
    // evaluate the tangent
    tens4dss c = m_pBase->Tangent(pt);
    
    // return damaged tangent
    return c*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEDamageMaterial::StrainEnergyDensity(FEMaterialPoint& pt)
{
    // evaluate the damage
    double d = Damage(pt);
    
    // evaluate the strain energy density
    double sed = m_pBase->StrainEnergyDensity(pt);

    // return damaged sed
    return sed*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate damage at material point
double FEDamageMaterial::Damage(FEMaterialPoint& pt)
{
    // get the damage material point data
    FEDamageMaterialPoint& pd = *pt.ExtractData<FEDamageMaterialPoint>();

    // evaluate the trial value of the damage criterion
    // this must be done before evaluating the damage
    pd.m_Etrial = m_pCrit->DamageCriterion(pt);
    
    // evaluate and set the damage
    double d = m_pDamg->Damage(pt);
    pd.m_D = d;
    
    return d;
}

//-----------------------------------------------------------------------------
void FEDamageMaterial::SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp)
{
	FEElasticMaterial::SetLocalCoordinateSystem(el, n, mp);
	m_pBase->SetLocalCoordinateSystem(el, n, mp);
}
