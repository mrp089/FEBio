//
//  FENewtonianFluid.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 9/6/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "FENewtonianFluid.h"
#include "FEFluid.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FENewtonianFluid, FEViscousFluid)
    ADD_PARAMETER2(m_kappa, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "kappa");
    ADD_PARAMETER2(m_mu, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FENewtonianFluid::FENewtonianFluid(FEModel* pfem) : FEViscousFluid(pfem)
{
    m_kappa = 0;
    m_mu = 0;
}

//-----------------------------------------------------------------------------
//! viscous stress
mat3ds FENewtonianFluid::Stress(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    mat3ds D = vt.RateOfDeformation();
    
    mat3ds s = mat3dd(1.0)*(D.tr()*(m_kappa - 2.*m_mu/3.)) + D*(2*m_mu);
    
    return s;
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to strain J
mat3ds FENewtonianFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    return mat3ds(0,0,0,0,0,0);
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to rate of deformation tensor D
tens4ds FENewtonianFluid::Tangent_RateOfDeformation(FEMaterialPoint& mp)
{
    mat3dd I(1.0);
    tens4ds c = dyad1s(I)*(m_kappa - 2.*m_mu/3.) + dyad4s(I)*(2*m_mu);
    return c;
}

//-----------------------------------------------------------------------------
//! dynamic shear viscosity
double FENewtonianFluid::ShearViscosity(FEMaterialPoint& mp)
{
    return m_mu;
}

//-----------------------------------------------------------------------------
//! bulke viscosity
double FENewtonianFluid::BulkViscosity(FEMaterialPoint& mp)
{
    return m_kappa;
}