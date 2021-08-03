/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FEPermExpIso.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEPermExpIso, FEHydraulicPermeability)
    ADD_PARAMETER(m_perm , FE_RANGE_GREATER_OR_EQUAL(0.0), "perm" );
    ADD_PARAMETER(m_M    , FE_RANGE_GREATER_OR_EQUAL(0.0), "M"    );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEPermExpIso::FEPermExpIso(FEModel* pfem) : FEHydraulicPermeability(pfem)
{
    m_perm = 1;
    m_M = 0;
}

//-----------------------------------------------------------------------------
//! Permeability tensor.
mat3ds FEPermExpIso::Permeability(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
    FEBiphasicMaterialPoint* pt = mp.ExtractData<FEBiphasicMaterialPoint>();
    FEBiphasicFSIMaterialPoint* bpt = mp.ExtractData<FEBiphasicFSIMaterialPoint>();
    
    // relative volume
    double J = et.m_J;
    // referential solid volume fraction
    double phi0 = 0.0;
    if(pt)
        phi0 = pt->m_phi0;
    else if (bpt)
        phi0 = bpt->m_phi0;
    
    // --- strain-dependent isotropic permeability ---
    double k0 = m_perm*exp(m_M*(J-1)/(J-phi0));
    
    return mat3dd(k0);
}

//-----------------------------------------------------------------------------
//! Tangent of permeability
tens4dmm FEPermExpIso::Tangent_Permeability_Strain(FEMaterialPoint &mp)
{
    FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
    FEBiphasicMaterialPoint* pt = mp.ExtractData<FEBiphasicMaterialPoint>();
    FEBiphasicFSIMaterialPoint* bpt = mp.ExtractData<FEBiphasicFSIMaterialPoint>();
    
    // relative volume
    double J = et.m_J;
    
    // referential solid volume fraction
    double phi0 = 0.0;
    if(pt)
        phi0 = pt->m_phi0;
    else if (bpt)
        phi0 = bpt->m_phi0;

    mat3dd I(1);    // Identity
    
    double k0 = m_perm*exp(m_M*(J-1)/(J-phi0));
    double k0prime = m_M*(1-phi0)/pow(J-phi0,2)*k0;

    mat3ds k0hat = I*(k0 + J*k0prime);
    
    return dyad1mm(I,k0hat)-dyad4s(I)*(2*k0);
}
