//
//  FEFluidNormalTraction.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 8/22/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#ifndef FEFluidNormalTraction_hpp
#define FEFluidNormalTraction_hpp

#include "FECore/FESurfaceLoad.h"
#include <FECore/FESurfaceMap.h>

//-----------------------------------------------------------------------------
//! FEFluidNormalTraction is a fluid surface that has a normal
//! viscous traction prescribed on it.
//!
class FEFluidNormalTraction : public FESurfaceLoad
{
public:
    //! constructor
    FEFluidNormalTraction(FEModel* pfem);
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps) override;
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override {}
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R) override;
    
    //! Unpack surface element data
    void UnpackLM(FEElement& el, vector<int>& lm);
    
private:
    double			m_scale;	//!< magnitude of traction load
    FESurfaceMap	m_TC;		//!< traction boundary cards
    
    int		m_dofWX;
    int		m_dofWY;
    int		m_dofWZ;
    
    DECLARE_PARAMETER_LIST();
};

#endif /* FEFluidNormalTraction_hpp */