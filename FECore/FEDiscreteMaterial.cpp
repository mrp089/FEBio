#include "stdafx.h"
#include "FEDiscreteMaterial.h"

//=============================================================================
FEMaterialPoint* FEDiscreteMaterialPoint::Copy()
{
	FEDiscreteMaterialPoint* pt = new FEDiscreteMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//=============================================================================
FEDiscreteMaterial::FEDiscreteMaterial(FEModel* pfem) : FEMaterial(pfem) {}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEDiscreteMaterial::CreateMaterialPointData() { return new FEDiscreteMaterialPoint; }