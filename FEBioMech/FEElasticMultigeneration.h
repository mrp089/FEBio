#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material defining a single generation of a multi-generation material
class FEGenerationMaterial : public FEElasticMaterial
{
public:
	FEGenerationMaterial(FEModel* pfem);

	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;
		
	//! calculate tangent stiffness at material point
	tens4dss Tangent(FEMaterialPoint& pt) override;

	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
    // returns a pointer to a new material point object
    FEMaterialPoint* CreateMaterialPointData() override {
        return m_pMat->CreateMaterialPointData();
    }
    
    //! Get the elastic component
    FEElasticMaterial* GetElasticMaterial() override { return m_pMat; }
    
public:
	double	btime;	//!< generation birth time

	FEPropertyT<FEElasticMaterial>	m_pMat;	//!< pointer to elastic material

	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
// forward declaration of material class
class FEElasticMultigeneration;

//-----------------------------------------------------------------------------
//! Multigenerational material point.
//! First generation exists at t=0. Second, third, etc. generations appear at t>0.
//! This material point stores the inverse of the relative deformation gradient of
//! second, third, etc. generations.  These relate the reference configuration of 
//! each generation relative to the first generation.

class FEMultigenerationMaterialPoint : public FEMaterialPoint
{
public:
    FEMultigenerationMaterialPoint();
		
	FEMaterialPoint* Copy();

	//! Add a child material point
	void AddMaterialPoint(FEMaterialPoint* pt);
		
	//! data serialization
	void Serialize(DumpStream& ar);

	void Init();

	void Update(const FETimeInfo& timeInfo);

    FEMaterialPoint* GetPointData(int i) { return m_mp[i]; }
    
public:
	// multigenerational material data
    vector<FEMaterialPoint*>    m_mp;   //!< material point data for multigeneration components
	double	m_tgen;		//!< last generation time
    int     m_ngen;     //!< number of active generations
	FEElasticMultigeneration*	m_pmat;
};

//-----------------------------------------------------------------------------
//! Multigenerational solid

class FEElasticMultigeneration : public FEElasticMaterial
{
public:
	FEElasticMultigeneration(FEModel* pfem);
		
	// returns a pointer to a new material point object
    FEMaterialPoint* CreateMaterialPointData();

    // return number of materials
    int Materials() { return (int)m_MG.size(); }
    
    // return a generation material component
    FEGenerationMaterial* GetMaterial(int i) { return m_MG[i]; }
    
	void AddMaterial(FEElasticMaterial* pmat);
	
public:
    //! Set the local coordinate system for a material point (overridden from FEMaterial)
    void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp);
    
public:
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);
		
	//! calculate tangent stiffness at material point
	tens4dss Tangent(FEMaterialPoint& pt);
		
	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt);
    
	int CheckGeneration(const double t);

public:
	FEVecPropertyT<FEGenerationMaterial>	m_MG;		//!< multigeneration data

	// declare the parameter list
//	DECLARE_PARAMETER_LIST();
};
