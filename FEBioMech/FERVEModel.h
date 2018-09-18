#pragma once
#include "FECore/FEModel.h"
#include <FECore/tens4d.h>

class FEBCPrescribedDeformation;

//-----------------------------------------------------------------------------
// Class describing the RVE model.
// This is used by the homogenization code.
class FERVEModel : public FEModel
{
public:
	enum RVE_TYPE
	{
		DISPLACEMENT,		// prescribed displacement
		PERIODIC_LC,		// periodic, linear constraints
		PERIODIC_AL			// periodic, augmented Lagrangian (NOTE: obsolete, should probably delete)
	};

public:
	FERVEModel();
	~FERVEModel();

	//! one time initialization
	bool InitRVE(int rveType, const char* szbc);

	//! Return the initial volume (calculated in Init)
	double InitialVolume() const { return m_V0; }

	//! return current volume (calculated each time)
	double CurrentVolume();

	// scale the geometry
	void ScaleGeometry(double scale);

	//! see if node is boundary node
	bool IsBoundaryNode(int i) const { return (m_BN[i]==1); }

	//! Update the RVE (before it is solved)
	void Update(const mat3d& F);

	// copy from the master RVE
	void CopyFrom(FERVEModel& rve);

	//! Calculate the stress average
	mat3ds StressAverage(FEMaterialPoint& mp);

	//! Calculate the stiffness average
	tens4ds StiffnessAverage(FEMaterialPoint &mp);

protected:
	//! Calculate the initial volume
	void EvalInitialVolume();

	//! find the list of boundary nodes
	void FindBoundaryNodes(vector<int>& BN);

	//! Center the RVE
	void CenterRVE();

	bool PrepDisplacementBC(const FENodeSet& set);
	bool PrepPeriodicBC(const char* szbc);
	bool PrepPeriodicLC();

private:
	double			m_V0;				//!< initial volume
	int				m_bctype;			//!< RVE type
	FEBoundingBox	m_bb;				//!< bounding box of mesh
	vector<int>		m_BN;				//!< boundary node flags
};