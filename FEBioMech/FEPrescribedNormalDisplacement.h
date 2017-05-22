#pragma once
#include <FECore/BC.h>

class FEPrescribedNormalDisplacement : public FEPrescribedBC
{
	struct NODE
	{
		int		nodeId;		// node ID
		vec3d	normal;		// initial normal at node
	};

public:
	// constructor
	FEPrescribedNormalDisplacement(FEModel* fem);

	// initialization
	bool Init();

	// activation
	void Activate();

	// deactivation
	void Deactivate();

public:
	// assign a node set to the prescribed BC
	void AddNodes(const FESurface& surf);

	// This function is called when the solver needs to know the 
	// prescribed dof values. The brel flag indicates wheter the total 
	// value is needed or the value with respect to the current nodal dof value
	void PrepStep(std::vector<double>& ui, bool brel = true);

	// This is called during nodal update and should be used to enforce the 
	// nodal degrees of freedoms
	void Update();

	// copy data from another class
	void CopyFrom(FEPrescribedBC* pbc);

private:
	vector<NODE>	m_node;

	double	m_scale;

	DECLARE_PARAMETER_LIST();
};