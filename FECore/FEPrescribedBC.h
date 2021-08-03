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



#pragma once
#include "FEBoundaryCondition.h"
#include "FEDofList.h"

//-----------------------------------------------------------------------------
class FENodeSet;
class FESurface;

//-----------------------------------------------------------------------------
// base class for prescribed boundary conditions
class FECORE_API FEPrescribedBC : public FEBoundaryCondition
{
public:
	// constructor
	FEPrescribedBC(FEModel* pfem);

	// set the relative flag
	void SetRelativeFlag(bool br);

	// This function is called when the solver needs to know the 
	// prescribed dof values. The brel flag indicates wheter the total 
	// value is needed or the value with respect to the current nodal dof value
	void PrepStep(std::vector<double>& ui, bool brel = true) override;

	// serialization
	void Serialize(DumpStream& ar) override;

	// initialization
	bool Init() override;

	// activation
	void Activate() override;

	// deactivation
	void Deactivate() override;

	// Set the node list
	void SetNodeList(const FENodeList& nodeList);

	void Update() override;
    
    void Repair() override;

public:
	// Derived classes must implement this to set the dof list.
	// This function is called during initialization
	virtual bool SetDofList(FEDofList& dofs) = 0;

	//! Derived classes need to override this function.
	//! return the value for node i, dof j (i is index into nodeset, j is index into doflist)
	virtual void GetNodalValues(int nodelid, std::vector<double>& val) = 0;

private:
	bool	m_brelative;		//!< relative flag

private:
	FENodeList			m_nodeList;	//!< list of nodes to apply bc too
	std::vector<double>	m_rval;		//!< values used for relative BC

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Base class for prescribed BCs on a nodeset
class FECORE_API FEPrescribedNodeSet : public FEPrescribedBC
{
public:
	FEPrescribedNodeSet(FEModel* fem);

	void SetNodeSet(FENodeSet* nodeSet);

	const FENodeSet* GetNodeSet();

	bool Init() override;

	void Activate() override;

private:
	FENodeSet*	m_nodeSet;

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Base class for prescribed BCs on a surface
class FECORE_API FEPrescribedSurface : public FEPrescribedBC
{
public:
	FEPrescribedSurface(FEModel* fem);

	void SetSurface(FESurface* surface);

	const FESurface* GetSurface();

	bool Init() override;

	void Activate() override;

private:
	FESurface*	m_surface;

	DECLARE_FECORE_CLASS();
};
