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
#include "FELinearConstraint.h"
#include "FEMesh.h"
#include "FEModel.h"
#include "DumpStream.h"

//-----------------------------------------------------------------------------
FELinearConstraint::FELinearConstraint() : FEModelComponent(nullptr)
{
	m_off = 0.0;
}

//-----------------------------------------------------------------------------
FELinearConstraint::FELinearConstraint(FEModel* pfem) : FEModelComponent(pfem) 
{
	m_off = 0.0;
}

//-----------------------------------------------------------------------------
FELinearConstraint::FELinearConstraint(const FELinearConstraint& LC) : FEModelComponent(LC.GetFEModel())
{
	m_parentDof = LC.m_parentDof;
	m_off = LC.m_off;
	int n = (int) LC.m_childDof.size();
	vector<DOF>::const_iterator it = LC.m_childDof.begin();
	for (int i=0; i<n; ++i, ++it) m_childDof.push_back(*it);
}

//-----------------------------------------------------------------------------
void FELinearConstraint::CopyFrom(const FELinearConstraint& LC)
{
	m_parentDof = LC.m_parentDof;
	m_off = LC.m_off;
	int n = (int)LC.m_childDof.size();
	vector<DOF>::const_iterator it = LC.m_childDof.begin();
	for (int i = 0; i<n; ++i, ++it) m_childDof.push_back(*it);
}

//-----------------------------------------------------------------------------
void FELinearConstraint::SetParentDof(int dof, int node)
{
	m_parentDof.dof = dof;
	m_parentDof.node = node;
}

//-----------------------------------------------------------------------------
void FELinearConstraint::AddChildDof(int dof, int node, double v)
{
	DOF d;
	d.dof = dof;
	d.node = node;
	d.val = v;
	m_childDof.push_back(d);
}

//-----------------------------------------------------------------------------
// Initialization.
// Make sure the parent dof does not appear as a child dof
bool FELinearConstraint::Init()
{
	int n = (int)m_childDof.size();
	for (int i=0; i<n; ++i)
	{
		DOF& childNode = m_childDof[i];
		if ((childNode.node == m_parentDof.node) && (childNode.dof == m_parentDof.dof)) return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
// This is called during model activation (i.e. at the start of an analysis step)
// The parent dof is fixed in order to make sure that they are not assigned an equation number.
void FELinearConstraint::Activate()
{
	FEModelComponent::Activate();
	FEMesh& mesh = GetFEModel()->GetMesh();

	// we need the parent node to be fixed so that no equation is allocated
	FENode& node = mesh.Node(m_parentDof.node);
	node.set_bc(m_parentDof.dof, DOF_FIXED);
}

//-----------------------------------------------------------------------------
void FELinearConstraint::Deactivate()
{
	FEModelComponent::Deactivate();
	FEMesh& mesh = GetFEModel()->GetMesh();

	FENode& node = mesh.Node(m_parentDof.node);
	node.set_bc(m_parentDof.dof, DOF_OPEN);
}

//-----------------------------------------------------------------------------
void FELinearConstraint::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	if (ar.IsSaving())
	{
		ar.write(&m_parentDof, sizeof(DOF), 1);
		int n = (int)m_childDof.size();
		ar << n;
		vector<DOF>::iterator it = m_childDof.begin();
		for (int i=0; i<n; ++i, ++it) ar << it->val << it->node << it->dof;
	}
	else
	{
		m_childDof.clear();
		ar.read(&m_parentDof, sizeof(DOF), 1);
		int n;
		ar >> n;
		for (int i=0; i<n; ++i)
		{
			DOF dof;
			ar >> dof.val >> dof.node >> dof.dof;
			m_childDof.push_back(dof);
		}
	}
}
