#pragma once
#include "FECore/LinearSolver.h"
#include "BlockMatrix.h"
#include "PardisoSolver.h"

//-----------------------------------------------------------------------------
// This class implements solution strategies for solving linear systems by taking
// advantage of their block structure.
class BlockSolver : public LinearSolver
{
public:
	//! constructor
	BlockSolver();

	//! destructor
	~BlockSolver();

	//! Preprocess 
	bool PreProcess();

	//! Factor matrix
	bool Factor();

	//! Backsolve the linear system
	bool BackSolve(vector<double>& x, vector<double>& b);

	//! Clean up
	void Destroy();

	//! Create a sparse matrix
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype);

public:
	// Set the relative convergence tolerance
	void SetRelativeTolerance(double tol);

	// get the iteration count
	int GetIterations() const;

	// set the print level
	void SetPrintLevel(int n);

private:
	BlockMatrix*			m_pA;		//!< block matrices
	vector<PardisoSolver*>	m_solver;	//!< solvers for solving diagonal blocks

private:
	double	m_tol;			//!< convergence tolerance
	int		m_maxiter;		//!< max number of iterations
	int		m_iter;			//!< nr of iterations of last solve
	int		m_printLevel;	//!< set print level
};