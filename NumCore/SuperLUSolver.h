#pragma once
#include "FECore/LinearSolver.h"
#include "CompactUnSymmMatrix.h"

//-----------------------------------------------------------------------------
//! Implements a wrapper class for the SuperLU library

//! This solver can only be used on systems where it is available.
//! This solver also uses some of the BLAS routines so this package also needs
//! to be available on the system. Although SuperLU comes with a stripped down
//! version of BLAS.

#ifdef SUPERLU
		#include "slu_ddefs.h"
#endif

class SuperLUSolver : public LinearSolver
{
public:
	SuperLUSolver();
	bool PreProcess();
	bool Factor();
	bool BackSolve(vector<double>& x, vector<double>& b);
	void Destroy();
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype);

	void print_cnorm(bool b) { m_bcond = b; }

#ifdef SUPERLU
protected:
	double norm(SparseMatrix& K); // calculates the 1-norm of the matrix A
#endif

private:

	bool m_bsymm;	// use symmetric mode or not
	bool m_balloc;
	bool m_bfact;
	bool m_bcond;	// calculate condition numbers

	CCSSparseMatrix*	m_pA;

#ifdef SUPERLU

	SuperMatrix A, L, U, B, X;
	vector<int>	perm_c;
	vector<int>	perm_r;
	vector<int>	etree;

	superlu_options_t	options;
	SuperLUStat_t	stat;
	mem_usage_t	mem_usage;

	double	rpg, rcond;
	double	ferr, berr;
	int		info;
	char	equed[1];

#endif // SUPERLU
};
