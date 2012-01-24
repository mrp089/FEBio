#pragma once
#include "FEBiphasicSolute.h"

//-----------------------------------------------------------------------------
//! Base class for solute materials.

class FESolute : public FEMaterial
{
public:
	FESolute();

public:
	void Init();
	
	//! solute density
	double Density() { return m_rhoT; }
	
	//! solute molecular weight
	double MolarMass() { return m_M; }
	
	//! solute valence
	double ChargeNumber() { return m_z; }
	
	//! Serialization
	void Serialize(DumpFile& ar);

	//! set solute ID
	void SetSoluteID(const int ID) {m_ID = ID;}
	
	//! get solute ID
	int GetSoluteID() {return m_ID;}
	
private:
	int						m_ID;		//!< solute ID
	
public:
	double					m_rhoT;		//!< true solute density
	double					m_M;		//!< solute molecular weight
	int						m_z;		//!< charge number of solute
	FESoluteDiffusivity*	m_pDiff;	//!< pointer to diffusivity material
	FESoluteSolubility*		m_pSolub;	//!< pointer to solubility material
	
	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! Global solute data
//! This structure uniquely identifies a solute in multiphasic problems

struct FESoluteData {
	int					m_nID;			//!< solute ID
	char				m_szname[128];	//!< solute name
	double				m_z;			//!< solute charge number
};
