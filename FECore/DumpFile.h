// DumpFile.h: interface for the DumpFile class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ARCHIVE_H__B95A81B1_BBFB_46E5_B9B3_7675ED8A6029__INCLUDED_)
#define AFX_ARCHIVE_H__B95A81B1_BBFB_46E5_B9B3_7675ED8A6029__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>
#include <string.h>
#include "vec3d.h"
#include "mat3d.h"
#include "quatd.h"
#include <vector>

class FEModel;

//-----------------------------------------------------------------------------
//! Class for serializing data to a binary archive.

//! This class is used to read data from or write
//! data to a binary file. The class defines several operators to 
//! simplify in- and output.
//! \sa FEM::Serialize()

class DumpFile  
{
public:
	// This class is thrown when an error occurs reading the dumpfile
	class ReadError{};

public:
	DumpFile(FEModel* pfem);
	virtual ~DumpFile();

	//! Open archive for reading
	bool Open(const char* szfile);

	//! Open archive for writing
	bool Create(const char* szfile);

	//! Open archive for appending
	bool Append(const char* szfile);

	//! Close archive
	void Close();

	//! Check mode
	bool IsSaving() { return m_bsave; }

	//! See if the archive is valid
	bool IsValid() { return (m_fp != 0); }

	//! Flush the archive
	void Flush() { fflush(m_fp); }

	//@{ 
	//! output operators
	DumpFile& operator << (const char* sz);

	DumpFile& operator << (char* sz);

	DumpFile& operator << (const double a[3][3]);

	template <class T> DumpFile& operator << (const T& o) { write(&o, sizeof(T), 1); return (*this); }

	template <class T> DumpFile& operator << (std::vector<T>& v)
	{
		int n = v.size();
		write(&n, sizeof(int), 1);
		if (n>0) write((T*) &v[0], sizeof(T), v.size());
		return (*this);
	}

	DumpFile& operator << (std::vector<bool>& v);
	//@}


	//@{
	//! input operators
	DumpFile& operator >> (char* sz);

	DumpFile& operator >> (double a[3][3]);

	template <class T> DumpFile& operator >> (T& o) { read(&o, sizeof(T), 1); return (*this); }

	template <class T> DumpFile& operator >> (std::vector<T>& v)
	{
		int n;
		read(&n, sizeof(int), 1);
		if (n>0)
		{
			v.resize(n);
			read((T*) &v[0], sizeof(T), n);
		}
		else v.clear();
		return (*this);
	}

	DumpFile& operator >> (std::vector<bool>& v);
	//@}


	//! write buffer to archive
	size_t write(const void* pd, size_t size, size_t count);

	//! read buffer from archive
	size_t read(void* pd, size_t size, size_t count);

	//! get FEM model
	FEModel* GetFEModel() { return m_pfem; }

	//! get the current index
	int GetDataIndex() const { return m_nindex; }

protected:
	FILE*		m_fp;		//!< The actual file pointer
	FEModel*	m_pfem;		//!< FEM data that will be serialized
	bool		m_bsave;	//!< Save flag
	int			m_nindex;	//!< file index (gives amount of bytes written or read in so far)
};

#endif // !defined(AFX_ARCHIVE_H__B95A81B1_BBFB_46E5_B9B3_7675ED8A6029__INCLUDED_)
