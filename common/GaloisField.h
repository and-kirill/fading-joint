/* 
 * This file is part of the faing-joint distribution
 * https://github.com/and-kirill/fading-joint.
 * Copyright (c) 2019 Alexey Frolov al.frolov@skoltech.ru.
 * 
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef GALOIS_FIELD_H
#define GALOIS_FIELD_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ====================== Template helpers ============================ */
template<unsigned int N> class Binary
{
public:
	enum
	{
		value = (N % 10) + 2 * Binary<N / 10>::value 
	};
};

template<> class Binary<0>
{
public:
	enum
	{ 
		value = 0 
	};
};

/* ====================== Galois field ================================ */
template<unsigned int iPower, unsigned int iPrimitivePolinomial>
class GaloisField
{
public:
	static GaloisField& instance()
	{
		static GaloisField<iPower, iPrimitivePolinomial> gf;
		return gf;
	}
	
	~GaloisField()
	{
		if(m_pDegreeToElementTable)
		{
			delete[] m_pDegreeToElementTable;
		}
		if(m_pElementToDegreeTable)
		{
			delete[] m_pElementToDegreeTable;
		}
	}

	
	unsigned char getPower() const
	{
		return iPower;
	}

	unsigned int getCardinality() const
	{
		return m_iCardinality;
	}

	int getElementByDegree(int iDegree) const
	{
		int iModulo = (int)(m_iCardinality - 1);
		iDegree %= iModulo;
		if (iDegree < 0)
		{
			iDegree += iModulo;
		}
		return m_pDegreeToElementTable[iDegree];
	}

	int getDegreeByElement(unsigned int iElement) const
	{
		return m_pElementToDegreeTable[moduloCast(iElement)];
	}

	unsigned int moduloCast(unsigned int iElement) const
	{
		int iShiftedModulo = iElement;
		int iShift = 0;

		while(iElement > (m_iCardinality - 1))
		{
			iShift = (int)(log((float)iElement)/log((float)2)) - m_iPower;
			iShiftedModulo = m_iPrimitivePolinomial << iShift;
			iElement ^= iShiftedModulo;
		}
		return iElement;
	}

private:
	void initializeTables()
	{
		unsigned int iPolinomial = Binary<10>::value; // alpha
	
		m_pDegreeToElementTable[0] = 1; 
		m_pElementToDegreeTable[0] = -1;

		m_pDegreeToElementTable[1] = Binary<10>::value; // alpha
		m_pElementToDegreeTable[1] = 0;

		m_pElementToDegreeTable[Binary<10>::value] = 1;



		for(unsigned int iDegree = 2; iDegree < (m_iCardinality - 1); ++iDegree)
		{
			iPolinomial <<= 1;
			if(iPolinomial >= (m_iCardinality - 1))
			{
				iPolinomial ^= m_iPrimitivePolinomial;
			}
			m_pDegreeToElementTable[iDegree] = iPolinomial;
			m_pElementToDegreeTable[iPolinomial] = iDegree;
		}

		m_pDegreeToElementTable[m_iCardinality - 1] = 1;
	}

	GaloisField()
		: m_iPower(iPower)
		, m_iCardinality(1 << iPower)
		, m_iPrimitivePolinomial(iPrimitivePolinomial)
		, m_pDegreeToElementTable(0)
		, m_pElementToDegreeTable(0)
	{
		m_pDegreeToElementTable = new int[m_iCardinality];
		m_pElementToDegreeTable = new int[m_iCardinality];
		initializeTables();
	}

	GaloisField(const GaloisField& gf);
private:
	unsigned char m_iPower;
	unsigned int  m_iCardinality;
	unsigned int  m_iPrimitivePolinomial;
	int*		  m_pDegreeToElementTable;
	int*		  m_pElementToDegreeTable;
};

// Spesialization for the binary case
template<unsigned int iPrimitivePolinomial>
class GaloisField<1, iPrimitivePolinomial>
{
public:
	static GaloisField& instance()
	{
		static GaloisField<1, iPrimitivePolinomial> gf;
		return gf;
	}
	
	~GaloisField()
	{
		if(m_pDegreeToElementTable)
		{
			delete[] m_pDegreeToElementTable;
		}
		if(m_pElementToDegreeTable)
		{
			delete[] m_pElementToDegreeTable;
		}
	}

	
	unsigned char getPower() const
	{
		return 1;
	}

	unsigned int getCardinality() const
	{
		return m_iCardinality;
	}

	int getElementByDegree(int iDegree) const
	{
		return m_pDegreeToElementTable[0];
	}

	int getDegreeByElement(unsigned int iElement) const
	{
		return m_pElementToDegreeTable[moduloCast(iElement)];
	}

	unsigned int moduloCast(unsigned int iElement) const
	{
		int iModulo = (int)(m_iCardinality);
		iElement %= 2;
		if (iElement < 0)
		{
			iElement += iModulo;
		}
		return iElement;
	}

private:
	void initializeTables()
	{
		m_pDegreeToElementTable[0] = 1; 
		m_pElementToDegreeTable[0] = -1;
		m_pElementToDegreeTable[1] = 0;
	}

	GaloisField()
		: m_iCardinality(2)
		, m_pDegreeToElementTable(0)
		, m_pElementToDegreeTable(0)
	{
		m_pDegreeToElementTable = new int[m_iCardinality];
		m_pElementToDegreeTable = new int[m_iCardinality];
		initializeTables();
	}

	GaloisField(const GaloisField& gf);
private:
	unsigned int  m_iCardinality;
	int*		  m_pDegreeToElementTable;
	int*		  m_pElementToDegreeTable;
};


/* =========================== Galois field element ===========================*/
template<typename T>
class GaloisFieldElement
{
public:
	GaloisFieldElement()
		: m_gf(T::instance())
		, m_iElement(0)
		, m_iDegree(-1)
		, m_bErased(false)
	{
	}
	
	GaloisFieldElement(unsigned int iElement)
		: m_gf(T::instance())
		, m_iElement(-1)
		, m_iDegree(-1)
		, m_bErased(false)
	{
		m_iElement = m_gf.moduloCast(iElement);
		m_iDegree = m_gf.getDegreeByElement(m_iElement);
	}

	GaloisFieldElement(const GaloisFieldElement& element)
		: m_gf(element.m_gf)
		, m_iElement(element.m_iElement)
		, m_iDegree(element.m_iDegree)
		, m_bErased(element.m_bErased)
	{
	}

	~GaloisFieldElement()
	{
	}

	const T& getGaloisField() const
	{
		return m_gf;
	}
	
	int getElement() const
	{
		return m_iElement;
	}
	
	int getDegree() const
	{
		return m_iDegree;
	}

	const GaloisFieldElement& operator=(const GaloisFieldElement& element)
	{
		m_iElement = element.m_iElement;
		m_iDegree = element.m_iDegree;
		m_bErased = element.m_bErased;
		return *this;
	}

	const GaloisFieldElement& operator=(unsigned int iElement)
	{
		m_iElement = m_gf.moduloCast(iElement);
		m_iDegree = m_gf.getDegreeByElement(m_iElement);
		m_bErased = false;
		return *this;
	}

	bool operator==(const GaloisFieldElement& element)  const
	{
		return (m_iElement == element.m_iElement);
	}

	GaloisFieldElement operator+(const GaloisFieldElement& element) const
	{
		GaloisFieldElement result(m_gf.moduloCast(m_iElement ^ element.m_iElement));
		return result;
	}

	GaloisFieldElement operator+(unsigned int iElement) const
	{
		GaloisFieldElement result(m_gf.moduloCast(m_iElement ^ iElement));
		return result;
	}

	const GaloisFieldElement& operator+=(const GaloisFieldElement& element)
	{
		m_iElement = m_gf.moduloCast(m_iElement ^ element.m_iElement);
		m_iDegree = m_gf.getDegreeByElement(m_iElement);
		return *this;
	}
	
	const GaloisFieldElement& operator+=(unsigned int iElement)
	{
		m_iElement = m_gf.moduloCast(m_iElement ^ iElement);
		m_iDegree = m_gf.getDegreeByElement(m_iElement);
		return *this;
	}

	GaloisFieldElement operator*(const GaloisFieldElement& element) const
	{
		GaloisFieldElement result;
		if(m_iElement && element.m_iElement)
		{
			result = m_gf.getElementByDegree(m_iDegree + element.m_iDegree);
		}
		return result;
	}

	GaloisFieldElement operator*(unsigned int iElement) const
	{
		GaloisFieldElement result;
		if(m_iElement && iElement)
		{
			result = m_gf.getElementByDegree(m_iDegree + m_gf.getDegreeByElement(iElement));
		}
		return result;
	}
	
	GaloisFieldElement operator/(const GaloisFieldElement& element) const
	{
		GaloisFieldElement result;
		if (!element.m_iElement)
		{
			printf("Divide by zero\n");
			exit(1);
		}
		if(m_iElement)
		{
			result = m_gf.getElementByDegree(m_iDegree - element.m_iDegree);
		}
		return result;
	}

	GaloisFieldElement operator^(unsigned int iDegree) const
	{
		GaloisFieldElement result(m_gf.getElementByDegree(m_iDegree * iDegree));
		return result;
	}

	void setErased(bool bState=true)
	{
		m_bErased = bState;
	}

	bool getErased() const
	{
		return m_bErased;
	}

private:
	const T& m_gf;
	int		 m_iElement;
	int		 m_iDegree;
	bool	 m_bErased;	
};

#endif
