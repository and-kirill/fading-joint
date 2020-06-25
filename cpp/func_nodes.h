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

#ifndef FUNC_CONNECTIONS_H
#define FUNC_CONNECTIONS_H

#include "Matrix.h"
#include <complex>
using namespace std;

class FuncConnections
{
public:
    FuncConnections(const int T, const int N, int* rowWeights, Matrix<int>& S_positions, Matrix<complex<double> >& S_values, int M, complex<double>* modTable)
		: m_T(T)
		, m_N(N)
        , m_rowWeights(0)
		, m_S_positions(S_positions)
		, m_S_values(S_values)
        , m_M(M)
        , m_modTable(0)
	{
        m_rowWeights = new int[T];
        m_dmax = 0;
        for (int i = 0; i < T; ++i)
        {
            m_rowWeights[i] = rowWeights[i];
            if (rowWeights[i] > m_dmax)
            {
                m_dmax = rowWeights[i];
            }
        }

        m_modTable = new complex<double>[M];
        for (int i = 0; i < M; ++i)
        {
            m_modTable[i] = modTable[i];
        }
	}

	~FuncConnections()
	{
	    delete[] m_rowWeights;
	    delete[] m_modTable;
	}
	
public:
	int					            m_T;
	int					            m_N;
	int                             m_dmax;
	int*                            m_rowWeights;
	Matrix<int>	                    m_S_positions;
	Matrix<complex<double> >	    m_S_values;
	int                             m_M;
	complex<double>*                m_modTable;
};


#endif
