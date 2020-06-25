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
#include <stdio.h>
#include <vector>
#include <mex.h>
#include "common.h"
#include "func_nodes.h"
#include "ObjectHandle.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 6)
    {
        mexErrMsgTxt("usage: initialize_func_nodes(M, N, row_weights, S_positions, S_values, mod_table)");
    }

    int T = (int)mxGetScalar(prhs[0]);
    int N = (int)mxGetScalar(prhs[1]);

    double* input = (double*)mxGetData(prhs[2]);
    int* rowWeights = new int[T];
    int dmax = 0;
    for (int i = 0; i < T; ++i)
    {
        rowWeights[i] = (int)input[i];
        if (rowWeights[i] > dmax)
        {
            dmax = rowWeights[i];
        }
    }

    Matrix<int> S_positions(T, dmax);
    input = (double*)mxGetData(prhs[3]);
    for (int i = 0; i < T; ++i)
    {
        for (int j = 0; j < dmax; ++j)
        {
            S_positions(i, j) = (int)input[T*j + i];
        }
    }

    Matrix<complex<double> > S_values(T, dmax);
    double* inputR = (double*)mxGetPr(prhs[4]);
    double* inputI = (double*)mxGetPi(prhs[4]);

    for (int i = 0; i < T; ++i)
    {
        for (int j = 0; j < dmax; ++j)
        {
            if (inputI)
            {
                S_values(i, j) = complex<double>(inputR[T*j + i], inputI[T*j + i]);
            }
            else
            {
                S_values(i, j) = complex<double>(inputR[T*j + i], 0);
            }
        }
    }

    int M = mxGetN(prhs[5]);
    complex<double>* modTable = new complex<double>[M];
    inputR = (double*)mxGetPr(prhs[5]);
    inputI = (double*)mxGetPi(prhs[5]);

    for (int i = 0; i < M; ++i)
    {
        if (inputI)
        {
            modTable[i] = complex<double>(inputR[i], inputI[i]);
        }
        else
        {
            modTable[i] = complex<double>(inputR[i], 0);
        }
    }

    FuncConnections* pFC = new FuncConnections(T, N, rowWeights, S_positions, S_values, M, modTable);
    plhs[0] = create_handle(pFC);

    delete[] rowWeights;
}
