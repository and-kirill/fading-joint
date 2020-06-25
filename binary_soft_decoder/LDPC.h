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

#ifndef LDPC_H
#define LDPC_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <vector>

#include "../common/Matrix.h"
#include "../common/common.h"

using namespace std;

class LDPC
{
public:
    int n;
    int m;
    int cmax;
    int rmax;
    vector<int> col_weight;
    vector<int> row_weight;
    Matrix<int> row_col; // positions of ones in rows
    Matrix<int> row_N; // positions of ones in col_row
    Matrix<int> col_row; // positions of ones in columns
    Matrix<int> col_N; // positions of ones in row_col
    Matrix<FieldElement> H; // Sparse parity-check matrix
    Matrix<std::pair<int, int> > msgs_col; // pair shows the position of msg in R_msg and Q_msg

public:
    LDPC()
    {}
    ~LDPC()
    {}
    void init(char *s)
    {
        int q, *count;
        FILE *fp = fopen(s, "rt");
        if (fp == NULL)
        {
            mexPrintf("ERROR: Cannot open file '%s'\n", s);
            return;
        }
        char buf[256];
        fgets(buf, sizeof(buf), fp);
        if (sscanf(buf, "%d%d%d", &n, &m, &q) == 2)
        {
            q = 2;
        }
        if (q != Q)
        {
            mexPrintf("ERROR: GF order(%d) does not match with matrix file (%d)\n", Q, q);
            return;
        }

        fscanf(fp, "%d%d", &cmax, &rmax);
        col_weight.resize(n);
        for (int i = 0; i < n; i++)
        {
            fscanf(fp, "%d", &col_weight[i]);
        }
        row_weight.resize(m);
        for (int j = 0; j < m; j++)
        {
            fscanf(fp, "%d", &row_weight[j]);
        }
        
        row_col = Matrix<int>(m, rmax);
        row_N = Matrix<int>(m, rmax);
        col_row = Matrix<int>(n, cmax);
        col_N = Matrix<int>(n, cmax);
        H = Matrix<FieldElement>(m, rmax);
        msgs_col = Matrix<std::pair<int, int> >(n, cmax);
        
        count = (int*)malloc(sizeof(int) * n);
        memset(count, 0, sizeof(int) * n);

        //skip n lines
        if (q > 2)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < cmax; j++)
                {
                    fscanf(fp, "%*d%*d");
                }
            }
        }
        else
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < cmax; j++)
                {
                    fscanf(fp, "%*d");
                }
            }
        }

        if (q > 2)
        {
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < row_weight[i]; ++j)
                {
                    int v, u;
                    fscanf(fp, "%d%d", &v, &u);
                    H(i,j) = FieldElement(u);
                    v--;
                    row_col(i,j) = v;
                    row_N(i,j) = count[v];
                    col_row(v,count[v]) = i;
                    col_N(v,count[v]) = j;
                    count[v]++;
                }

                //skip fillers
                int a,b;
                for (int j = row_weight[i]; j < rmax; j++)
                {
                  fscanf(fp, "%d%d", &a, &b);
                  if (a!=0 || b!=0)
                  {
                      mexPrintf("ERROR: fillers at row %d, %d %d\n", j, a, b);
                      return;
                  }
                }
            }
        }
        else
        {
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < row_weight[i]; ++j)
                {
                    int v;
                    fscanf(fp, "%d", &v);
                    H(i,j) = 1;
                    v--;
                    row_col(i,j) = v;
                    row_N(i,j) = count[v];
                    col_row(v,count[v]) = i;
                    col_N(v,count[v]) = j;
                    count[v]++;
                }

                //skip fillers
                int a;
                for (int j = row_weight[i]; j < rmax; j++)
                {
                  fscanf(fp, "%d", &a);
                  if (a!=0)
                  {
                      mexPrintf("ERROR: fillers at row %d, %d\n", j, a);
                      return;
                  }
                }
            }
        }

        free(count);
        fclose(fp);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < col_weight[i]; j++)
            {
                msgs_col(i,j) = std::pair<int, int>(col_row(i,j), col_N(i,j));
            }
        }
    }
    bool syndrome(const Matrix<FieldElement>& r, Matrix<FieldElement>& s)
    {
        s = Matrix<FieldElement>(1, m);
        bool bCodeword = true;
        
        for (int i = 0; i < m; ++i)
        {
            s(i) = 0;
            for (int j = 0; j < row_weight[j]; ++j)
            {
                s(i) += H(i, j)*r(row_col(i,j));
            }
            if (s(i).getElement())
            {
                bCodeword = false;
            }
        }

        return bCodeword;
    }
};

#endif