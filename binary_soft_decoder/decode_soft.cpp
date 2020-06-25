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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <vector>
#include <list>

#include <mex.h>
#include <matrix.h>

#include "Matrix.h"
#include "common.h"
#include "LDPC.h"
#include "ObjectHandle.h"

#define ABS(A) (((A) >= 0) ? (A) : -(A))
#define HARD(A) (((A) < 0)?1:0)
#define OFFSET(A, B) (((A) > (B)) ? ((A)-(B)) : 0)
#define INFTY 1000000

using namespace std;

double safe_exp(double x)
{
    if (x > 700)
    {
        return 1.014232054735005e+304;
    }
    if (x < -700)
    {
        return 9.859676543759770e-305;
    }
    return exp(x);
}

double safe_log(double x)
{
    if (x < 9.859676543759770e-305)
    {
        return -700;
    }
    return log(x);
}


#if 0
double phi(double x)
{
    static const double lim_tanh = 31.0;
    static const double min_tanh = 6.883382752676208e-14; //log( (exp((double)lim) + 1)/(exp((double)lim) - 1));
    
    if (x > lim_tanh)
    {
      return 2*safe_exp(-x);
    }
    
    //if (x < min_tanh)
    //{
    //  return -safe_log(x/2);
    //}
    
    return -safe_log(tanh(x/2));    
}
#else
static float phi(float x)
{
	//	const union {
	//		long long i;
	//		double d;
	//	} lim = { 0x3d335fffffffff44 };
	static const float lim = 31.;
	static const float inv_lim = 
		log( (exp((double)lim) + 1)/(exp((double)lim) - 1) );
	
	if (x > lim)
		return( 0 );
	else if (x < inv_lim)
		return( lim );
	else {
		double t = exp( (double) x );
		return (float) log( (t+1)/(t-1) ); 
	}
}

static float phi_sign(float x)
{
	if (x < 0)
	{
		return -phi(fabs(x));
	}
	return phi(x);
}
#endif

//=======================================================================================================================================
//============================================================= Decoders (Log Domain) ===================================================
//=======================================================================================================================================


// Sum-Product Decoder
bool SumProduct(LDPC& ldpc, vector<double>& in_llr, int max_iter, vector<double>& out_llr, vector<int>& y, int* number_of_iter)
{
    // Auxiliary matrices
    Matrix<double>  R_msgs(ldpc.m, ldpc.rmax); // messages from check to variable nodes
    Matrix<double>  Q_tanhs(ldpc.m, ldpc.rmax);
    Matrix<int>  Q_signs(ldpc.m, ldpc.rmax);
    double sum_tanhs = 0;
    int sum_sign = 0;
    int sign = 0;
    double temp = 0;


    // Initialization
    for (int i = 0; i < ldpc.n; ++i)
    {
        out_llr[i] = in_llr[i];
        y[i] = HARD(out_llr[i]);
	  
        for (int j = 0; j < ldpc.col_weight[i]; ++j)
        {
            Q_signs(ldpc.msgs_col(i,j)) = 0;
            if (in_llr[i] < 0)
            {
                Q_signs(ldpc.msgs_col(i,j)) = 1;
            }
            Q_tanhs(ldpc.msgs_col(i,j)) = phi(fabs(in_llr[i]));
        }
    }

    for (int loop = 0; loop < max_iter; ++loop)
    {
        // Update R
        for (int j = 0; j < ldpc.m; j++)
        {
            sum_tanhs = 0;
            sum_sign = 0;

            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                sum_tanhs += Q_tanhs(j, k);
                sum_sign ^= Q_signs(j, k);
            }
            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                sign = sum_sign^Q_signs(j, k);
                R_msgs(j,k) = (1-2*sign)*phi(sum_tanhs - Q_tanhs(j, k));
            }
        }
        // Update Q
        for (int i = 0; i < ldpc.n; i++)
        {
            out_llr[i] = in_llr[i];
        
            for (int k = 0; k < ldpc.col_weight[i]; k++)
            {
                out_llr[i] += R_msgs(ldpc.msgs_col(i,k));
            }

            y[i] = HARD(out_llr[i]);
            
	        for (int k = 0; k < ldpc.col_weight[i]; k++)
		    {
	            temp = out_llr[i] - R_msgs(ldpc.msgs_col(i,k));
		        Q_signs(ldpc.msgs_col(i,k)) = 0;
                if (temp < 0)
                {
                    Q_signs(ldpc.msgs_col(i,k)) = 1;
                }
                Q_tanhs(ldpc.msgs_col(i,k)) = phi(fabs(temp));
		    }
        }
    }

    if (number_of_iter)
    {
        *number_of_iter = max_iter;
    }
    return 1;
}

// Layered Sum-Product Decoder
bool LayeredSumProduct(LDPC& ldpc, vector<double>& in_llr, int max_iter, vector<double>& out_llr, vector<int>& y, int* number_of_iter, vector<int>& row_sequence)
{
    Matrix<double>  R_msgs(ldpc.m, ldpc.rmax);
    double* Q_tanhs = new double[ldpc.rmax];
    int* Q_signs = new int[ldpc.rmax];
    double sum_tanhs = 0;
    int sign = 0;
    int sum_sign = 0;
    double Q_msg = 0;

    for (int i = 0; i < ldpc.n; ++i)
    {
        out_llr[i] = in_llr[i];
        y[i] = HARD(out_llr[i]);
    }

    for (int i = 0; i < ldpc.m; ++i)
    {
        for (int j = 0; j < ldpc.rmax; ++j)
        {
            R_msgs(i,j) = 0; // very important
        }
    }

    for (int loop = 0; loop < max_iter; ++loop)
    {
        for (int t = 0; t < ldpc.m; t++)
        {
            int j = row_sequence[t];
            sum_tanhs = 0;
            sum_sign = 0;

            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                Q_msg = out_llr[ldpc.row_col(j, k)] - R_msgs(j,k);
                Q_signs[k] = 0;
                if (Q_msg < 0)
                {
                    Q_signs[k] = 1;
                    sum_sign ^= 1;
                }
                Q_tanhs[k] = phi(fabs(Q_msg));
                sum_tanhs += Q_tanhs[k];

            }
            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                sign = sum_sign^Q_signs[k];
                out_llr[ldpc.row_col(j, k)] -= R_msgs(j,k);
                R_msgs(j,k) = (1-2*sign)*phi(sum_tanhs - Q_tanhs[k]);
                out_llr[ldpc.row_col(j, k)] += R_msgs(j,k);
            }
        }
    }

    delete[] Q_tanhs;
    delete[] Q_signs;

    for (int i = 0; i < ldpc.n; i++)
    {
        y[i] = HARD(out_llr[i]);
    }

    if (number_of_iter)
    {
        *number_of_iter = max_iter;
    }
    return 1;
}

// Min-Sum Decoder
bool MinSum(LDPC& ldpc, vector<double>& in_llr, int max_iter, vector<double>& out_llr, vector<int>& y, int* number_of_iter, vector<double>& scale_array, vector<double>& offset_array)
{
    // Auxiliary matrices
    Matrix<double>  R_msgs(ldpc.m, ldpc.rmax); // messages from check to variable nodes
    Matrix<double>  Q_msgs(ldpc.m, ldpc.rmax); // messages from variable to check nodes

    // Initialization
    for (int i = 0; i < ldpc.n; ++i)
    {
        out_llr[i] = in_llr[i];
        y[i] = HARD(out_llr[i]);

        for (int j = 0; j < ldpc.col_weight[i]; ++j)
        {
            Q_msgs(ldpc.msgs_col(i,j)) = in_llr[i];
        }
    }

    for (int loop = 0; loop < max_iter; ++loop)
    {
        // Update R
        for (int j = 0; j < ldpc.m; j++)
        {
            double first_min = INFTY;
            double second_min = INFTY;
            int min_index = 0;
            int sum_sign = 0;

            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                double temp = fabs(Q_msgs(j, k));
                if (temp < first_min)
                {
                    second_min = first_min;
                    first_min = temp;
                    min_index = k;
                }
                else if (temp < second_min)
                {
                    second_min = temp;
                }

                if (Q_msgs(j, k) < 0)
                {
                    sum_sign ^= 1;
                }
            }
            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                int sign = sum_sign;
                if (Q_msgs(j, k) < 0)
                {
                    sign ^= 1;
                }
                if (k == min_index)
                {
                    R_msgs(j,k) = second_min;
                }
                else
                {
                    R_msgs(j,k) = first_min;
                }
                double scale = scale_array[ldpc.row_col(j, k)];
                double offset = offset_array[ldpc.row_col(j, k)];
                R_msgs(j,k) = (1-2*sign)*scale*OFFSET(R_msgs(j,k), offset);
            }
        }
        // Update Q
        for (int i = 0; i < ldpc.n; i++)
        {
            out_llr[i] = in_llr[i];

            for (int k = 0; k < ldpc.col_weight[i]; k++)
            {
                out_llr[i] += R_msgs(ldpc.msgs_col(i,k));
            }

            y[i] = HARD(out_llr[i]);

            for (int k = 0; k < ldpc.col_weight[i]; k++)
            {
                Q_msgs(ldpc.msgs_col(i,k)) = out_llr[i] - R_msgs(ldpc.msgs_col(i,k));
            }
        }
    }

    if (number_of_iter)
    {
        *number_of_iter = max_iter;
    }
    return 1;
}

// Layered Min-Sum Decoder
bool LayeredMinSum(LDPC& ldpc, vector<double>& in_llr, int max_iter, vector<double>& out_llr, vector<int>& y, int* number_of_iter, vector<int>& row_sequence, vector<double>& scale_array, vector<double>& offset_array)
{
    Matrix<double>  R_msgs(ldpc.m, ldpc.rmax);

    for (int i = 0; i < ldpc.n; ++i)
    {
        out_llr[i] = in_llr[i];
        y[i] = HARD(out_llr[i]);
    }

    for (int i = 0; i < ldpc.m; ++i)
    {
        for (int j = 0; j < ldpc.rmax; ++j)
        {
            R_msgs(i,j) = 0; // very important
        }
    }

    for (int loop = 0; loop < max_iter; ++loop)
    {
        for (int t = 0; t < ldpc.m; t++)
        {
            int j = row_sequence[t];
            double first_min = INFTY;
            double second_min = INFTY;
            int min_index = 0;
            int sum_sign = 0;

            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                double Q_msg = out_llr[ldpc.row_col(j, k)] - R_msgs(j,k);
                double temp = fabs(Q_msg);
                if (temp < first_min)
                {
                    second_min = first_min;
                    first_min = temp;
                    min_index = k;
                }
                else if (temp < second_min)
                {
                    second_min = temp;
                }

                if (Q_msg < 0)
                {
                    sum_sign ^= 1;
                }
            }
            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                int sign = sum_sign;
                double Q_msg = out_llr[ldpc.row_col(j, k)] - R_msgs(j,k);
                if (Q_msg < 0)
                {
                    sign ^= 1;
                }
                out_llr[ldpc.row_col(j, k)] = Q_msg;
                if (k == min_index)
                {
                    R_msgs(j,k) = second_min;
                }
                else
                {
                    R_msgs(j,k) = first_min;
                }
                double scale = scale_array[ldpc.row_col(j, k)];
                double offset = offset_array[ldpc.row_col(j, k)];
                R_msgs(j,k) = (1-2*sign)*scale*OFFSET(R_msgs(j,k), offset);
                out_llr[ldpc.row_col(j, k)] += R_msgs(j,k);
            }
        }
    }

    for (int i = 0; i < ldpc.n; i++)
    {
        y[i] = HARD(out_llr[i]);
    }

    if (number_of_iter)
    {
        *number_of_iter = max_iter;
    }
    return 1;
}

// Adjusted Min-Sum Decoder
bool AdjMinSum(LDPC& ldpc, vector<double>& in_llr, int max_iter, vector<double>& out_llr, vector<int>& y, int* number_of_iter, bool bEnhanced = 0)
{
    // Auxiliary matrices
    Matrix<double>  R_msgs(ldpc.m, ldpc.rmax); // messages from check to variable nodes
    Matrix<double>  Q_abs(ldpc.m, ldpc.rmax); // messages from variable to check nodes
    Matrix<int>     Q_signs(ldpc.m, ldpc.rmax);
    double sum_tanhs = 0;
    int sum_sign = 0;
    int sign = 0;
    double temp = 0;
    double Q_min = INFTY;
    double Q_max = 0;
    int min_index = 0;
    int max_index = 0;    

    // Initialization
    for (int i = 0; i < ldpc.n; ++i)
    {
        out_llr[i] = in_llr[i];
        y[i] = HARD(out_llr[i]);
	  
        for (int j = 0; j < ldpc.col_weight[i]; ++j)
        {
            Q_signs(ldpc.msgs_col(i,j)) = 0;
            if (in_llr[i] < 0)
            {
                Q_signs(ldpc.msgs_col(i,j)) = 1;
            }
            Q_abs(ldpc.msgs_col(i,j)) = fabs(in_llr[i]);
        }
    }

    for (int loop = 0; loop < max_iter; ++loop)
    {
        // Update R
        for (int j = 0; j < ldpc.m; j++)
        {
            sum_tanhs = 0;
            sum_sign = 0;
            Q_min = INFTY;
            Q_max = 0;
            min_index = 0;
            max_index = 0;

            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                sum_tanhs += phi(Q_abs(j, k));
                sum_sign ^= Q_signs(j, k);
                if (Q_abs(j, k) < Q_min)
                {
                    Q_min = Q_abs(j, k);
                    min_index = k;
                }
                if (Q_abs(j, k) >= Q_max)
                {
                    Q_max = Q_abs(j, k);
                    max_index = k;
                }
            }
            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                sign = sum_sign^Q_signs(j, k);
                if (!bEnhanced)
                {
                    /* Ordinar Adjusted Min-Sum. Along the edge carrying
                     * the minimum incoming magnitude LLR the outgoing 
                     * message is the one computed by the SP rule. And 
                     * along the rest of the edges (second stored message) 
                     * the outgoing message magnitude used is the one 
                     * obtained by applying the SP rule to all the incoming 
                     * messages
                     */ 
                    if (k == min_index)
                    {
                       R_msgs(j,k) = (1-2*sign)*phi(sum_tanhs - phi(Q_min)); 
                    }
                    else
                    {
                       R_msgs(j,k) = (1-2*sign)*phi(sum_tanhs); 
                    }
                }
                else
                {
                    /* Enhanced Adjusted Min-Sum. Better performance is 
                     * obtained when for the second stored message, the SP 
                     * outgoing message along the maximum incoming 
                     * reliability (or in other words the outgoing SP 
                     * message with minimum reliability) is used instead. 
                     * This decoder is a slight modification of [5] and is 
                     * observed to yield better performance than the one 
                     * provided in [5]. 
                     */
                    if (k == min_index)
                    {
                       R_msgs(j,k) = (1-2*sign)*phi(sum_tanhs - phi(Q_min)); 
                    }
                    else
                    {
                       R_msgs(j,k) = (1-2*sign)*phi(sum_tanhs - phi(Q_max)); 
                    }
                }
            }
        }
        // Update Q
        for (int i = 0; i < ldpc.n; i++)
        {
            out_llr[i] = in_llr[i];
        
            for (int k = 0; k < ldpc.col_weight[i]; k++)
            {
                out_llr[i] += R_msgs(ldpc.msgs_col(i,k));
            }

            y[i] = HARD(out_llr[i]);
            
	        for (int k = 0; k < ldpc.col_weight[i]; k++)
		    {
	            temp = out_llr[i] - R_msgs(ldpc.msgs_col(i,k));
		        Q_signs(ldpc.msgs_col(i,k)) = 0;
                if (temp < 0)
                {
                    Q_signs(ldpc.msgs_col(i,k)) = 1;
                }
                Q_abs(ldpc.msgs_col(i,k)) = fabs(temp);
		    }
        }
    }

    if (number_of_iter)
    {
        *number_of_iter = max_iter;
    }
    return 1;
}

// Layered Adjusted Min-Sum Decoder
bool LayeredAdjMinSum(LDPC& ldpc, vector<double>& in_llr, int max_iter, vector<double>& out_llr, vector<int>& y, int* number_of_iter, vector<int>& row_sequence, bool bEnhanced = 0)
{
    Matrix<double>  R_msgs(ldpc.m, ldpc.rmax);
    int* Q_signs = new int[ldpc.rmax];
    double sum_tanhs = 0;
    int sign = 0;
    int sum_sign = 0;
    double Q_msg = 0;
    double Q_min = INFTY;
    double Q_max = 0;
    int min_index = 0;
    int max_index = 0; 

    for (int i = 0; i < ldpc.n; ++i)
    {
        out_llr[i] = in_llr[i];
        y[i] = HARD(out_llr[i]);
    }

    for (int i = 0; i < ldpc.m; ++i)
    {
        for (int j = 0; j < ldpc.rmax; ++j)
        {
            R_msgs(i,j) = 0; // very important
        }
    }

    for (int loop = 0; loop < max_iter; ++loop)
    {
        for (int t = 0; t < ldpc.m; t++)
        {
            int j = row_sequence[t];
            sum_tanhs = 0;
            sum_sign = 0;
            Q_min = INFTY;
            Q_max = 0;
            min_index = 0;
            max_index = 0;

            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                Q_msg = out_llr[ldpc.row_col(j, k)] - R_msgs(j,k);
                Q_signs[k] = 0;
                if (Q_msg < 0)
                {
                    Q_signs[k] = 1;
                    sum_sign ^= 1;
                }
                Q_msg = fabs(Q_msg);
                sum_tanhs += phi(Q_msg);
                
                if (Q_msg < Q_min)
                {
                    Q_min = Q_msg;
                    min_index = k;
                }
                if (Q_msg >= Q_max)
                {
                    Q_max = Q_msg;
                    max_index = k;
                }
            }
            
            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                sign = sum_sign^Q_signs[k];
                out_llr[ldpc.row_col(j, k)] -= R_msgs(j,k);
                
                if (!bEnhanced)
                {
                    /* Ordinar Adjusted Min-Sum. Along the edge carrying
                     * the minimum incoming magnitude LLR the outgoing 
                     * message is the one computed by the SP rule. And 
                     * along the rest of the edges (second stored message) 
                     * the outgoing message magnitude used is the one 
                     * obtained by applying the SP rule to all the incoming 
                     * messages
                     */ 
                    if (k == min_index)
                    {
                       R_msgs(j,k) = (1-2*sign)*phi(sum_tanhs - phi(Q_min)); 
                    }
                    else
                    {
                       R_msgs(j,k) = (1-2*sign)*phi(sum_tanhs); 
                    }
                }
                else
                {
                    /* Enhanced Adjusted Min-Sum. Better performance is 
                     * obtained when for the second stored message, the SP 
                     * outgoing message along the maximum incoming 
                     * reliability (or in other words the outgoing SP 
                     * message with minimum reliability) is used instead. 
                     * This decoder is a slight modification of [5] and is 
                     * observed to yield better performance than the one 
                     * provided in [5]. 
                     */
                    if (k == min_index)
                    {
                       R_msgs(j,k) = (1-2*sign)*phi(sum_tanhs - phi(Q_min)); 
                    }
                    else
                    {
                       R_msgs(j,k) = (1-2*sign)*phi(sum_tanhs - phi(Q_max)); 
                    }
                }
                
                out_llr[ldpc.row_col(j, k)] += R_msgs(j,k);
            }
        }
    }

    delete[] Q_signs;

    for (int i = 0; i < ldpc.n; i++)
    {
        y[i] = HARD(out_llr[i]);
    }

    if (number_of_iter)
    {
        *number_of_iter = max_iter;
    }
    return 1;
}


// Parallel Bit-Flipping
bool BitFlipping(LDPC& ldpc, vector<int>& x, int max_iter, vector<int>& y, int* number_of_iter)
{
    Matrix<int>  R_msgs(ldpc.m, ldpc.rmax);
    Matrix<int>  Q_msgs(ldpc.m, ldpc.rmax);
    int badCodeword = 0;
    
    for (int i = 0; i < ldpc.n; ++i)
    {
        y[i] = x[i];
        
        for (int j = 0; j < ldpc.col_weight[i]; ++j)
        {
            Q_msgs(ldpc.msgs_col(i,j)) = y[i];
        }
    }

    for (int loop = 0; loop < max_iter; ++loop)
    {
        // Update R
        badCodeword = 0;
        for (int j = 0; j < ldpc.m; j++)
        {
            int sum = 0;
            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                sum ^= Q_msgs(j, k);
            }
            if (sum)
            {
            	badCodeword = 1;
            }
            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                R_msgs(j,k) = sum; // ???? ^ Q_msgs(j, k)
            }
        }
        if (!badCodeword) // zero syndrome -> decoding is successful
        {
        	if (number_of_iter)
		    {
		        *number_of_iter = max_iter;
		    }
		    return 1;
        }
        // Update Q
        for (int i = 0; i < ldpc.n; i++)
        {
            int unsatisfied = 0;
            for (int k = 0; k < ldpc.col_weight[i]; k++)
            {
                if (R_msgs(ldpc.msgs_col(i,k)))
                {
                	unsatisfied += 1;
                }
            }
            if (unsatisfied > ldpc.col_weight[i]/2)
            {
            	y[i] ^= 1; // flip the bit
            }
            for (int k = 0; k < ldpc.col_weight[i]; k++)
		    {
		        Q_msgs(ldpc.msgs_col(i,k)) = y[i];
	        }
        }
    }

    if (number_of_iter)
    {
        *number_of_iter = max_iter;
    }
    return 0;
}


// Sum-Product FFT Decoder
bool SumProduct_FFT(LDPC& ldpc, vector<double>& in_llr_fft, int max_iter, vector<double>& out_llr, vector<int>& y, int* number_of_iter)
{
    // Auxiliary matrices
    Matrix<double>  R_msgs(ldpc.m, ldpc.rmax); // messages from check to variable nodes
    Matrix<int>  R_signs(ldpc.m, ldpc.rmax);
    Matrix<double>  Q_msgs(ldpc.m, ldpc.rmax);
    Matrix<int>  Q_signs(ldpc.m, ldpc.rmax);
    double sum = 0;
    int sum_sign = 0;
    double temp = 0;

    // Initialization
    for (int i = 0; i < ldpc.n; ++i)
    {
        out_llr[i] = phi_sign(in_llr_fft[i]);
        y[i] = HARD(out_llr[i]);
	  
        for (int j = 0; j < ldpc.col_weight[i]; ++j)
        {
            Q_signs(ldpc.msgs_col(i,j)) = 0;
            if (in_llr_fft[i] < 0)
            {
                Q_signs(ldpc.msgs_col(i,j)) = 1;
            }
            Q_msgs(ldpc.msgs_col(i,j)) = fabs(in_llr_fft[i]);
        }
    }

    for (int loop = 0; loop < max_iter; ++loop)
    {
        // Update R
        for (int j = 0; j < ldpc.m; j++)
        {
            sum = 0;
            sum_sign = 0;

            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                sum += Q_msgs(j, k);
                sum_sign ^= Q_signs(j, k);
            }
            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                R_signs(j,k) = sum_sign^Q_signs(j, k);
                R_msgs(j,k) = sum - Q_msgs(j, k);
            }
        }
        // Update Q
        for (int i = 0; i < ldpc.n; i++)
        {
            out_llr[i] = phi_sign(in_llr_fft[i]);
        
            for (int k = 0; k < ldpc.col_weight[i]; k++)
            {
                out_llr[i] += (1-2*R_signs(ldpc.msgs_col(i,k)))*phi(R_msgs(ldpc.msgs_col(i,k)));
            }

            y[i] = HARD(out_llr[i]);
            
	        for (int k = 0; k < ldpc.col_weight[i]; k++)
		    {
	            temp = out_llr[i] - (1-2*R_signs(ldpc.msgs_col(i,k)))*phi(R_msgs(ldpc.msgs_col(i,k)));
		        Q_signs(ldpc.msgs_col(i,k)) = 0;
                if (temp < 0)
                {
                    Q_signs(ldpc.msgs_col(i,k)) = 1;
                }
                Q_msgs(ldpc.msgs_col(i,k)) = phi(fabs(temp));
		    }
        }
    }

    if (number_of_iter)
    {
        *number_of_iter = max_iter;
    }
    return 1;
}

// Min-Sum FFT Decoder
bool MinSum_FFT(LDPC& ldpc, vector<double>& in_llr_fft, int max_iter, vector<double>& out_llr, vector<int>& y, int* number_of_iter, vector<double>& scale_array, vector<double>& offset_array)
{
    // Auxiliary matrices
    Matrix<double>  R_msgs(ldpc.m, ldpc.rmax); // messages from check to variable nodes
    Matrix<int>  R_signs(ldpc.m, ldpc.rmax);
    Matrix<double>  Q_msgs(ldpc.m, ldpc.rmax); // messages from variable to check nodes
    Matrix<int>  Q_signs(ldpc.m, ldpc.rmax);
    double sum = 0;
    int sum_sign = 0;

    // Initialization
    for (int i = 0; i < ldpc.n; ++i)
    {
        out_llr[i] = phi_sign(in_llr_fft[i]);
        y[i] = HARD(out_llr[i]);
	  
        for (int j = 0; j < ldpc.col_weight[i]; ++j)
        {
            Q_signs(ldpc.msgs_col(i,j)) = 0;
            if (in_llr_fft[i] < 0)
            {
                Q_signs(ldpc.msgs_col(i,j)) = 1;
            }
            Q_msgs(ldpc.msgs_col(i,j)) = fabs(in_llr_fft[i]);
        }
    }

    for (int loop = 0; loop < max_iter; ++loop)
    {
        // Update R
        for (int j = 0; j < ldpc.m; j++)
        {
            sum = 0;
            sum_sign = 0;

            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                sum += Q_msgs(j, k);
                sum_sign ^= Q_signs(j, k);
            }
            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                R_signs(j,k) = sum_sign^Q_signs(j, k);
                R_msgs(j,k) = sum - Q_msgs(j, k);
            }
        }
        // Update Q
#if 1
        for (int i = 0; i < ldpc.n; i++)
		{
			int in_sign = (in_llr_fft[i]<0);
			double first_min = INFTY;
			double second_min = INFTY;
            int min_index = 0;
            int sign1 = 0;
			int sign2 = 0;

            for (int k = 0; k < ldpc.col_weight[i]; k++)
            {
            	double temp = R_msgs(ldpc.msgs_col(i,k));
                if (temp < first_min)
                {
                    second_min = first_min;
                    sign2 = sign1;
                    first_min = temp;
                    sign1 = R_signs(ldpc.msgs_col(i,k));
                    min_index = k;
                }
                else if (temp < second_min)
                {
                    second_min = temp;
                    sign2 = R_signs(ldpc.msgs_col(i,k));
                }
            }

            for (int k = 0; k < ldpc.col_weight[i]; k++)
            {
            	if (fabs(in_llr_fft[i]) < first_min)
            	{
            		out_llr[i] = (1-2*in_sign)*phi(fabs(in_llr_fft[i]));
            	}
            	else
            	{
            		out_llr[i] = (1-2*sign1)*phi(first_min);		
            	}
                
            }

            y[i] = HARD(out_llr[i]);

            for (int k = 0; k < ldpc.col_weight[i]; k++)
            {
            	double scale = scale_array[i];
                double offset = offset_array[i];
            	if (k == min_index)
                {
                	if (fabs(in_llr_fft[i]) >= first_min)
            		{
                		Q_signs(ldpc.msgs_col(i,k)) = sign2;
                    	Q_msgs(ldpc.msgs_col(i,k)) = scale*OFFSET(second_min, offset);
                    }
                    else
                    {
                    	Q_signs(ldpc.msgs_col(i,k)) = in_sign;
                    	Q_msgs(ldpc.msgs_col(i,k)) = scale*OFFSET(fabs(in_llr_fft[i]), offset);
                    }
                }
                else
                {
                	if (fabs(in_llr_fft[i]) >= first_min)
            		{
                		Q_signs(ldpc.msgs_col(i,k)) = sign1;
                    	Q_msgs(ldpc.msgs_col(i,k)) = scale*OFFSET(first_min, offset);
                    }
                    else
                    {
                    	Q_signs(ldpc.msgs_col(i,k)) = in_sign;
                    	Q_msgs(ldpc.msgs_col(i,k)) = scale*OFFSET(fabs(in_llr_fft[i]), offset);
                    }
                }
            }
        }
#else
        double temp = 0;
        for (int i = 0; i < ldpc.n; i++)
        {
            out_llr[i] = phi_sign(in_llr_fft[i]);
        
            for (int k = 0; k < ldpc.col_weight[i]; k++)
            {
                out_llr[i] += (1-2*R_signs(ldpc.msgs_col(i,k)))*phi(R_msgs(ldpc.msgs_col(i,k)));
            }

            y[i] = HARD(out_llr[i]);
            
	        for (int k = 0; k < ldpc.col_weight[i]; k++)
		    {
	            temp = out_llr[i] - (1-2*R_signs(ldpc.msgs_col(i,k)))*phi(R_msgs(ldpc.msgs_col(i,k)));
		        Q_signs(ldpc.msgs_col(i,k)) = 0;
                if (temp < 0)
                {
                    Q_signs(ldpc.msgs_col(i,k)) = 1;
                }
                Q_msgs(ldpc.msgs_col(i,k)) = phi(fabs(temp));
		    }
        }
#endif        
    }

    if (number_of_iter)
    {
        *number_of_iter = max_iter;
    }
    return 1;
}


//=======================================================================================================================================
//================================================================== Mex function =======================================================
//=======================================================================================================================================


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int command = (int)*(mxGetPr(prhs[0]));
    double* input = 0;

    switch(command)
    {
        case 0:
        {
            char buf[256];
            mxGetString(prhs[1], buf, 256); 
             
            LDPC* pLdpc = new LDPC();
            pLdpc->init(buf);
            
            if (nlhs > 0)
            {
	            plhs[0] = create_handle(pLdpc);
            }
            if (nlhs > 1)
            {
                plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
                *((double*)(mxGetPr(plhs[1]))) = pLdpc->n;
            }
            if (nlhs > 2)
            {	
                plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
                *((double*)(mxGetPr(plhs[2]))) = pLdpc->m;
            }
            break;   
        }
        case -1:
        {
            // LDPC object
            destroy_object<LDPC>(prhs[1]);
            break;   
        }
        /* Sum-Product */
        case 1:
        {
            // LDPC object
            LDPC& ldpc = get_object<LDPC>(prhs[1]);

            vector<double> in_llr(ldpc.n);
            vector<double> out_llr(ldpc.n);
            vector<int> y(ldpc.n);

            // in LLR
            input = mxGetPr(prhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                in_llr[i] = input[i];
            }

            // iterations
            int iMaxNumberOfIterations = (int)*(mxGetPr(prhs[3]));
            int number_of_iter = 0;

            bool result = true;
            result = SumProduct(ldpc, in_llr, iMaxNumberOfIterations, out_llr, y, &number_of_iter);

            /* denial flag */
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[0]))) = result;

            /* number of iterations */
            plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[1]))) = number_of_iter;

            double* data = NULL;

            /* hard desicion */
            plhs[2] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = y[i];
            }

            /* out LLR */
            plhs[3] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[3]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = out_llr[i];
            }
            break;
        }
        /* Layered Sum-Product */
        case 2:
        {
            // LDPC object
            LDPC& ldpc = get_object<LDPC>(prhs[1]);

            vector<double> in_llr(ldpc.n);
            vector<double> out_llr(ldpc.n);
            vector<int> row_sequence(ldpc.m);
            vector<int> y(ldpc.n);

            // in LLR
            input = mxGetPr(prhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                in_llr[i] = input[i];
            }

            // iterations
            int iMaxNumberOfIterations = (int)*(mxGetPr(prhs[3]));
            int number_of_iter = 0;

            // Row sequence (for layered decoder the row order is important)
            input = mxGetPr(prhs[4]);
            for(int i = 0; i < ldpc.m; ++i)
            {
                row_sequence[i] = (int)input[i];
            }

            bool result = true;
            result = LayeredSumProduct(ldpc, in_llr, iMaxNumberOfIterations, out_llr, y, &number_of_iter, row_sequence);

            /* denial flag */
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[0]))) = result;

            /* number of iterations */
            plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[1]))) = number_of_iter;

            double* data = NULL;

            /* hard desicion */
            plhs[2] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = y[i];
            }

            /* out LLR */
            plhs[3] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[3]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = out_llr[i];
            }
            break;
        }
        /* Min-Sum */
        case 3:
        {
            // LDPC object
            LDPC& ldpc = get_object<LDPC>(prhs[1]);

            vector<double> in_llr(ldpc.n);
            vector<double> out_llr(ldpc.n);
            vector<double> scale_array(ldpc.n);
            vector<double> offset_array(ldpc.n);
            vector<int> y(ldpc.n);

            // in LLR
            input = mxGetPr(prhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                in_llr[i] = input[i];
            }

            // iterations
            int iMaxNumberOfIterations = (int)*(mxGetPr(prhs[3]));
            int number_of_iter = 0;

            // scale array
            input = mxGetPr(prhs[4]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                scale_array[i] = input[i];
            }

            // offset array
            input = mxGetPr(prhs[5]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                offset_array[i] = input[i];
            }

            bool result = true;
            result = MinSum(ldpc, in_llr, iMaxNumberOfIterations, out_llr, y, &number_of_iter, scale_array, offset_array);

            /* denial flag */
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[0]))) = result;

            /* number of iterations */
            plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[1]))) = number_of_iter;

            double* data = NULL;

            /* hard desicion */
            plhs[2] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = y[i];
            }

            /* out LLR */
            plhs[3] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[3]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = out_llr[i];
            }
            break;
        }
        /* Layered Min-Sum */
        case 4:
        {
            // LDPC object
            LDPC& ldpc = get_object<LDPC>(prhs[1]);

            vector<double> in_llr(ldpc.n);
            vector<double> out_llr(ldpc.n);
            vector<int> row_sequence(ldpc.m);
            vector<double> scale_array(ldpc.n);
            vector<double> offset_array(ldpc.n);
            vector<int> y(ldpc.n);

            // in LLR
            input = mxGetPr(prhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                in_llr[i] = input[i];
            }

            // iterations
            int iMaxNumberOfIterations = (int)*(mxGetPr(prhs[3]));
            int number_of_iter = 0;

            // row sequence (for layered decoder the row order is important)
            input = mxGetPr(prhs[4]);
            for(int i = 0; i < ldpc.m; ++i)
            {
                row_sequence[i] = (int)input[i];
            }

            // scale array
            input = mxGetPr(prhs[5]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                scale_array[i] = input[i];
            }

            // offset array
            input = mxGetPr(prhs[6]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                offset_array[i] = input[i];
            }

            bool result = true;
            result = LayeredMinSum(ldpc, in_llr, iMaxNumberOfIterations, out_llr, y, &number_of_iter, row_sequence, scale_array, offset_array);

            /* denial flag */
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[0]))) = result;

            /* number of iterations */
            plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[1]))) = number_of_iter;

            double* data = NULL;

            /* hard desicion */
            plhs[2] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = y[i];
            }

            /* out LLR */
            plhs[3] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[3]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = out_llr[i];
            }
            break;
        }
        /* Adjusted Min-Sum */
        case 5:
        {
            // LDPC object
            LDPC& ldpc = get_object<LDPC>(prhs[1]);

            vector<double> in_llr(ldpc.n);
            vector<double> out_llr(ldpc.n);
            vector<int> y(ldpc.n);

            // in LLR
            input = mxGetPr(prhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                in_llr[i] = input[i];
            }

            // iterations
            int iMaxNumberOfIterations = (int)*(mxGetPr(prhs[3]));
            int number_of_iter = 0;

            bool result = true;
            result = AdjMinSum(ldpc, in_llr, iMaxNumberOfIterations, out_llr, y, &number_of_iter, 0);

            /* denial flag */
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[0]))) = result;

            /* number of iterations */
            plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[1]))) = number_of_iter;

            double* data = NULL;

            /* hard desicion */
            plhs[2] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = y[i];
            }

            /* out LLR */
            plhs[3] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[3]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = out_llr[i];
            }
            break;
        }
        /* Layered Adjusted Min-Sum */
        case 6:
        {
            // LDPC object
            LDPC& ldpc = get_object<LDPC>(prhs[1]);

            vector<double> in_llr(ldpc.n);
            vector<double> out_llr(ldpc.n);
            vector<int> row_sequence(ldpc.m);
            vector<int> y(ldpc.n);

            // in LLR
            input = mxGetPr(prhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                in_llr[i] = input[i];
            }

            // iterations
            int iMaxNumberOfIterations = (int)*(mxGetPr(prhs[3]));
            int number_of_iter = 0;

            // Row sequence (for layered decoder the row order is important)
            input = mxGetPr(prhs[4]);
            for(int i = 0; i < ldpc.m; ++i)
            {
                row_sequence[i] = (int)input[i];
            }

            bool result = true;
            result = LayeredAdjMinSum(ldpc, in_llr, iMaxNumberOfIterations, out_llr, y, &number_of_iter, row_sequence, 0);

            /* denial flag */
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[0]))) = result;

            /* number of iterations */
            plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[1]))) = number_of_iter;

            double* data = NULL;

            /* hard desicion */
            plhs[2] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = y[i];
            }

            /* out LLR */
            plhs[3] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[3]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = out_llr[i];
            }
            break;
        }
        /* Adjusted Min-Sum (Enhanced) */
        case 7:
        {
            // LDPC object
            LDPC& ldpc = get_object<LDPC>(prhs[1]);

            vector<double> in_llr(ldpc.n);
            vector<double> out_llr(ldpc.n);
            vector<int> y(ldpc.n);

            // in LLR
            input = mxGetPr(prhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                in_llr[i] = input[i];
            }

            // iterations
            int iMaxNumberOfIterations = (int)*(mxGetPr(prhs[3]));
            int number_of_iter = 0;

            bool result = true;
            result = AdjMinSum(ldpc, in_llr, iMaxNumberOfIterations, out_llr, y, &number_of_iter, 1);

            /* denial flag */
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[0]))) = result;

            /* number of iterations */
            plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[1]))) = number_of_iter;

            double* data = NULL;

            /* hard desicion */
            plhs[2] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = y[i];
            }

            /* out LLR */
            plhs[3] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[3]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = out_llr[i];
            }
            break;
        }
        /* Layered Adjusted Min-Sum (Enhanced) */
        case 8:
        {
            // LDPC object
            LDPC& ldpc = get_object<LDPC>(prhs[1]);

            vector<double> in_llr(ldpc.n);
            vector<double> out_llr(ldpc.n);
            vector<int> row_sequence(ldpc.m);
            vector<int> y(ldpc.n);

            // in LLR
            input = mxGetPr(prhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                in_llr[i] = input[i];
            }

            // iterations
            int iMaxNumberOfIterations = (int)*(mxGetPr(prhs[3]));
            int number_of_iter = 0;

            // Row sequence (for layered decoder the row order is important)
            input = mxGetPr(prhs[4]);
            for(int i = 0; i < ldpc.m; ++i)
            {
                row_sequence[i] = (int)input[i];
            }

            bool result = true;
            result = LayeredAdjMinSum(ldpc, in_llr, iMaxNumberOfIterations, out_llr, y, &number_of_iter, row_sequence, 1);

            /* denial flag */
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[0]))) = result;

            /* number of iterations */
            plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[1]))) = number_of_iter;

            double* data = NULL;

            /* hard desicion */
            plhs[2] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = y[i];
            }

            /* out LLR */
            plhs[3] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[3]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = out_llr[i];
            }
            break;
        }
        // Bit-flipping
        case 9:
        {
            LDPC& ldpc = get_object<LDPC>(prhs[1]);

            vector<int> x(ldpc.n);
            vector<int> y(ldpc.n);

            double* input = mxGetPr(prhs[2]);
            int iMaxNumberOfIterations = (int)*(mxGetPr(prhs[3]));
            
            for(int i = 0; i < ldpc.n; ++i)
            {
                x[i] = input[i];
            }
         
            int number_of_iter = 0;
         
            bool result = true;
            result = BitFlipping(ldpc, x, iMaxNumberOfIterations, y, &number_of_iter);
                     
            /* denial flag */
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[0]))) = result;
         
            /* number of iterations */
		    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
		    *((double*)(mxGetPr(plhs[1]))) = number_of_iter;
		 
		    double* data = NULL;
		 
		     /* hard desicion */
		     plhs[2] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
		     data = mxGetPr(plhs[2]);
		     for(int i = 0; i < ldpc.n; ++i)
		     {
			     data[i] = y[i];
		     }
    		 
		     break;
        }
        /* Sum-Product FFT */
        case 10:
        {
            // LDPC object
            LDPC& ldpc = get_object<LDPC>(prhs[1]);

            vector<double> in_llr(ldpc.n);
            vector<double> out_llr(ldpc.n);
            vector<int> y(ldpc.n);

            // in LLR
            input = mxGetPr(prhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                in_llr[i] = input[i];
            }

            // iterations
            int iMaxNumberOfIterations = (int)*(mxGetPr(prhs[3]));
            int number_of_iter = 0;

            bool result = true;
            result = SumProduct_FFT(ldpc, in_llr, iMaxNumberOfIterations, out_llr, y, &number_of_iter);

            /* denial flag */
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[0]))) = result;

            /* number of iterations */
            plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[1]))) = number_of_iter;

            double* data = NULL;

            /* hard desicion */
            plhs[2] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = y[i];
            }

            /* out LLR */
            plhs[3] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[3]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = out_llr[i];
            }
            break;
        }
        /* Min-Sum FFT*/
        case 11:
        {
            // LDPC object
            LDPC& ldpc = get_object<LDPC>(prhs[1]);

            vector<double> in_llr(ldpc.n);
            vector<double> out_llr(ldpc.n);
            vector<double> scale_array(ldpc.n);
            vector<double> offset_array(ldpc.n);
            vector<int> y(ldpc.n);

            // in LLR
            input = mxGetPr(prhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                in_llr[i] = input[i];
            }

            // iterations
            int iMaxNumberOfIterations = (int)*(mxGetPr(prhs[3]));
            int number_of_iter = 0;

            // scale array
            input = mxGetPr(prhs[4]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                scale_array[i] = input[i];
            }

            // offset array
            input = mxGetPr(prhs[5]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                offset_array[i] = input[i];
            }

            bool result = true;
            result = MinSum_FFT(ldpc, in_llr, iMaxNumberOfIterations, out_llr, y, &number_of_iter, scale_array, offset_array);

            /* denial flag */
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[0]))) = result;

            /* number of iterations */
            plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[1]))) = number_of_iter;

            double* data = NULL;

            /* hard desicion */
            plhs[2] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = y[i];
            }

            /* out LLR */
            plhs[3] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[3]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = out_llr[i];
            }
            break;
        }
        default:
        {
            break;
        }
    }
}
