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
#ifdef OMP
#include <omp.h>
#endif
#include "common.h"
#include "func_nodes.h"
#include "ObjectHandle.h"
#include "Matrix.h"
#include "sampleGenerator.hpp"

const double MAX_EXP = exp((double)700);
const double MIN_EXP = exp((double)-700);

void dec_func_nodes_x(FuncConnections& fc, 
                      double* outLLR, 
                      double* inLLR, 
                      complex<double>* y, 
                      double sigma, 
                      unsigned int samples,
                      Matrix<double>& h_Q_re,
                      Matrix<double>& h_Q_im,
                      double* pdf_interval,
                      unsigned int seed);
void value2cwd(int M, int k, int value, vector<int>& cwd);
double safe_exp(double x);
double safe_log(double x);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 10)
    {
        mexErrMsgTxt("usage: decode_func_nodes_x(<fc handle>, <in_llr>, <rx>, <sigma>, <num_samples>, <bins>, <pmf_re>, <pmf_im>, <pmf_interval>, <seed>)");
    }
    
    FuncConnections& fc = get_object<FuncConnections>(prhs[0]);
	int T = fc.m_T;
    int N = fc.m_N;
    int dmax = fc.m_dmax;
    
    double* inLLR = (double*)mxGetData(prhs[1]);

    double* yR = (double*)mxGetPr(prhs[2]);
    double* yI = (double*)mxGetPi(prhs[2]);
	
    complex<double>* y = new complex<double>[T];
    for (int i = 0; i < T; ++i)
    {
        if (yI)
        {
            y[i] = complex<double>(yR[i], yI[i]);
        }
        else
        {
            y[i] = complex<double>(yR[i], 0);
        }
    }

    double sigma = mxGetScalar(prhs[3]);
    
    int samples = (int)*(mxGetPr(prhs[4]));
    
    unsigned int bins = (int)*(mxGetPr(prhs[5]));
    
    Matrix<double> h_Q_re(N, bins);
    double* input = (double*)mxGetPr(prhs[6]);
   
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < bins; ++j)
        {
            h_Q_re(i, j) = input[N*j + i];
        }
    }
    
    Matrix<double> h_Q_im(N, bins);
    input = (double*)mxGetPr(prhs[7]);
   
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < bins; ++j)
        {
            h_Q_im(i, j) = input[N*j + i];
        }
    }
    
    input = (double*)mxGetPr(prhs[8]);
   
    double* pdf_interval = new double[bins];
    for (int i = 0; i < bins; ++i)
    {
        pdf_interval[i] = input[i];
    }
    
    unsigned int seed = (int)*(mxGetPr(prhs[8]));
    
    double* outLLR = new double[N];
    dec_func_nodes_x(fc, outLLR, inLLR, y, sigma, samples, h_Q_re, h_Q_im, pdf_interval, seed);


    plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
	double* output = (double*)mxGetData(plhs[0]);

	for (int i=0; i < N; ++i)
	{
	    output[i] = outLLR[i];
	}

	delete[] y;
	delete[] outLLR;
    delete[] pdf_interval;
}

void dec_func_nodes_x(FuncConnections& fc, 
                      double* outLLR, 
                      double* inLLR, 
                      complex<double>* y, 
                      double sigma, 
                      unsigned int samples,
                      Matrix<double>& h_Q_re,
                      Matrix<double>& h_Q_im,
                      double* pdf_interval,
                      unsigned int seed)
{
    int N = fc.m_N; // N = K*n
    int T = fc.m_T; // number of functional nodes
    int* rowWeights = fc.m_rowWeights;
    Matrix<int> S_positions(fc.m_S_positions);
    Matrix<complex<double> > S_values(fc.m_S_values);
    int dmax = fc.m_dmax;
    int M = fc.m_M;
    complex<double>* modTable = fc.m_modTable;
    unsigned int bins = h_Q_re.getColumnsNumber();
    double max_value = pdf_interval[bins-1];
    SampleGenerator sampleGenerator(seed);
    double h_sample_re = 0;
    double h_sample_im = 0;
    
    
    Matrix<double> inLogProb(2, N);
    for (int j = 0; j < N; ++j)
    {
        inLogProb(0, j) = inLLR[j]/2;
        inLogProb(1, j) = -inLLR[j]/2;
    }

    Matrix<double> outProb(2, N);
    for (int j = 0; j < N; ++j)
    {
        outProb(0, j) = 0;
        outProb(1, j) = 0;
    }

#ifdef OMP    
    unsigned int CPU_N = omp_get_max_threads();
    //mexPrintf("%d\n", CPU_N);
#else
    unsigned int CPU_N = 1;
#endif
    
#ifdef OMP
    #pragma omp parallel num_threads(CPU_N) shared(N, T, inLogProb, y, sigma, samples, pdf_interval, rowWeights, S_positions, S_values, dmax, M, modTable, bins, max_value, outProb, seed)
#endif
    {
#ifdef OMP
        int threadID = omp_get_thread_num();
#else
        int threadID = 0;
#endif
        vector<int> cwd(dmax);
        double h_sample_re = 0;
        double h_sample_im = 0;
        SampleGenerator sampleGenerator(seed + threadID);
        
        Matrix<double> outProb_local(outProb);
        
        int s = 0; // current number of samples

        while (s < samples/CPU_N + 1) 
        {
            // we need to average over fading samples
            for (int i = 0; i < T; ++i)
            {
                // decode one functional node
                int dc = rowWeights[i];
                for (int value = 0; value < pow((double)M, (double)dc); ++value)
                {
                    // check all the codewords
                    value2cwd(M, dc, value, cwd);
                    complex<double> signal_sum(0, 0);
                    for (int t = 0; t < dc; ++t)
                    {
                        h_sample_re = sampleGenerator.gen(bins, h_Q_re.getRowElements(S_positions(i, t)), pdf_interval);
                        h_sample_im = sampleGenerator.gen(bins, h_Q_im.getRowElements(S_positions(i, t)), pdf_interval);
                        complex<double> h_sample(h_sample_re, h_sample_im);
                        signal_sum += S_values(i, t)*h_sample*modTable[cwd[t]];
                    }
                    // P(y | sum)
                    double dist = abs(y[i]-signal_sum);

                    double log_P = - dist*dist/2/sigma/sigma;
                    // Prod p
                    double log_p = 0;
                    for (int t = 0; t < dc; ++t)
                    {
                        log_p += inLogProb(cwd[t], S_positions(i, t));
                    }
                    double Pp = safe_exp(log_P + log_p);

                    for (int t = 0; t < dc; ++t)
                    {
                        outProb_local(cwd[t], S_positions(i, t)) += Pp;
                    }
                }
            }
            ++s;
        }
#ifdef OMP            
        #pragma omp critical
#endif
        {
            for (int j = 0; j < N; ++j)
            {
                outProb(0, j) += outProb_local(0, j);
                outProb(1, j) += outProb_local(1, j);
            }
        }
    }
    for (int i = 0; i < N; ++i)
    {
        outLLR[i] = safe_log(outProb(0, i)) - safe_log(outProb(1, i));
    }
}

void value2cwd(int M, int k, int value, vector<int>& cwd)
{
    int temp = value;
    for (int i = 0; i < k; ++i)
    {
        cwd[i] = temp % M;
        temp = temp/M;
    }
}

double safe_exp(double x)
{
    if (x > 700)
    {
        return MAX_EXP;
    }
    if (x < -700)
    {
        return MIN_EXP;
    }
    return exp(x);
}
double safe_log(double x)
{
    if (x < MIN_EXP)
    {
        return -700;
    }
    return log(x);
}
