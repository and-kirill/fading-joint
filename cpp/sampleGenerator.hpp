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
#ifndef SAMPLE_GENERATOR_H
#define SAMPLE_GENERATOR_H

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
using namespace boost;


class SampleGenerator
{
public:
    SampleGenerator(unsigned int seed)
    : m_gen(seed)
    {}
    double gen(unsigned int bins, double* pmf, double* values)
    {
        double sum = 0;
        double p = m_dis(m_gen);
        for (int i = 0; i < bins; ++i)
        {
            sum += pmf[i];
            if (p < sum)
            {
                return values[i];
            }
        }
        return pmf[bins-1];
    }
private:
    random::mt19937 m_gen;
    random::uniform_real_distribution<double> m_dis;
};

#endif
