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

#pragma once

#ifndef LOG_Q
#define LOG_Q					3
#endif
#define Q 						(1<<LOG_Q)

#include "GaloisField.h"

//using namespace std;

#if LOG_Q==1
#define POLYNOM Binary<11>::value
#elif LOG_Q==2
#define POLYNOM Binary<111>::value
#elif LOG_Q==3
#define POLYNOM Binary<1011>::value
#elif LOG_Q==4
#define POLYNOM Binary<10011>::value
#elif LOG_Q==5
#define POLYNOM Binary<100101>::value
#elif LOG_Q==6
#define POLYNOM Binary<1000011>::value
#elif LOG_Q==8
#define POLYNOM Binary<100011101>::value
#endif

typedef GaloisField<LOG_Q, POLYNOM> Field;
typedef GaloisFieldElement<Field> 	FieldElement;
