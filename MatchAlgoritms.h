/*******************************************************************************
 * This file is part of SAP.
 * 
 * SAP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * SAP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SAP.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/




#ifndef MATCHALGORITHMS_H
#define MATCHALGORITHMS_H

#include "MatchStructures.h"

namespace MatchAlgorithms{
	template <class T> T **create2DimArray(int dim1, int dim2, T initVal);
	
	template <class T> void erase2DimArray(T **a, int dim1);
	
	static int calcMaxMatch(char *s1, int len1, char *s2, int len2);
	
	static int calcMatch(const char *s1, const char *s2, int len);
};

#include "MatchAlgorithms.cpp"

#endif
