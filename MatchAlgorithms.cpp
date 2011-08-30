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


#include <algorithm>

namespace MatchAlgorithms{
	template <class T>
	T **create2DimArray (int dim1, int dim2, T initVal){
		T **ret = new T*[dim1];
		for (int i = 0; i < dim1; i ++){
			ret[i] = new T[dim2];
			for (int j = 0; j < dim2; j ++) ret[i][j] = initVal;
		}
		return ret;
	}
	
	template <class T>
	void erase2DimArray(T **a, int dim1){
		for (int i = 0; i < dim1; i ++) delete[] a[i];
		delete[] a;
	}
	
	int calcMaxMatch(char *s1, int len1, char *s2, int len2){
		int i, j;
		int **f = create2DimArray(len1 + 1, len2 + 1, 0);
		for (i = 0; i < len1; i ++){
			for (j = 0; j < len2; j ++){
				f[i + 1][j] = std::max(f[i + 1][j], f[i][j]);
				f[i][j + 1] = std::max(f[i][j + 1], f[i][j]);
				if (s1[i] == s2[j]) f[i + 1][j + 1] = std::max(f[i + 1][j + 1], f[i][j] + 1);
			}
		}
		int ret = std::max(f[len1][len2], std::max(f[len1 - 1][len2], f[len1][len2 - 1]));
		erase2DimArray(f, len1 + 1);
		return ret;
	}
	
	int calcMatch(const char *s1, const char *s2, int len){
		int ret = 0;
		for (int i = 0; i < len; i ++, s1 ++, s2 ++)
			ret += (*s1 == *s2);
		return ret;
	}
};
