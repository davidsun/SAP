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





#ifndef DYNAMICARRAY_H
#define DYNAMICARRAY_H

template <class T>
class DynamicArray{
	public:
		DynamicArray(int baseSize = 100);
		DynamicArray(const DynamicArray <T> &array);
		~DynamicArray();
		
		T *data() const;
		void expand();
		void resize(unsigned size);
		unsigned size() const;
		T &operator [](int x);
		T operator [](int x) const;
		
	private:
		T *_arrayNodes;
		unsigned _size;
};

#include "DynamicArray.cpp"

#endif
