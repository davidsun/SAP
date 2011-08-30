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





#include <string.h>
#include <algorithm>

/*
* DynamicArray <T>
*/
template <class T>
DynamicArray <T>::DynamicArray(int baseSize) : _size(baseSize){
	_arrayNodes = new T[_size];
	memset(_arrayNodes, 0, sizeof(T) * _size);
}

template <class T>
DynamicArray <T>::DynamicArray(const DynamicArray <T> &array) : _size(array._size){
	memcpy(_arrayNodes, array._arrayNodes, sizeof(T) * _size);
}

template <class T>
DynamicArray <T>::~DynamicArray(){
	delete[] _arrayNodes;
}

template <class T>
T *DynamicArray <T>::data() const{
	return _arrayNodes;
}

template <class T>
void DynamicArray <T>::expand(){
	T *newNodes = new T[_size << 1];
	memcpy(newNodes, _arrayNodes, sizeof(T) * _size);
	memset(newNodes + _size, 0, sizeof(T) * _size);
	delete[] _arrayNodes;
	_arrayNodes = newNodes;
	_size <<= 1;
}

template <class T>
void DynamicArray <T>::resize(unsigned int size){
	T *newNodes = new T[size];
	memcpy(newNodes, _arrayNodes, sizeof(T) * std::min(size, _size));
	if (size > _size) memset(newNodes + _size, 0, sizeof(T) * (size - _size));
	delete[] _arrayNodes;
	_arrayNodes = newNodes;
	_size = size;
}


template <class T>
T &DynamicArray <T>::operator [](int x){
	return _arrayNodes[x];
}

template <class T>
T DynamicArray <T>::operator [](int x) const{
	return _arrayNodes[x];
}

template <class T>
unsigned DynamicArray <T>::size() const{
	return _size;
}
