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




#ifndef STRING_H
#define STRING_H

#include "DynamicArray.h"

namespace String{
	/*
	* change invalid char in DNA into 'n', and return the valid length of a DNA
	*/
	int dnaFormat(char *s);
	int dnaFormat(DynamicArray <char> &s);
	bool isDna(char *s);
	bool isDna(DynamicArray <char> &s);
	
	void reverseComplement(char *s, int size);
	
	void toLower(char *s);
};

#endif
