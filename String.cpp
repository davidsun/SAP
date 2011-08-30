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




#include "String.h"

int String::dnaFormat(char *s){
	int ret = 0;
	toLower(s);
	for (int i = 0; s[i]; i ++){
		if (s[i] == 'a' || s[i] == 't' || s[i] == 'g' || s[i] == 'c') ret ++;
		else if (s[i] != 'n') return -1;
	}
	return ret;
}

int String::dnaFormat(DynamicArray <char> &s){
	int ret = 0;
	toLower(s.data());
	for (int i = 0; s[i]; i ++){
		if (s[i] == 'a' || s[i] == 't' || s[i] == 'g' || s[i] == 'c') ret ++;
		else if (s[i] != 'n') return -1;
	}
	return ret;
}

bool String::isDna(char *s){
	for (int i = 0; s[i]; i ++)
		if (s[i] != 'a' && s[i] != 't' && s[i] != 'g' && s[i] != 'c' && s[i] != 'n' && s[i] != 'A' && s[i] != 'T' && s[i] != 'G' && s[i] != 'C' && s[i] != 'N') return 0;
	return 1;
}

bool String::isDna(DynamicArray <char> &s){
	for (int i = 0; s[i]; i ++)
		if (s[i] != 'a' && s[i] != 't' && s[i] != 'g' && s[i] != 'c' && s[i] != 'n' && s[i] != 'A' && s[i] != 'T' && s[i] != 'G' && s[i] != 'C' && s[i] != 'N') return 0;
	return 1;
}

void String::reverseComplement(char *s, int size){
	char *p = s, *q = s + size - 1;
	for (; p < q; p ++, q --){
		char c = *p;
		*p = *q;
		*q = c;
	}
	for (int i = 0; i < size; i ++){
		if (s[i] == 'a') s[i] = 't';
		else if (s[i] == 't') s[i] = 'a';
		else if (s[i] == 'g') s[i] = 'c';
		else if (s[i] == 'c') s[i] = 'g';
	}
}


void String::toLower(char *s){
	for (int i = 0; s[i]; i ++)
		if (s[i] <= 'Z' && s[i] >= 'A') s[i] = s[i] - 'A' + 'a';
}
