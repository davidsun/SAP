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


#include <stdio.h>
#include "String.h"
#include "IO.h"
#include "DynamicArray.h"

void trim(char *s){
	for (unsigned i = 0; s[i]; i ++)
		if (s[i] <= 'Z' && s[i] >= 'A') s[i] = s[i] - 'A' + 'a';
}

int main(int argc, char **argv){
	DynamicArray <char> name(100), read(100), sign(100), quality(100);
	IO::FileReader f(argv[1]);
	while (1){
		if (IO::readLine(f, name) == EOF) break;
		if (IO::readLine(f, read) == EOF) break;
		if (IO::readLine(f, sign) == EOF) break;
		if (IO::readLine(f, quality) == EOF) break;
		trim(name.data());
		trim(read.data());
		unsigned len = std::min(strlen(quality.data()), strlen(read.data()));
		read.data()[len] = 0;
		quality.data()[len] = 0;
		puts(name.data());
		puts(read.data());
		puts(quality.data());
	}
	return 0;
}
