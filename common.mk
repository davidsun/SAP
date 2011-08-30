#-------------------------------------------------------------------------------
# This file is part of SAP.
# 
# SAP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# SAP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with SAP.  If not, see <http://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------------
CC=gcc
CPP=g++
MACHTYPE=x86_64
ifeq (${COPT},)
    COPT=-O3
endif
CFLAGS=-g
HG_DEFS=-DMACHTYPE_${MACHTYPE}
HG_WARN=-Wformat -Wimplicit -Wuninitialized -Wreturn-type

ifeq (${BINDIR},)
    BINDIR = bin/${MACHTYPE}
endif
MKDIR=mkdir -p
ifeq (${STRIP},)
   STRIP=strip
endif

# location of stringify program
STRINGIFY = ${BINDIR}/stringify

%.o: %.c
	${CC} ${COPT} ${CFLAGS} ${HG_DEFS} ${HG_WARN} ${HG_INC} ${XINC} -o $@ -c $<

%.o: %.cpp
	${CPP} ${COPT} ${CFLAGS} ${HG_DEFS} ${HG_WARN} ${HG_INC} ${XINC} -o $@ -c $<
