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
include common.mk

L = -g -lm -lpthread

MapperO = IO.o MatchStructures.o main.o String.o MatchHash.o MatchTrie.o 
PredictorO = Predictor.o MatchStructures.o IO.o String.o
FastqToFDQO = FastqToFDQ.o String.o IO.o MatchStructures.o
FastaToFDAO = FastaToFDA.o String.o IO.o MatchStructures.o
SNPFilterO = SNPFilter.o
IndelFilterO = IndelFilter.o

main:   ${MapperO} ${FastqToFDQO} ${FastaToFDAO} ${PredictorO} ${SNPFilterO} ${IndelFilterO}
	${CPP} ${COPT} ${CFLAGS} -o ${BINDIR}/Mapper ${MapperO} $(MYLIBS) $L
	${STRIP} ${BINDIR}/Mapper
	${CPP} ${COPT} ${CFLAGS} -o ${BINDIR}/FastqToFDQ ${FastqToFDQO} $(MYLIBS) $L
	${STRIP} ${BINDIR}/FastqToFDQ
	${CPP} ${COPT} ${CFLAGS} -o ${BINDIR}/FastaToFDA ${FastaToFDAO} $(MYLIBS) $L
	${STRIP} ${BINDIR}/FastaToFDA
	${CPP} ${COPT} ${CFLAGS} -o ${BINDIR}/Predictor ${PredictorO} $(MYLIBS) $L
	${STRIP} ${BINDIR}/Predictor
	${CPP} ${COPT} ${CFLAGS} -o ${BINDIR}/SNPFilter ${SNPFilterO} $(MYLIBS) $L
	${STRIP} ${BINDIR}/SNPFilter
	${CPP} ${COPT} ${CFLAGS} -o ${BINDIR}/IndelFilter ${IndelFilterO} $(MYLIBS) $L
	${STRIP} ${BINDIR}/IndelFilter

clean:
	rm -f *.o
	rm -f bin/${MACHTYPE}/Mapper
	rm -f bin/${MACHTYPE}/FastaToFDA
	rm -f bin/${MACHTYPE}/Predictor
	rm -f bin/${MACHTYPE}/PredictorBeyes

all:
	cd ../lib && ${MAKE}
	make


