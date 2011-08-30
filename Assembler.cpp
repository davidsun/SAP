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





#include <map>
#include <vector>
#include <string>
#include <algorithm>

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MatchStructures.h"
#include "MatchAlgoritms.h"
#include "IO.h"
#include "String.h"
#include "MatchHash.h"

#define DEFAULT_MIN_READ_QUALITY .3

using namespace std;

string inputFileName;
string referenceFileName;
string outputFileName;
double minReadQuality = DEFAULT_MIN_READ_QUALITY;

void showWelcome(){
	fprintf(stderr, "DNA Assembler - 0.99.85\n");
}

void showUsage(){
	fprintf(stderr, "	-i	Set input file name.\n");
	fprintf(stderr, "	-r	Set reference file name.\n");
	fprintf(stderr, "	-o	Set output file name.\n");
	fprintf(stderr, "	-q	Minimum quality to validate a read.\n");
	fprintf(stderr, "	-h	Show this help.\n");
}

bool processArguments(int argc, char **argv){
	char c;
	while ((c = getopt(argc, argv, "i:r:o:q:h")) != EOF){
		switch (c){
			case 'i':
				inputFileName = optarg;
				break;
			case 'r':
				referenceFileName = optarg;
				break;
			case 'o':
				outputFileName = optarg;
				break;
			case 'q':
				minReadQuality = atof(optarg);
				break;
			case 'h':
				return 1;
		}
	}
	
	fprintf(stderr, "	Input file name: %s\n", inputFileName.c_str());
	fprintf(stderr, "	Reference file name: %s\n", referenceFileName.c_str());
	fprintf(stderr, "	Output file name: %s\n", outputFileName.c_str());
	fprintf(stderr, "	Minimum read quality: %lf\n", minReadQuality);
	return 0;
}

struct mappingInfo{
	int exonId;
	double score;
	std::string mappingString;
	int strLoc, exonLoc;
	char isReversed;
	
	bool operator < (const mappingInfo &m) const{
		return score > m.score;
	}
};

void sortReads(string filename){
	vector <FILE *> files;
	
	int dnaLen, bufferLen;
	DynamicArray <char> dna, buffer, quality;
	IO::BufferedFileReader *reader = IO::BufferedFileReader::newBufferedFileReader(filename.c_str());
	double maxScore = -1;
	dnaLen = reader -> readLine(dna);
	reader -> readLine(quality);
	if (dnaLen == EOF) return;
	while ((bufferLen = reader -> readLine(buffer)) != EOF){
		if (bufferLen == 0){
			double totalScore = 0.0, totalQuality = 0.0;
			for (std::list <mappingInfo>::iterator it = infos.begin(); it != infos.end(); it ++)
				if (it -> exonId != -1 && it -> score >= maxScore * 0.9) totalScore += it -> score;
			for (unsigned i = 0; i < dnaLen; i ++){
				if (quality[i] >= 93) quality[i] = 93;
				quality[i] -= 33;
			}
			if (1.0 / totalScore * maxScore >= minReadQuality){
				for (std::list <mappingInfo>::iterator it = infos.begin(); it != infos.end(); it ++){
					if (it -> exonId != -1 && it -> score < maxScore * 0.9) continue;
					int s = it -> strLoc, e = it -> exonLoc;
					std::string mappingString = it -> mappingString;
					ExonList::iterator itx = exonList -> exonById(it -> exonId);
					if (itx.isEnd()) continue;
					MatchExon *exon = (MatchExon *)itx.exon();
					if (it -> isReversed){
						String::reverseComplement(dna.data(), dnaLen);
						
					}
					
					if (it -> isReversed) String::reverseComplement(dna.data(), dnaLen);
				}
			}
			infos.clear();
			dnaLen = reader -> readLine(dna);
			if (dnaLen == EOF) break;
			reader -> readLine(quality);
			maxScore = -1;
		}  else {
			mappingInfo info;
			int lastLoc = 0, matchPlace = 0;
			for (int i = 0; i < bufferLen; i ++){
				if (buffer[i] == '\t'){
				}
			}
			info.mappingString = (buffer.data() + lastLoc);
			infos.push_back(info);
		}
	}
}

int main(int argc, char **argv){
	showWelcome();
	if (processArguments(argc, argv)){
		showUsage();
		return 0;
	}
	
	
	
	return 0;
}
