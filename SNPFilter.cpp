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
#include <stdlib.h>
#include <string.h>
#include <string>

#define MIN_READ_DEPTH 1
#define MAX_READ_DEPTH 10000
#define MAX_SCORE 5.00

using namespace std;

string inputFileName;
string outputFileName;
int minReadDepth = MIN_READ_DEPTH;
int maxReadDepth = MAX_READ_DEPTH;
double maxScore = MAX_SCORE;

void showWelcome(){
	fprintf(stderr, "SNP Filter - 0.99.85\n");
}

void showUsage(){
	fprintf(stderr, "	-i	Set input file name.\n");
	fprintf(stderr, "	-o	Set output file name.\n");
	fprintf(stderr, "	-m	Minimum depth of reads to call a SNP.\n");
	fprintf(stderr, "	-M	Maximum depth of reads to call a SNP.\n");
	fprintf(stderr, "	-s	Maximum score to call a SNP.\n");
	fprintf(stderr, "	-h	Show this help.\n");
}

bool processArguments(int argc, char **argv){
	char c;
	while ((c = getopt(argc, argv, "i:o:m:M:s:h")) != EOF){
		switch (c){
			case 'i':
				inputFileName = optarg;
				break;
			case 'o':
				outputFileName = optarg;
				break;
			case 'm':
				minReadDepth = atoi(optarg);
				break;
			case 'M':
				maxReadDepth = atoi(optarg);
				break;
			case 's':
				maxScore = atof(optarg);
				break;
			case 'h':
				return 1;
		}
	}
	
	if (inputFileName.empty()){
		fprintf(stderr, "ERROR: input file (-i) missing.\n");
		return 1;
	}
	
	if (outputFileName.empty()){
		fprintf(stderr, "ERROR: output file (-o) missing.\n");
		return 1;
	}
	
	fprintf(stderr, "	Input file name: %s\n", inputFileName.c_str());
	fprintf(stderr, "	Output file name: %s\n", outputFileName.c_str());
	fprintf(stderr, "	Minimum read depth: %d\n", minReadDepth);
	fprintf(stderr, "	Maximum read depth: %d\n", maxReadDepth);
	fprintf(stderr, "	Maximum score (lower score is better): %lf\n", maxScore);
	return 0;
}

int main(int argc, char **argv){
	showWelcome();
	if (processArguments(argc, argv)){
		showUsage();
		return 0;
	}
	FILE *f = fopen(inputFileName.c_str(), "r");
	if (!f){
		printf("ERROR: cannot open input file.\n");
		return 0;
	}
	FILE *fout = fopen(outputFileName.c_str(), "w");
	char s[1000], name[1000], type[1000];
	while (fgets(s, 1000, f)){
		s[strlen(s) - 1] = 0;
		int depth;
		double score;
		sscanf(s, "%s%s%lf%d", name, type, &score, &depth);
		if (type[0] <= '9' && type[0] >= '0' && score <= maxScore && depth >= minReadDepth && depth <= maxReadDepth) 
			fprintf(fout, "%s\n", s);
	}
	fclose(f);
	return 0;
}
