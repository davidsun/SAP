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




#include "MatchStructures.h"
#include "IO.h"
#include "String.h"

#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <string.h>
#include <string>

#define DEFAULT_MIN_VALID_MATCH_COUNT 6
#define DEFAULT_MIN_QUALITY .4
#define DEFAULT_INSERTION_PREDICTION_SCORE .4
#define DEFAULT_DELETION_PREDICTION_SCORE .4
#define DEFAULT_MIN_READ_QUALITY .3
#define THETA 0.85
#define ETA 0.03

std::string inputFileName;
std::string referenceFileName;
std::string outputFileName;
int minValidMatchCount = DEFAULT_MIN_VALID_MATCH_COUNT;
double insertionPredictionScore = DEFAULT_INSERTION_PREDICTION_SCORE;
double deletionPredictionScore = DEFAULT_DELETION_PREDICTION_SCORE;
double minQuality = DEFAULT_MIN_QUALITY;
double minReadQuality = DEFAULT_MIN_READ_QUALITY;

#ifdef __cplusplus
extern "C" {
#endif
	long double expl(long double);
	long double logl(long double);
#ifdef __cplusplus
}
#endif

double LHet[256][256];
double Coef[64][256][256];

void calcHet(){
	for (int n1 = 0; n1 < 256; n1 ++)
		for (int n2 = 0; n2 < 256; n2 ++)
			LHet[n1][n2] = lgamma(n1 + n2 + 1) - lgamma(n1 + 1) - lgamma(n2 + 1) + logl(.5) * (n1 + n2);
}

void calcCoef(){
	long double sum_a[257], b[256], q_c[256], tmp[256], fk[256], fk2[256];
	double lC[256][256];

	fk[0] = fk2[0] = 1.0;
	for (int n = 1; n < 256; n ++){
		fk[n] = pow(THETA, n) * (1.0 - ETA) + ETA;
		fk2[n] = fk[n >> 1];
	}
	for (int n = 1; n < 256; n ++)
		for (int k = 1; k <= n; k ++)
			lC[n][k] = lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1);
	for (int q = 1; q < 64; q ++){
		double e = pow(10.0, - q / 10.0);
		double le = log(e);
		double le1 = log(1.0 - e);
		for (int n = 1; n < 256; n ++){
			sum_a[n + 1] = 0.0;
			for (int k = n; k >= 0; k --){
				sum_a[k] = sum_a[k + 1] + expl(lC[n][k] + k * le + (n - k) * le1);
				b[k] = sum_a[k + 1] / sum_a[k];
				if (b[k] > 0.99) b[k] = 0.99;
			}
			for (int k = 0; k < n; k ++) q_c[k] = - fk2[k] * logl(b[k] / e);
			for (int k = 1; k < n; k ++) q_c[k] += q_c[k - 1];
			for (int k = 0; k <= n; k ++){
				tmp[k] = - 4.343 *  logl(1.0 - expl(fk2[k] * logl(b[k])));
				Coef[q][n][k] = (k ? q_c[k - 1] : 0) + tmp[k];
			}
		}
	}
}

void showWelcome(){
	fprintf(stderr, "DNA Mutation Predictor - 0.99.85\n");
}

void showUsage(){
	fprintf(stderr, "	-i	Set input file name.\n");
	fprintf(stderr, "	-r	Set reference file name.\n");
	fprintf(stderr, "	-o	Set output file name.\n");
	fprintf(stderr, "	-m	Minimum read depth to call a SNP.\n");
	fprintf(stderr, "	-M	Maximum read depth to call a SNP.\n");
	fprintf(stderr, "	-q	Minimum quality.\n");
	fprintf(stderr, "	-Q	Minimum quality to validate a read.\n");
	fprintf(stderr, "	-I	Minimum score to call an insertion.\n");
	fprintf(stderr, "	-D	Minimum score to call a deletion.\n");
	fprintf(stderr, "	-h	Show this help.\n");
}

bool processArguments(int argc, char **argv){
	char c;
	while ((c = getopt(argc, argv, "i:r:o:m:S:I:D:q:Q:h")) != EOF){
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
			case 'm':
				minValidMatchCount = atoi(optarg);
				break;
			case 'I':
				insertionPredictionScore = atof(optarg);
				break;
			case 'D':
				deletionPredictionScore = atof(optarg);
				break;
			case 'q':
				minQuality = atof(optarg);
				break;
			case 'Q':
				minReadQuality = atof(optarg);
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
	
	if (referenceFileName.empty()){
		fprintf(stderr, "ERROR: reference file (-r) missing.\n");
		return 1;
	}
	
	fprintf(stderr, "	Input file name: %s\n", inputFileName.c_str());
	fprintf(stderr, "	Output file name: %s\n", outputFileName.c_str());
	fprintf(stderr, "	Reference file name: %s\n", referenceFileName.c_str());
	fprintf(stderr, "	Minimum valid match count: %d\n", minValidMatchCount);
	fprintf(stderr, "	Minimum quality: %lf\n", minQuality);
	fprintf(stderr, "	Minimum read quality: %lf\n", minReadQuality);
	fprintf(stderr, "	Deletion prediction score: %lf\n", deletionPredictionScore);
	fprintf(stderr, "	Insertion prediction score: %lf\n", insertionPredictionScore);
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

void processMatching(ExonList *exonList, const std::map <std::string, int> &exonNameToId){
	int dnaLen, bufferLen;
	DynamicArray <char> dna, buffer, quality;
	IO::BufferedFileReader *reader = IO::BufferedFileReader::newBufferedFileReader(inputFileName.c_str());
	std::list <mappingInfo> infos;
	double maxScore = -1;
	dnaLen = reader -> readLine(dna);
	reader -> readLine(quality);
	if (dnaLen == EOF) return;
	int valid = 0, invalid = 0;
	while ((bufferLen = reader -> readLine(buffer)) != EOF){
		if (bufferLen == 0){
			double totalScore = 0.0, totalQuality = 0.0;
			for (std::list <mappingInfo>::iterator it = infos.begin(); it != infos.end(); it ++)
				if (it -> exonId != -1 && it -> score >= maxScore * 0.9) totalScore += it -> score;
			for (unsigned i = 0; i < dnaLen; i ++){
				if (quality[i] >= 93) quality[i] = 93;
				quality[i] -= 33;
				totalQuality += quality[i] / 60.0;
			}
			if (1.0 / totalScore * maxScore >= minReadQuality && totalQuality / dnaLen >= minReadQuality){
				for (std::list <mappingInfo>::iterator it = infos.begin(); it != infos.end(); it ++){
					if (it -> exonId != -1 && it -> score < maxScore * 0.9) continue;
					int s = it -> strLoc, e = it -> exonLoc;
					std::string mappingString = it -> mappingString;
					ExonList::iterator itx = exonList -> exonById(it -> exonId);
					if (itx.isEnd()) continue;
					MatchExon *exon = (MatchExon *)itx.exon();
					if (it -> isReversed) String::reverseComplement(dna.data(), dnaLen);
					for (unsigned i = 0; i < mappingString.size(); i ++){
						if (e >= exon -> size()){
							//fprintf(stderr, "NO\n");
							break;
						}
						if (mappingString[i] == 'n' || mappingString[i] == 'c'){
							if (!it -> isReversed) exon -> updateMatchValue(e, dna[s], quality[s] / 60.0);
							else exon -> updateMatchValue(e, dna[s], quality[dnaLen - s - 1] / 60.0);
							s ++; e ++;
						}  else if (mappingString[i] == 'i'){
							int r = i, sp = s;
							while (r + 1 < mappingString.size() && mappingString[r + 1] == 'i') r ++, sp ++;
							if (!it -> isReversed) exon -> insert(e, dna.data() + sp, r - i + 1, quality[s] / 60.0);
							else exon -> insert(e, dna.data() + sp, r - i + 1, quality[dnaLen - s - 1] / 60.0);
							i = r;
							s = sp + 1;
						}  else {
							exon -> updateDeletionValue(e, totalQuality / dnaLen);
							e ++;
						}  
					}
					if (it -> isReversed) String::reverseComplement(dna.data(), dnaLen);
				}
				valid ++;
			}  else {
				invalid ++;
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
					if (matchPlace == 0){
						buffer[i] = 0;
						std::map <std::string, int>::const_iterator it = exonNameToId.find(buffer.data() + lastLoc);
						if (it == exonNameToId.end()) info.exonId = -1;
						else info.exonId = it -> second;
					}  else if (matchPlace == 1){
						if (buffer[lastLoc] == 'N') info.isReversed = 0;
						else info.isReversed = 1;
					}  else if (matchPlace == 2){
						sscanf(buffer.data() + lastLoc, "%d", &info.strLoc);
					}  else if (matchPlace == 3){
						sscanf(buffer.data() + lastLoc, "%d", &info.exonLoc);
					}  else if (matchPlace == 4){
						sscanf(buffer.data() + lastLoc, "%lf", &info.score);
						maxScore = std::max(maxScore, info.score);
					}  
					matchPlace ++;
					lastLoc = i + 1;
				}
			}
			info.mappingString = (buffer.data() + lastLoc);
			infos.push_back(info);
		}
	}
	fprintf(stderr, "Valid Reads: %d/%d\n", valid, valid + invalid);
}

void removeUnmappedExons(ExonList *exonList){
	double minScore = minValidMatchCount / 1.5;
	for (ExonList::iterator it = exonList -> begin(); !it.isEnd();){
		MatchExon *exon = (MatchExon *)it.exon();
		int matchLen = 0;
		for (int i = 0; i < exon -> size(); i ++){
			int insertionCount = 0;
			double insertionScore = 0.0;
			for (MatchExon::Insertion *ins = exon -> insertion(i); ins; ins = ins -> next()){
				insertionCount ++;
				insertionScore += ins -> score();
			}
			matchLen += ((exon -> deleteCount(i) >= minValidMatchCount && exon -> deleteScore(i) >= minScore) || 
						 (exon -> matchCount(i) >= minValidMatchCount && exon -> matchScore(i) >= minScore) ||
						 (insertionCount >= minValidMatchCount && insertionScore >= minScore));
		}
		if ((double)matchLen / exon -> size() < .7) it = exonList -> removeExon(it);
		else it ++;
	}
}

int main(int argc, char **argv){
	showWelcome();
	if (processArguments(argc, argv)){
		showUsage();
		return 0;
	}
	
	calcHet();
	calcCoef();
	
	ExonList *exonList = new ExonList;
	exonList -> readMatchExon(referenceFileName.c_str());
	std::map <std::string, int> exonNameToId;
	for (ExonList::iterator it = exonList -> begin(); !it.isEnd(); it ++)
		exonNameToId[it.exon() -> name()] = it.exon() -> id();
	processMatching(exonList, exonNameToId);
	//removeUnmappedExons(exonList);
	
	int snpCount = 0, deletionCount = 0, insertionCount = 0;
	FILE *fout = fopen(outputFileName.c_str(), "w");
	for (ExonList::iterator it = exonList -> begin(); !it.isEnd(); it ++){
		MatchExon *info = (MatchExon *)it.exon();
		bool found = 0;
		//printf("%s %d\n", info -> name(), info -> size());
		for (int i = 0; i < info -> size(); i ++){
			if (info -> matchCount(i) < minValidMatchCount) continue;
			printf("OK\n");
			char *dnaString = "atgc";
			int max[2] = {-1, -1}, c[3];
			double maxScore[2] = {0.0, 0.0}, q[3];
			for (int j = 0; j < 4; j ++){
				if (info -> matchScore(i, dnaString[j]) >= maxScore[0]){
					max[1] = max[0];
					max[0] = j;
					maxScore[1] = maxScore[0];
					maxScore[0] = info -> matchScore(i, dnaString[j]);
				}  else if (info -> matchScore(i, dnaString[j]) >= maxScore[1]){
					max[1] = j;
					maxScore[1] = info -> matchScore(i, dnaString[j]);
				}  
			}
			c[0] = max[0] >= 0 ? info -> matchCount(i, dnaString[max[0]]) : 0;
			c[1] = max[1] >= 0 ? info -> matchCount(i, dnaString[max[1]]) : 0;
			if ((c[2] = c[0] + c[1]) > 255){
				c[0] = int(255.0 * c[0] / c[2] + 0.5);
				c[1] = int(255.0 * c[1] / c[2] + 0.5);
				c[2] = 255;
			}
			/*double sum = log(0.001 * expl(LHet[c[1]][c[0]]) + (expl(Coef[c[2]][c[1]]) + expl(Coef[c[2]][c[0]])) * (1.0 - 0.001) / 2.0);
			q[0] = log((1.0 - 0.001) / 2.0) + Coef[c[2]][c[1]] - sum;
			q[1] = log((1.0 - 0.001) / 2.0) + Coef[c[2]][c[0]] - sum;*/
			int q0 = c[0] > 0 ? int(info -> matchScore(i, dnaString[max[0]]) / c[0] * 60.0 + .5) : 1;
			int q1 = c[1] > 0 ? int(info -> matchScore(i, dnaString[max[1]]) / c[1] * 60.0 + .5) : 1;
			q[0] = c[0] > 0 ? log((1.0 - 0.001) / 2.0) + Coef[q0][c[2]][c[1]] : -1e20;
			q[1] = c[1] > 0 ? log((1.0 - 0.001) / 2.0) + Coef[q1][c[2]][c[0]] : -1e20;
			q[2] = - 4.343 * log(2.0 * 0.001 / (1.0 - 0.001)) - 4.343 * LHet[c[1]][c[0]];
			printf("%lf %lf %lf\n", q[0], q[1], q[2]);
			//if (q[0] < 0.0) q[0] = 0.0;
			//if (q[1] < 0.0) q[1] = 0.0;
			/*printf("%d\t", info -> matchCount(i));
			for (int j = 0; j < 4; j ++){
				printf("%c/%d/%.4lf\t", dnaString[j], info -> matchCount(i, dnaString[j]), info -> matchScore(i, dnaString[j]));
			}
			printf("\n");*/
			char co = info -> dna()[i];
			if (q[0] >= q[1] && q[0] >= q[2]){
				if (info -> matchScore(i, dnaString[max[0]]) / info -> matchCount(i, dnaString[max[0]]) < minQuality) continue;
				if (co != dnaString[max[0]]){
					if (!found){
						found = 1;
						fprintf(fout, "%s\n", ((MatchExon *)info) -> name());
					}
					fprintf(fout, "%d\t%.3lf(%.3lf/%.3lf)\t%d/%d\t%c\t%c\n", i, 
							info -> matchScore(i, dnaString[max[0]]) / (info -> matchScore(i) + info -> deleteScore(i)), 
							info -> matchScore(i, dnaString[max[0]]), info -> matchScore(i) + info -> deleteScore(i), 
							info -> matchCount(i, dnaString[max[0]]), info -> matchCount(i),
							co, dnaString[max[0]]);
					snpCount ++;
				}
			}  else if (q[1] >= q[0] && q[1] >= q[2]){
				if (info -> matchScore(i, dnaString[max[1]]) / info -> matchCount(i, dnaString[max[1]]) < minQuality) continue;
				if (co != dnaString[max[1]]){
					//printf("OK\n");
					if (!found){
						found = 1;
						fprintf(fout, "%s\n", ((MatchExon *)info) -> name());
					}
					fprintf(fout, "%d\t%.3lf(%.3lf/%.3lf)\t%d/%d\t%c\t%c\n", i, 
							info -> matchScore(i, dnaString[max[1]]) / (info -> matchScore(i) + info -> deleteScore(i)), 
							info -> matchScore(i, dnaString[max[1]]), info -> matchScore(i) + info -> deleteScore(i), 
							info -> matchCount(i, dnaString[max[1]]), info -> matchCount(i),
							co, dnaString[max[1]]);
					snpCount ++;
				}
			}  else {
				if (info -> matchScore(i) / info -> matchCount(i) < minQuality) continue;
				if (!found){
					found = 1;
					fprintf(fout, "%s\n", ((MatchExon *)info) -> name());
				}
				fprintf(fout, "SNP\tLOC=%4d\tSCR=%.3lf(%.3lf/%.3lf)\tCHG=%c->%c/%c\n", i, 
						info -> matchScore(i, dnaString[max[1]]) / (info -> matchScore(i) + info -> deleteScore(i)), 
						info -> matchScore(i, dnaString[max[1]]), info -> matchScore(i) + info -> deleteScore(i), co, dnaString[max[0]], dnaString[max[1]]);
				snpCount ++;
			}
		}
	}
	//fprintf(stderr, "OK\n");
	
	for (ExonList::iterator it = exonList -> begin(); !it.isEnd(); it ++){
		MatchExon *info = (MatchExon *)it.exon();
		bool found = 0;
		for (int i = 0; i < info -> size(); i ++){
			double ms = info -> matchScore(i), ds = info -> deleteScore(i) * 1.3;
			if (info -> deleteCount(i) >= minValidMatchCount && ds / (ms + ds) >= deletionPredictionScore
				&& info -> deleteCount(i) + info -> matchCount(i) >= minValidMatchCount){
					if (!found){
						found = 1;
						fprintf(fout, "%s\n", ((MatchExon *)info) -> name());
					}
					fprintf(fout, "DEL\tLOC=%4d\tSCR=%.3lf(%.3lf/%.3lf)\n", i, (double)ds / (ms + ds), ds, ms + ds);
					deletionCount ++;
			}
		}
	}
	//fprintf(stderr, "OK\n");
	
	for (ExonList::iterator it = exonList -> begin(); !it.isEnd(); it ++){
		MatchExon *info = (MatchExon *)it.exon();
		bool found = 0;
		for (int i = 0; i < info -> size(); i ++){
			int total = 0;
			double totalScore = 0.0;
			std::map <std::string, double> insertion;
			std::map <int, double> lenPos;
			//fprintf(stderr, "OK1\n");
			for (MatchExon::Insertion *ins = info -> insertion(i); ins; ins = ins -> next()){
				insertion[ins -> dna()] += ins -> score();
				totalScore += ins -> score();
				lenPos[strlen(ins -> dna())] += ins -> score();
				total ++;
			}
			//fprintf(stderr, "OK2\n");
			totalScore *= 1.3;
			
			double scoreNear = 0.0;
			int matchNear = 0;
			if (i > 0){
				scoreNear = std::max(scoreNear, info -> matchScore(i - 1));
				matchNear = std::max(matchNear, info -> matchCount(i - 1));
			}
			if (i + 1 < info -> size()){
				scoreNear = std::max(scoreNear, info -> matchScore(i + 1));
				matchNear = std::max(matchNear, info -> matchCount(i + 1));
			}
			totalScore = std::min(totalScore, scoreNear);
			int mostProbLen = 1;
			for (std::map <int, double>::iterator it = lenPos.begin(); it != lenPos.end(); it ++)
				if (it -> second > lenPos[mostProbLen]) mostProbLen = it -> first;
			if (total >= minValidMatchCount && totalScore / scoreNear >= insertionPredictionScore && 
				total + matchNear >= minValidMatchCount){
					if (!found){
						found = 1;
						fprintf(fout, "%s\n", ((MatchExon *)info) -> name());
					}
					fprintf(fout, "INS\tLOC=%4d\tLEN=%d\tSCR=%.3lf(%.3lf/%.3lf)\tCHG=", i, mostProbLen,
							(double)totalScore / scoreNear, totalScore, scoreNear);
					for (std::map <std::string, double>::iterator it = insertion.begin(); it != insertion.end(); it ++)
						fprintf(fout, "%s(%.3lf) ", it -> first.c_str(), it -> second);
					fprintf(fout, "\n");
					insertionCount ++;
			}
		}
	}
	//fprintf(stderr, "OK\n");
	fclose(fout);
	fprintf(stderr, "Found SNP:%d, Deletion:%d, Insertion:%d\n", snpCount, deletionCount, insertionCount);
	return 0;
}
