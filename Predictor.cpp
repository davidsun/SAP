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
#include <string>
#include <vector>
#include <algorithm>

#define DEFAULT_MIN_VALID_MATCH_COUNT 5
#define DEFAULT_INSERTION_PREDICTION_SCORE .4
#define DEFAULT_DELETION_PREDICTION_SCORE .4
#define DEFAULT_MIN_READ_QUALITY .3
#define PR .0001

std::string inputFileName;
std::string referenceFileName;
std::string outputFileName;
int minValidMatchCount = DEFAULT_MIN_VALID_MATCH_COUNT;
double minReadQuality = DEFAULT_MIN_READ_QUALITY;

double LHet[256][256];
double KM1[256], KM2[256];

std::pair <double, double> LogP(double q){
	double val = pow10(- q / 10.0);
	return std::make_pair(log(1.0 - val), log(val));
}

void preCalc(){
	for (int n1 = 0; n1 < 256; n1 ++)
		for (int n2 = 0; n2 < 256; n2 ++)
			LHet[n1][n2] = lgamma(n1 + n2 + 1) - lgamma(n1 + 1) - lgamma(n2 + 1) + logl(.5) * (n1 + n2);
	for (int i = 0; i < 256; i ++){
		KM1[i] = LogP(i).first;
		KM2[i] = LogP(i).second;
	}
}

void showWelcome(){
	fprintf(stderr, "DNA Mutation Predictor - 0.99.85\n");
}

void showUsage(){
	fprintf(stderr, "	-i	Set input file name.\n");
	fprintf(stderr, "	-r	Set reference file name.\n");
	fprintf(stderr, "	-o	Set output file name.\n");
	fprintf(stderr, "	-Q	Minimum quality to validate a read.\n");
	fprintf(stderr, "	-m	Mimimum depth of reads to validate a insertion/deletion/SNP.\n");
	fprintf(stderr, "	-h	Show this help.\n");
}

bool processArguments(int argc, char **argv){
	char c;
	while ((c = getopt(argc, argv, "i:r:o:m:I:D:Q:h")) != EOF){
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
				if (quality[i] <= 0) quality[i] = 1;
				totalQuality += quality[i];
			}
			if (totalQuality / dnaLen >= minReadQuality){
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
							break;
						}
						if (mappingString[i] == 'n' || mappingString[i] == 'c'){
							if (!it -> isReversed) exon -> updateMatchValue(e, dna[s], KM1[quality[s]] - KM2[quality[s]], KM2[quality[s]]);
							else exon -> updateMatchValue(e, dna[s], KM1[quality[dnaLen - s - 1]] - KM2[quality[dnaLen - s - 1]], KM2[quality[dnaLen - s - 1]]);
							s ++; e ++;
						}  else if (mappingString[i] == 'i'){
							int r = i, sp = s;
							double totalQ = (it -> isReversed) ? quality[dnaLen - s - 1] : quality[s];
							while (r + 1 < mappingString.size() && mappingString[r + 1] == 'i'){
								r ++, sp ++;
								totalQ += (it -> isReversed) ? quality[dnaLen - sp - 1] : quality[sp];
							}
							if (!it -> isReversed) exon -> insert(e, dna.data() + s, r - i + 1, LogP(totalQ / (r - i + 1)).first);
							else exon -> insert(e, dna.data() + s, r - i + 1, LogP(totalQ / (r - i + 1)).first);
							i = r;
							s = sp + 1;
						}  else {
							if (!it -> isReversed) exon -> updateDeletionValue(e, KM1[quality[s]] - KM2[quality[s]], KM2[quality[s]]);
							else exon -> updateDeletionValue(e, KM1[quality[dnaLen - s - 1]] - KM2[quality[dnaLen - s - 1]], KM2[quality[dnaLen - s - 1]]);
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
	preCalc();
	
	showWelcome();
	if (processArguments(argc, argv)){
		showUsage();
		return 0;
	}
	
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
		//printf("%s %d\n", info -> name(), info -> size());
		for (int i = 0; i < info -> size(); i ++){
			if (info -> matchCount(i) < minValidMatchCount) continue;
			std::vector <std::pair <double, int> > score;
			for (int j = 0; j < 4; j ++)
				if (info -> matchCount(i, dnaString[j]) > 0)
					score.push_back(std::make_pair(- info -> matchCount(i, dnaString[j]), j));
			std::sort(score.begin(), score.end());
			/*printf("%d\t%d\t\t", i, info -> matchCount(i));
			for (int j = 0; j < score.size(); j ++){
				int cur = dnaString[score[j].second];
				printf("%c/%.4lf/%d/%.4lf\t", cur, - score[j].first, info -> matchCount(i, cur), info -> matchScore(i, cur));
			}
			printf("\n");*/
			if (score.size() == 0) continue;
			char co = info -> dna()[i];
			int c1 = dnaString[score[0].second], c2 = -1;
			double tc = 0.0;
			if (score.size() >= 2){
				c2 = dnaString[score[1].second];
				int count1 = - score[0].first, count2 = - score[1].first;
				int sum = - score[0].first - score[1].first;
				/*if (sum > 255){
					count1 = int((double)count1 / sum * 255.0 + .5);
					count2 = int((double)count2 / sum * 255.0 + .5);
					sum = 255;
				}*/
				//printf("OK\n");
				//printf("%d %d %d\n", sum, count1, count2);
				double pp1 = info -> matchScore(i, c1) + info -> totalQ(i);
				double pp2 = info -> matchScore(i, c2) + info -> totalQ(i);
				double pp3 = lgamma(sum + 1) - lgamma(count1 + 1) - lgamma(count2 + 1) + log(.5) * (sum);
				double div = log(PR * exp(pp3) + (1.0 - PR) / 2.0 * (exp(pp1) + exp(pp2)));
				double p1 = pp1 + log((1.0 - PR) / 2.0) - div;
				double p2 = pp2 + log((1.0 - PR) / 2.0) - div;
				double p3 = log(PR) + pp3;
				//printf("%lf %lf %lf\n", pp1, pp2, pp3);
				//printf("%lf\n", PR * exp(pp3) + (1.0 - PR) / 2.0 * (exp(pp1) + exp(pp2)));
				//printf("%d %d %lf %lf %lf %lf %lf\n", count1, count2, div, info -> matchScore(i, c1), info -> matchScore(i, c2), info -> totalQ(i));
				if (p1 >= p2 && p1 >= p3){
					c2 = -1;
					tc = fabs(p1 * p1 / p2 / p3);
				}  else if (p2 >= p1 && p2 >= p3){
					c1 = c2, c2 = -1;
					tc = fabs(p2 * p2 / p1 / p3);
				}  else tc = fabs(p3 * p3 / p1 / p2);
			} 
			if (c2 == -1){
				int cur = c1;
				if (co != cur){
					fprintf(fout, "%s\t%d\t%.3lf\t%d\t%d\t%c\t%c\n", ((MatchExon *)info) -> name() + 1, i, 
							tc * 1000.0, info -> matchCount(i, cur), info -> matchCount(i),
							co, cur);
					snpCount ++;
				}
			}  else {
				fprintf(fout, "%s\t%d\t%.3lf\t%d\t%d\t%c\t%c%c\n", ((MatchExon *)info) -> name() + 1, i, 
						tc * 1000.0, info -> matchCount(i, c1) + info -> matchCount(i, c2), info -> matchCount(i),
						co, c1, c2);
				snpCount ++;
			}  
		}
	}
	
	for (ExonList::iterator it = exonList -> begin(); !it.isEnd(); it ++){
		MatchExon *info = (MatchExon *)it.exon();
		for (int i = 0; i < info -> size(); i ++){
			double ms = info -> matchScore(i) + info -> totalQ(i), ds = info -> deleteScore(i) + info -> totalQ(i);
			if (info -> deleteCount(i) >= minValidMatchCount && ds >= ms){
				fprintf(fout, "%s\tDEL\t%d\t%d\t%.3lf\t%.3lf\n", ((MatchExon *)info) -> name() + 1, i, info -> deleteCount(i), ds, ms);
				deletionCount ++;
			}
		}
	}
	
	for (ExonList::iterator it = exonList -> begin(); !it.isEnd(); it ++){
		MatchExon *info = (MatchExon *)it.exon();
		bool found = 0;
		for (int i = 0; i < info -> size(); i ++){
			int total = 0, len = -1;
			double totalScore = 0.0;
			std::map <int, double> scores;
			std::map <std::string, double> insertion;
			for (MatchExon::Insertion *ins = info -> insertion(i); ins; ins = ins -> next()){
				insertion[ins -> dna()] += ins -> score();
				totalScore += ins -> score();
				scores[strlen(ins -> dna())] += ins -> score();
				total ++;
			}
			
			double scoreNear = info -> matchScore(i) + info -> totalQ(i);
			if (i + 1 < info -> size()) scoreNear = std::min(scoreNear, info -> matchScore(i + 1) + info -> totalQ(i + 1));
			if (total >= minValidMatchCount && totalScore >= scoreNear){
				for (std::map <int, double>::iterator it = scores.begin(); it != scores.end(); it ++)
					if (len == -1 || it -> second > scores[len]) len = it -> first;
				fprintf(fout, "%s\tINS\t%d\t%d\t%.3lf\t%.3lf\tCHG=", ((MatchExon *)info) -> name() + 1, i, total, totalScore, scoreNear);
				for (std::map <std::string, double>::iterator it = insertion.begin(); it != insertion.end(); it ++)
					fprintf(fout, "%s(%.3lf) ", it -> first.c_str(), it -> second);
				fprintf(fout, "\n");
				insertionCount ++;
			}
		}
	}
	fclose(fout);
	fprintf(stderr, "Found SNP:%d, Deletion:%d, Insertion:%d\n", snpCount, deletionCount, insertionCount);
	return 0;
}
