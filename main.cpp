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
#include <list>
#include <cmath>
#include <string>
#include <algorithm>

#include <getopt.h>
#include <cstdlib>
#include <pthread.h>

#include "MatchStructures.h"
#include "MatchAlgoritms.h"
#include "IO.h"
#include "String.h"
#include "MatchHash.h"
#include <stdlib.h>

#define DEFAULT_MIN_QUALITY .90
#define DEFAULT_IS_FAST_MAP 0
#define DEFAULT_PIECE_SIZE 15
#define DEFAULT_THREAD_COUNT 1
#define DEFAULT_HASH_BINARY_SIZE 27
#define DEFAULT_CUT_COUNT 7
#define DEFAULT_MAXIMUM_GAP_RATIO 0.08
#define DEFAULT_INPUT_FILE_NAME "pieceOut.f"
#define DEFAULT_REFERENCE_FILE_NAME "templateOut.f"
#define DEFAULT_OUTPUT_FILE_NAME "result.out"

#define THREAD_OUTPUT_CACHE_SIZE 4194304
#define THREAD_OUTPUT_CACHE_BUFFER_SIZE 65536

const double FUNCTION_K = .99 / (log(.01) - log(1.01 - DEFAULT_MIN_QUALITY));
const double FUNCTION_B = 1.0 - FUNCTION_K * log(.01);
const int matchBonus = 20;
const int deletionPunishment = 13;

struct programParameter{
	std::string inputFileName;
	std::string referenceFileName;
	std::string outputFileName;
	float maximumGapRatio;
	bool isFastMap;
	int pieceSize;
	int threadCount;
	int hashBinarySize;
	int cutCount;
};

programParameter parameter = {DEFAULT_INPUT_FILE_NAME, DEFAULT_REFERENCE_FILE_NAME, DEFAULT_OUTPUT_FILE_NAME, 
								DEFAULT_MAXIMUM_GAP_RATIO,
								DEFAULT_IS_FAST_MAP, 
								DEFAULT_PIECE_SIZE, DEFAULT_THREAD_COUNT,
								DEFAULT_HASH_BINARY_SIZE, DEFAULT_CUT_COUNT};

/*
* static functions
*/
static inline int dnaToInt(char c){
	if (c >= 'A' && c <= 'Z') c = c - 'A' + 'a';
	switch (c){
		case 'a': return 0;
		case 't': return 1;
		case 'g': return 2;
		case 'c': return 3;
		default: return -1;
	}
}

/*
* Attention: 'n' in DNA is regarded as 'a'
*/
static inline unsigned long specialDnaToInt(char c){
	if (c >= 'A' && c <= 'Z') c = c - 'A' + 'a';
	switch (c){
		case 'n': return 0;
		case 'a': return 0;
		case 't': return 1;
		case 'g': return 2;
		case 'c': return 3;
		default: return -1;
	}
}

static inline int intOctLength(int x){
	if (x == 0) return 1;
	int ret = 0;
	for (; x; x /= 10) ret ++;
	return ret;
}

static inline int putInt(int x, char *s){
	int ret = intOctLength(x);
	for (int i = ret - 1; i >= 0; i --, x /= 10) s[i] = x % 10 + '0';
	return ret;
}

static inline int putUnitDouble(double x, char *s){
	if (x < 0.0) x = 0.0;
	s[0] = (int)x + '0';
	s[1] = '.';
	int t = (int)(x * 10000.0);
	if (t < 0) t = 0;
	for (int i = 5; i >= 2; i --, t /= 10) s[i] = t % 10 + '0';
	return 6;
}

static void addToHash(Hash *hash, Dna *dna, int l, int r, int segmentSize){
	unsigned long hashVal = 0;
	int nCount = 0;
	char *s = dna -> dna();
	if (segmentSize <= r - l + 1){
		for (int i = 0; i < segmentSize; i ++)
			if (s[l + i] == 'n') nCount ++;
		hashVal = Hash::calcHashValue(s + l, s + l + segmentSize, &specialDnaToInt);
		if (nCount <= 2) hash -> insert(dna, l, l + segmentSize, hashVal);
	}
	for (int i = l + 1; i + segmentSize <= r; i ++){
		if (s[i - 1] == 'n') nCount --;
		if (s[i + segmentSize - 1] == 'n') nCount ++;
		hashVal = (hashVal >> 2U) + (specialDnaToInt(s[i + segmentSize - 1]) << (segmentSize - 1 << 1ULL));
		if (nCount <= 2) hash -> insert(dna, i, i + segmentSize, hashVal);
	}
}

void processOneDnaDfs(char **ret, int **dp, const char *s1, const char *s2, int x, int y, int delta){
	if (x > 0 && dp[x - 1][y] == dp[x][y] && ret[x - 1][y] == -1){
		ret[x - 1][y] = 0;
		processOneDnaDfs(ret, dp, s1, s2, x - 1, y, delta);
	}
	if (x > 0 && s1[x - 1] == s2[x + y - delta - 1] && dp[x - 1][y] + matchBonus == dp[x][y] && ret[x - 1][y] == -1){
		ret[x - 1][y] = 0;
		processOneDnaDfs(ret, dp, s1, s2, x - 1, y, delta);
	}
	if (x > 0 && y < (delta << 1) && dp[x - 1][y + 1] - deletionPunishment == dp[x][y] && ret[x - 1][y + 1] == -1){
		ret[x - 1][y + 1] = 1;
		processOneDnaDfs(ret, dp, s1, s2, x - 1, y + 1, delta);
	}
	if (y > 0 && dp[x][y - 1] - deletionPunishment == dp[x][y] && ret[x][y - 1] == -1){
		ret[x][y - 1] = 2;
		processOneDnaDfs(ret, dp, s1, s2, x, y - 1, delta);
	}
}

bool processOneDna(ExonList *list, Hash *hash, bool isReversed, DynamicArray <char> &read, unsigned readSize, 
				   DynamicArray <char> &cache, int &cacheLoc, int **dp, char **next, bool fastMap){
	#define max(a, b) std::max(a, b)
	
	if (readSize < parameter.pieceSize) return 0;
	
	std::map <int, std::vector <int> > info;
	int lookUpLoc[parameter.cutCount];
	lookUpLoc[parameter.cutCount - 1] = readSize - parameter.pieceSize - 1;
	for (int i = 0; i < parameter.cutCount - 1; i ++) lookUpLoc[i] = (readSize - parameter.pieceSize) * i / (parameter.cutCount - 1);
	for (int i = 0; i < parameter.cutCount; i ++){
		int loc = lookUpLoc[i], nCount = 0;
		for (int j = 0; j < parameter.pieceSize; j ++)
			if (read[loc + j] == 'n') nCount ++;
		if (nCount > 2) continue;
		Hash::Result *res = hash -> exactFind(read.data() + loc, parameter.pieceSize);
		for (MatchHash::Result *p = res; p; p = p -> next) info[p -> seq].push_back(p -> start - loc);
		if (!fastMap && !res){
			MatchHash::Result *resx = hash -> oneMismatchFind(read.data() + loc, parameter.pieceSize);
			for (MatchHash::Result *p = resx; p; p = p -> next) info[p -> seq].push_back(p -> start - loc);
			MatchHash::deleteResult(resx);
		}
		Hash::deleteResult(res);
	}
	
	bool found = 0;
	double maxScore = 0.0;
	
	for (std::map <int, std::vector <int> >::iterator it = info.begin(); it != info.end(); it ++){
		std::vector <int> &dnaInfo = it -> second;
		std::sort(dnaInfo.begin(), dnaInfo.end());
		for (int i = 0; i < dnaInfo.size();){
			int r = i;
			int maxGapSize = (int)readSize * parameter.maximumGapRatio;
			while (r + 1 < dnaInfo.size() && dnaInfo[r + 1] - dnaInfo[i] < maxGapSize) r ++;
			if (r - i + 1 < 2){
				i = r + 1;
				continue;
			}
			
			Exon *dna = list -> exonById(it -> first).exon();
			unsigned dnaNameSize = strlen(dna -> name());
			if (dnaInfo[i] == dnaInfo[r] && dnaInfo[i] + readSize < dna -> size()){
				int matchLen = 0, left = dnaInfo[i];
				int bLeft = max(0, - dnaInfo[i]), bRight = std::min(readSize - 1, dna -> size() - dnaInfo[i] - 1);
				while (bLeft < readSize && read[bLeft] != (*dna)[left + bLeft]) bLeft ++;
				while (bRight >= bLeft && read[bRight] != (*dna)[left + bRight]) bRight --;
				for (int j = bLeft; j <= bRight; j ++) matchLen += (read[j] == (*dna)[left + j]);
				int length = bRight - bLeft + 1;
				if (matchLen / (double)readSize >= DEFAULT_MIN_QUALITY){
					double score = 1.0 - (1.0 - matchLen / (double)readSize) / (1.0 - DEFAULT_MIN_QUALITY);
					memcpy(cache.data() + cacheLoc, dna -> name(), dnaNameSize); cacheLoc += dnaNameSize;
					cache[cacheLoc ++] = '\t'; 
					if (isReversed) cache[cacheLoc ++] = 'R';
					else cache[cacheLoc ++] = 'N';
					cache[cacheLoc ++] = '\t'; 
					cacheLoc += putInt(bLeft, cache.data() + cacheLoc);
					cache[cacheLoc ++] = '\t';
					cacheLoc += putInt(left + bLeft, cache.data() + cacheLoc);
					cache[cacheLoc ++] = '\t';
					cacheLoc += putUnitDouble(score, cache.data() + cacheLoc);
					cache[cacheLoc ++] = '\t';
					for (int j = bLeft; j <= bRight; j ++){
						if (read[j] == (*dna)[left + j]) cache[cacheLoc ++] = 'n';
						else cache[cacheLoc ++] = 'c';
					}	
					cache[cacheLoc ++] = '\n';
					if (cacheLoc >= cache.size() - THREAD_OUTPUT_CACHE_BUFFER_SIZE)
						cache.resize(cache.size() + THREAD_OUTPUT_CACHE_BUFFER_SIZE);
					found = 1;
				}
			}  else {
				int left = dnaInfo[i], right = std::min(dnaInfo[r] + readSize, dna -> size() - 1) + 1;
				int bLeft = max(0, - dnaInfo[i]);
				int p1 = 0, p2 = 0, delta = dnaInfo[r] - dnaInfo[i], bd = delta << 1, maxK2 = right - left;
				for (int j = bLeft; j <= readSize; j ++){
					memset(dp[j], 128, sizeof(int) * ((maxGapSize << 1) + 5));
					memset(next[j], 255, (maxGapSize << 1) + 5);
				}
				for (int k = 0; k <= bd; k ++) dp[bLeft][k] = 0;
				for (int j = bLeft; j <= readSize; j ++, maxK2 --){
					for (int k = 0; k <= bd && k <= maxK2; k ++){
						int pv = dp[j][k], dnaPos = left + j + k - delta;
						if (k > 0) dp[j + 1][k - 1] = max(dp[j + 1][k - 1], pv - deletionPunishment);
						dp[j][k + 1] = max(dp[j][k + 1], pv - deletionPunishment);
						if (j < readSize && dnaPos < right && dnaPos >= 0){
							if (read[j] == (*dna)[dnaPos]) dp[j + 1][k] = max(dp[j + 1][k], pv + matchBonus);
							else dp[j + 1][k] = max(dp[j + 1][k], pv);
						}
						if (pv > dp[p1][p2]) p1 = j, p2 = k;
					}
				}
				
				while (p1 > 0 && dp[p1][p2] == dp[p1 - 1][p2]) p1 --;
				while (p2 > 0 && dp[p1][p2] - deletionPunishment == dp[p1][p2 - 1]) p2 --;
				while (p1 > 0 && p2 < bd && dp[p1][p2] - deletionPunishment == dp[p1 - 1][p2 + 1]) p1 --, p2 ++;
				
				processOneDnaDfs(next, dp, read.data(), dna -> dna() + left, p1, p2, delta);
				int s1 = bLeft, s2 = 0;
				while (next[s1][s2] == -1) s2 ++;
				while (s1 < readSize && dp[s1 + 1][s2] == dp[s1][s2] && next[s1 + 1][s2] != -1) s1 ++;
				while (s2 > 0 && dp[s1 + 1][s2 - 1] == dp[s1][s2] - deletionPunishment && next[s1 + 1][s2 - 1] != -1)
					s1 ++, s2 --;
				while (s2 < bd && dp[s1][s2 + 1] == dp[s1][s2] - deletionPunishment && next[s1][s2 + 1] != -1) s2 ++;
				
				int length = std::min(p1 - s1, p1 + p2 - (s1 + s2));
				double quality = dp[p1][p2] / (double)matchBonus / readSize;
				if (quality >= DEFAULT_MIN_QUALITY){
					double score = 1.0 - (1.0 - quality) / (1.0 - DEFAULT_MIN_QUALITY);
					memcpy(cache.data() + cacheLoc, dna -> name(), dnaNameSize); cacheLoc += dnaNameSize;
					cache[cacheLoc ++] = '\t'; 
					if (isReversed) cache[cacheLoc ++] = 'R';
					else cache[cacheLoc ++] = 'N';
					cache[cacheLoc ++] = '\t'; 
					cacheLoc += putInt(s1, cache.data() + cacheLoc);
					cache[cacheLoc ++] = '\t';
					cacheLoc += putInt(left + s1 + s2 - delta, cache.data() + cacheLoc);
					cache[cacheLoc ++] = '\t';
					cacheLoc += putUnitDouble(score, cache.data() + cacheLoc);
					cache[cacheLoc ++] = '\t';
					while (s1 != p1 || s2 != p2){
						if (next[s1][s2] == 0){
							if (dp[s1 + 1][s2] == dp[s1][s2]) cache[cacheLoc ++] = 'c';
							else cache[cacheLoc ++] = 'n';
							s1 ++;
						}  else if (next[s1][s2] == 1){
							cache[cacheLoc ++] = 'i';
							s1 ++; s2 --;
						}  else {
							cache[cacheLoc ++] = 'd';
							s2 ++;
						}
					}
					cache[cacheLoc ++] = '\n';
					if (cacheLoc >= cache.size() - THREAD_OUTPUT_CACHE_BUFFER_SIZE)
						cache.resize(cache.size() + THREAD_OUTPUT_CACHE_BUFFER_SIZE);
					found = 1;
				}
			}
			
			i = r + 1;
		}
	}
	return found;
}

struct threadedProcessDnaArg{
	IO::FileReader *reader;
	IO::FileWriter *writer;
	ExonList *list;
	Hash *hash;
};

struct threadedProcessResult{
	int dnaFound, dnaTotal;
};

void *threadedProcessDna(void *arg){
	threadedProcessDnaArg *args = (threadedProcessDnaArg *)arg;
	DynamicArray <char> name(1000), dna(1000), quality(1000);
	threadedProcessResult *ret = new threadedProcessResult;
	ret -> dnaFound = ret -> dnaTotal = 0;
	
	DynamicArray <char> cache(THREAD_OUTPUT_CACHE_SIZE + THREAD_OUTPUT_CACHE_BUFFER_SIZE);
	int cacheLoc = 0;
	int currentDnaMaxLength = 128;
	int maxGapSize = (int)(parameter.maximumGapRatio * currentDnaMaxLength);
	int **dp = MatchAlgorithms::create2DimArray(currentDnaMaxLength, (maxGapSize << 1) + 5, 0);
	char **next = MatchAlgorithms::create2DimArray(currentDnaMaxLength, (maxGapSize << 1) + 5, (char)0);
	
	while (1){
		std::pair <int, int> size;
		if ((size = args -> reader -> readExon(name, dna, quality)).first == EOF) break;
		
		while (size.second + 5 > currentDnaMaxLength){
			MatchAlgorithms::erase2DimArray(dp, currentDnaMaxLength);
			MatchAlgorithms::erase2DimArray(next, currentDnaMaxLength);
			currentDnaMaxLength <<= 1;
			int maxGapSize = (int)(parameter.maximumGapRatio * currentDnaMaxLength);
			dp = MatchAlgorithms::create2DimArray(currentDnaMaxLength, (maxGapSize << 1) + 5, 0);
			next = MatchAlgorithms::create2DimArray(currentDnaMaxLength, (maxGapSize << 1) + 5, (char)0);
		}
		bool dnaFound = 0;
		
		memcpy(cache.data() + cacheLoc, dna.data(), size.second); cacheLoc += size.second; 
		cache[cacheLoc ++] = '\n';
		memcpy(cache.data() + cacheLoc, quality.data(), size.second); cacheLoc += size.second; 
		cache[cacheLoc ++] = '\n';
		dnaFound |= processOneDna(args -> list, args -> hash, 0, dna, size.second, cache, cacheLoc, dp, next, parameter.isFastMap);
		String::reverseComplement(dna.data(), size.second);
		dnaFound |= processOneDna(args -> list, args -> hash, 1, dna, size.second, cache, cacheLoc, dp, next, parameter.isFastMap);
		String::reverseComplement(dna.data(), size.second);
		
		if (!dnaFound) cacheLoc -= ((size.second + 1) << 1);
		else cache[cacheLoc ++] = '\n';
		if (cacheLoc >= THREAD_OUTPUT_CACHE_SIZE){
			args -> writer -> putString(cache.data(), cacheLoc);
			cacheLoc = 0;
		}
		
		if (dnaFound) ret -> dnaFound ++;
		ret -> dnaTotal ++;
	}
	if (cacheLoc) args -> writer -> putString(cache, cacheLoc);
	MatchAlgorithms::erase2DimArray(dp, currentDnaMaxLength);
	MatchAlgorithms::erase2DimArray(next, currentDnaMaxLength);
	pthread_exit((void *)ret);
}

int processDna(ExonList *list, Hash *hash, const char *inputFileName, const char *outputFileName, int threadCount){
	IO::BufferedFileReader *reader = IO::BufferedFileReader::newBufferedFileReader(inputFileName);
	IO::BufferedFileWriter *writer = IO::BufferedFileWriter::newBufferedFileWriter(outputFileName);
	if (!reader -> isOpen()){
		delete reader;
		return 1;
	}
	
	pthread_t *threads = (pthread_t *)malloc(sizeof(pthread_t) * threadCount);
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	threadedProcessDnaArg *processDnaArg = new threadedProcessDnaArg[threadCount];
	for (int i = 0; i < threadCount; i ++){
		processDnaArg[i].reader = reader;
		processDnaArg[i].writer = writer;
		processDnaArg[i].hash = hash;
		processDnaArg[i].list = list;
	}
	
	for (int i = 0; i < threadCount; i ++)
		pthread_create(&threads[i], &attr, threadedProcessDna, (void *)(processDnaArg + i));
	
	int found = 0, total = 0;
	for (int i = 0; i < threadCount; i ++){
		void *status;
		pthread_join(threads[i], &status);
		
		threadedProcessResult *res = (threadedProcessResult *)status;
		found += res -> dnaFound;
		total += res -> dnaTotal;
		delete res;
	}
	delete []processDnaArg;
	
	fprintf(stderr, "\nProcessing finished. Found %d in %d (%lf).\n", found, total, (double)found / total);
	
	free(threads);
	pthread_attr_destroy(&attr);
	delete writer;
	delete reader;
	return 0;
}

void showWelcome(){
	fprintf(stderr, "DNA Mapper - 0.99.85\n");
}

void showUsage(){
	fprintf(stderr, "\t-f\tEnable FASTMAP mapping mode.\n");
	fprintf(stderr, "\t-t\tSet thread count. (Default: 1)\n");
	fprintf(stderr, "\t-H\tSet the binary size of hash, usually between 20 and 30 (Default: 27)\n");
	fprintf(stderr, "\t-C\tSet the number of pieces that a read is cut into, usually between 7 and 30 (Default: 7)\n");
	fprintf(stderr, "\t-G\tSet the maximum gap ratio, which is MaxGapLength/SequenceLength (Default: 0.1)\n");
	fprintf(stderr, "\t-i\tSet input file name. (Default: pieceOut.f)\n");
	fprintf(stderr, "\t-r\tSet reference file name. (Default: templateOut.f)\n");
	fprintf(stderr, "\t-o\tSet output file name. (Default: result.out)\n");
	fprintf(stderr, "\t-p\tSet the size of small pieces when mapping, usually between 10 and 16. (Default: 15)\n");
	fprintf(stderr, "\t-h\tShow this help.\n");
}

bool processArguments(int argc, char **argv){
	char c;
	while ((c = getopt(argc, argv, "H:C:G:t:fi:r:o:p:h")) != EOF){
		switch (c){
			case 'H':
				parameter.hashBinarySize = atoi(optarg);
				break;
			case 'C':
				parameter.cutCount = atoi(optarg);
				break;
			case 'G':
				parameter.maximumGapRatio = atof(optarg);
				break;
			case 't':
				parameter.threadCount = atoi(optarg);
				break;
			case 'f':
				parameter.isFastMap = 1;
				break;
			case 'i':
				parameter.inputFileName = optarg;
				break;
			case 'r':
				parameter.referenceFileName = optarg;
				break;
			case 'o':
				parameter.outputFileName = optarg;
				break;
			case 'p':
				parameter.pieceSize = atoi(optarg);
				break;
			case 'h':
				return 1;
		}
	}
	
	if (parameter.hashBinarySize < 20 || parameter.hashBinarySize > 32){
		printf("ERROR: Hash binary size should between 20 and 32.\n");
		return 1;
	}
	
	if (parameter.pieceSize > 20 || parameter.pieceSize < 10){
		printf("ERROR: Piece size should between 10 and 20.\n");
		return 1;
	}
	
	if (parameter.maximumGapRatio < 0.0 || parameter.maximumGapRatio > .5){
		printf("ERROR: Gap Ratio should between 0.0 and 0.5.\n");
		return 1;
	}
	
	if (parameter.cutCount > 30 || parameter.cutCount < 7){
		printf("ERROR: Cut count should between 7 and 15.\n");
		return 1;
	}
	
	fprintf(stderr, "\tInput file name: %s\n", parameter.inputFileName.c_str());
	fprintf(stderr, "\tReference file name: %s\n", parameter.referenceFileName.c_str());
	fprintf(stderr, "\tOutput file name: %s\n", parameter.outputFileName.c_str());
	fprintf(stderr, "\tHash size: %llu\n", 1ULL << parameter.hashBinarySize);
	fprintf(stderr, "\tCut count: %d\n", parameter.cutCount);
	fprintf(stderr, "\tThread count: %d\n", parameter.threadCount);
	fprintf(stderr, "\tPiece size: %d\n", parameter.pieceSize);
	fprintf(stderr, "\tMaximum gap ratio: %.4f\n", parameter.maximumGapRatio);
	if (parameter.isFastMap) fprintf(stderr, "	FASTMAP enabled.\n");
	return 0;
}

int main(int argc, char **argv){
	showWelcome();
	if (processArguments(argc, argv)){
		showUsage();
		return 0;
	}
	
	ExonList *exonList = new ExonList;
	if (exonList -> readExon(parameter.referenceFileName.c_str())){
		fprintf(stderr, "Cannot open reference file: %s.\n", parameter.referenceFileName.c_str());
		exit(1);
	}
	BufferedBinaryHash *hashExon = new BufferedBinaryHash(exonList -> totalExonSize(), parameter.hashBinarySize);
	for (ExonList::iterator it = exonList -> begin(); !it.isEnd(); it ++){
		Dna *dna = it.exon();
		addToHash(hashExon, dna, 0, dna -> size() - 1, parameter.pieceSize);
	}
	processDna(exonList, hashExon, parameter.inputFileName.c_str(), parameter.outputFileName.c_str(), parameter.threadCount);
	
	return 0;
}
