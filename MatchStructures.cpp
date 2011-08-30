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




#include <string.h>
#include <vector>
#include <map>
#include <stdio.h>
#include <stdlib.h>

#include "IO.h"
#include "MatchStructures.h"
#include "String.h"

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
static inline int specialDnaToInt(char c){
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

static inline unsigned calcHashValue(const char *start, const char *end, int (*funcDnaToInt)(char c) = 0){
	unsigned ret = 0;
	if (funcDnaToInt){
		for (const char *c = end - 1; c >= start; c --) ret = (ret << 2U) + funcDnaToInt(*c);
	}  else {
		for (const char *c = end - 1; c >= start; c --) ret = (ret << 2U) + dnaToInt(*c);
	}
	return ret;
}

extern "C"{
	static inline int atomic_exchange(int *dest, int value){
		int ret;
		__asm__ __volatile__("xchg %0,(%2)" : "=q"(ret) : "0"(value), "r"(dest) : "memory");
		return ret;
	}
}

/*
* Dna
*/
int Dna::totalDnaCount = 0;

Dna::Dna(const Dna &dna) : _size(dna._size), _id(dna._id){
	_dna = new char[_size + 1];
	memcpy(_dna, dna._dna, _size);
	_dna[_size] = 0;
}

Dna::Dna(char *dna, int size) : _size(size), _id(totalDnaCount ++){
	_dna = new char[_size + 1];
	memcpy(_dna, dna, _size);
	_dna[_size] = 0;
}

Dna::~Dna(){
	delete[] _dna;
}

int Dna::id() const{
	return _id;
}

unsigned Dna::size() const{
	return _size;
}

char *Dna::dna() const{
	return _dna;
}

char &Dna::operator[](int x){
	return _dna[x];
}

char Dna::operator[](int x) const{
	return _dna[x];
}


/*
* Exon
*/
Exon::Exon(const Exon& exon): Dna(exon){
	int len = strlen(exon._name);
	_name = new char[len + 1];
	memcpy(_name, exon._name, len);
	_name[len] = 0;
}

Exon::Exon(char *name, char *dna, int dnaSize) : Dna(dna, dnaSize){
	int len = strlen(name);
	_name = new char[len + 1];
	memcpy(_name, name, len);
	_name[len] = 0;
}

Exon::~Exon(){
	delete[] _name;
}

char* Exon::name() const{
	return _name;
}

/*
* MatchExon
*/
MatchExon::Insertion::Insertion(char *dna, unsigned short size, double score) : _size(size), _score(score){
	_dna = new char[_size + 1];
	memcpy(_dna, dna, _size);
	_dna[_size] = 0;
	_next = 0;
}

MatchExon::Insertion::~Insertion(){
	delete[] _dna;
}

char *MatchExon::Insertion::dna() const{
	return _dna;
}


MatchExon::Insertion *MatchExon::Insertion::next() const{
	return _next;
}

double MatchExon::Insertion::score() const{
	return _score;
}

MatchExon::MatchExon(const Exon &exon): Exon(exon){
	init();
}

MatchExon::MatchExon(char *name, char *exon, int size): Exon(name, exon, size){
	init();
}

MatchExon::~MatchExon(){
	delete[] _deletionCount;
	delete[] _deletionValue;
	
	for (int i = 0; i < _size; i ++){
		Insertion *next;
		for (Insertion *now = _insertion[i]; now; now = next){
			next = now -> _next;
			delete now;
		}
	}
	delete[] _insertion;
	
	for (int i = 0; i < 4; i ++){
		delete[] _matchCount[i];
		delete[] _matchValue[i];
	}
	delete[] _totalQ;
	
	delete[] _lock;
}

int MatchExon::deleteCount(int loc) const{
	return _deletionCount[loc];
}

double MatchExon::deleteScore(int loc) const{
	return _deletionValue[loc];
}

void MatchExon::insert(int loc, char* s, int len, double score){
	Insertion *insertion = new Insertion(s, len, score);
	insertion -> _next = _insertion[loc];
	_insertion[loc] = insertion;
}

MatchExon::Insertion *MatchExon::insertion(int loc){
	return _insertion[loc];
}

int MatchExon::matchCount(int loc) const{
	return _matchCount[0][loc] + _matchCount[1][loc] + _matchCount[2][loc] + _matchCount[3][loc];
}

int MatchExon::matchCount(int loc, char c) const{
	int v = dnaToInt(c);
	if (v == -1) return 0.0;
	return _matchCount[v][loc];
}

double MatchExon::matchScore(int loc) const{
	double ret = 0.0;
	for (int i = 0; i < 4; i ++) ret += _matchValue[i][loc];
	return ret;
}

double MatchExon::matchScore(int loc, char c) const{
	int v = dnaToInt(c);
	if (v == -1) return 0.0;
	return _matchValue[v][loc];
}

char MatchExon::mostProbableDna(int loc) const{
	char ret = 'n';
	double maxValue = -1.0;
	for (int i = 0; i < 4; i ++){
		if (_matchValue[i][loc] > maxValue){
			maxValue = _matchValue[i][loc];
			ret = dnaString[i];
		}
	}
	return ret;
}

void MatchExon::updateDeletionValue(int loc, double score, double quality){
	_deletionCount[loc] ++;
	_deletionValue[loc] += score;
	_totalQ[loc] += quality;
}

void MatchExon::updateMatchValue(int loc, char c, double score, double quality){
	int v = dnaToInt(c);
	if (v == -1) return;
	_matchCount[v][loc] ++;
	_matchValue[v][loc] += score;
	_totalQ[loc] += quality;
}

void MatchExon::init(){
	_deletionCount = new short[_size];
	memset(_deletionCount, 0, sizeof(short) * _size);
	
	_deletionValue = new double[_size];
	for (int i = 0; i < _size; i ++) _deletionValue[i] = 0.0;
	
	_insertion = new Insertion*[_size + 1];
	memset(_insertion, 0, sizeof(Insertion *) * (_size + 1));
	
	_totalQ = new double[_size];
	for (int i = 0; i < _size; i ++) _totalQ[i] = 0.0;
	
	for (int i = 0; i < 4; i ++){
		_matchCount[i] = new short[_size];
		_matchValue[i] = new double[_size];
		memset(_matchCount[i], 0, sizeof(short) * _size);
		for (int j = 0; j < _size; j ++) _matchValue[i][j] = 0.0;
	}
	
	_lock = new int[(_size >> 6) + 2];
	memset(_lock, 0, sizeof(int) * ((_size >> 6) + 2));
}

void MatchExon::lock(int loc, int len){
	for (int i = loc >> 6; (i << 6) <= loc + len; i ++)
		while (atomic_exchange(_lock + i, 1));
}

double MatchExon::totalQ(int loc) const{
	return _totalQ[loc];
}

void MatchExon::unlock(int loc, int len){
	for (int i = loc >> 6; (i << 6) <= loc + len; i ++) _lock[i] = 0;
}

/*
* ExonList
*/
ExonList::ExonList() : _totalExonSize(0){
}

ExonList::ExonList(const std::list <Exon *> &list){
	for (std::list <Exon *>::const_iterator it = list.begin(); it != list.end(); it ++){
		Exon *exon = new Exon((*it) -> name(), (*it) -> dna(), (*it) -> size());
		_infos[exon -> id()] = exon;
		_totalExonSize += exon -> size();
	}
}

ExonList::iterator ExonList::exonById(int id){
	return iterator(this, _infos.find(id));
}

ExonList::~ExonList(){
	for (std::map <int, Exon *>::iterator it = _infos.begin(); it != _infos.end(); it ++)
		delete it -> second;
}

bool ExonList::iterator::isEnd(){
	return _it == _parent -> _infos.end();
}

void ExonList::iterator::operator ++(){
	_it ++;
}

void ExonList::iterator::operator ++(int){
	_it ++;
}

ExonList::iterator::iterator(ExonList *parent, const std::map <int, Exon *>::iterator &it) : _parent(parent), _it(it){
}

ExonList::iterator ExonList::begin(){
	return iterator(this, _infos.begin());
}

Exon* ExonList::iterator::exon() const{
	return _it -> second;
}

void ExonList::readExon(const char* fileName){
	IO::BufferedFileReader *reader = IO::BufferedFileReader::newBufferedFileReader(fileName);
	if (!reader -> isOpen()) return;
	DynamicArray <char> bufferName(100), bufferDna(1000);
	int sizeName, sizeDna;
	while ((sizeName = reader -> readLine(bufferName) >= 0) && (sizeDna = reader -> readLine(bufferDna)) >= 0){
		String::toLower(bufferDna.data());
		Exon *exon = new Exon(bufferName.data(), bufferDna.data(), sizeDna);
		_infos[exon -> id()] = exon;
		_totalExonSize += exon -> size();
	}
	delete reader;
}

void ExonList::readMatchExon(const char* fileName){
	IO::BufferedFileReader *reader = IO::BufferedFileReader::newBufferedFileReader(fileName);
	if (!reader -> isOpen()) return;
	DynamicArray <char> bufferName(100), bufferDna(1000);
	int sizeName, sizeDna;
	while ((sizeName = reader -> readLine(bufferName) >= 0) && (sizeDna = reader -> readLine(bufferDna)) >= 0){
		String::toLower(bufferDna.data());
		Exon *exon = new MatchExon(bufferName.data(), bufferDna.data(), sizeDna);
		_infos[exon -> id()] = exon;
		_totalExonSize += exon -> size();
	}
	delete reader;
}

ExonList::iterator ExonList::removeExon(const ExonList::iterator &it){
	ExonList::iterator ret = it;
	_totalExonSize -= it._it -> second -> size();
	ret ++;
	delete it._it -> second;
	_infos.erase(it._it);
	return ret;
}

int ExonList::totalExonSize() const{
	return _totalExonSize;
}

