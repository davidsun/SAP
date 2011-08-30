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



#include "MatchHash.h"

/*
* static functions
*/
static inline unsigned long dnaToInt(char c){
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

/*
 * Hash
 */
Hash::Result::Result(const Hash::Result &res) : start(res.start), seq(res.seq){
}

Hash::Result::Result(int seq, int start) : start(start), seq(seq){
}

Hash::Hash(){
}

Hash::~Hash(){
}

unsigned long Hash::calcHashValue(const char *start, const char *end, unsigned long (*funcDnaToInt)(char c)){
	unsigned long ret = 0;
	if (funcDnaToInt){
		for (const char *c = end - 1; c >= start; c --) ret = (ret << 2U) + funcDnaToInt(*c);
	}  else {
		for (const char *c = end - 1; c >= start; c --) ret = (ret << 2U) + dnaToInt(*c);
	}
	return ret;
}

void Hash::deleteResult(MatchHash::Result *&res){
	Result *next;
	for (; res; res = next){
		next = res -> next;
		delete res;
	}
}

/*
* MatchHash
*/
MatchHash::MatchHash(){
}

MatchHash::~MatchHash(){
}

Hash::Result *MatchHash::exactFind(const char *s, unsigned len) const{
	Result *ret = 0;
	unsigned long hashVal = Hash::calcHashValue(s, s + len, &specialDnaToInt);
	find_p(hashVal, ret);
	return ret;
}

MatchHash::Result *MatchHash::oneMismatchFind(const char *s, unsigned len) const{
	Result *ret = 0;
	unsigned long hashVal = Hash::calcHashValue(s, s + len, &specialDnaToInt);
	for (unsigned i = 0, movLen = 0; i < len; i ++, movLen += 2){
		hashVal -= (specialDnaToInt(s[i]) << movLen);
		for (int j = 0;; j ++){
			if (dnaString[j] != s[i]) find_p(hashVal, ret);
			if (j < 3) hashVal += (1ULL << movLen);
			else break;
		}
		hashVal -= (3ULL - specialDnaToInt(s[i]) << movLen);
	}
	return ret;
}

void MatchHash::insert(Dna *seq, int start, int end){
	HashElement *element = new HashElement(seq, start, end);
	if (!insert_p(element)) delete element;
}

void MatchHash::insert(Dna *seq, int start, int end, unsigned long hashValue){
	HashElement *element = new HashElement(seq, start, end, hashValue);
	if (!insert_p(element)) delete element;
}

void MatchHash::remove(Dna *seq, int start, int end){
	unsigned hashVal = Hash::calcHashValue(seq -> dna() + start, seq -> dna() + end, &specialDnaToInt);
	remove_p(seq, start, hashVal);
}

void MatchHash::remove(Dna *seq, int start, int, unsigned long hashValue){
	remove_p(seq, start, hashValue);
}

MatchHash::HashElement::HashElement(Dna *seq, int start, int end) : next(0), seq(seq -> id()), start(start){
	hashValue = Hash::calcHashValue(seq -> dna() + start, seq -> dna() + end, &specialDnaToInt);
}

MatchHash::HashElement::HashElement(Dna *seq, int start, int end, unsigned long hashValue) : 
	hashValue(hashValue), next(0), seq(seq -> id()), start(start){
}

MatchHash::HashElement::~HashElement(){
}

bool MatchHash::HashElement::operator ==(const MatchHash::HashElement& element) const{
	return hashValue == element.hashValue;
}

/*
* BufferedMatchHash
*/
BufferedMatchHash::BufferedMatchHash(unsigned hashSize) : _hashSize(hashSize), _currentHashElement(0){
	_elements = new HashElement[_hashSize];
}

BufferedMatchHash::~BufferedMatchHash(){
	delete[] _elements;
}

Hash::Result *BufferedMatchHash::exactFind(const char *s, unsigned len) const{
	Result *ret = 0;
	unsigned long hashVal = Hash::calcHashValue(s, s + len, &specialDnaToInt);
	find_p(hashVal, ret);
	return ret;
}

Hash::Result *BufferedMatchHash::oneMismatchFind(const char *s, unsigned len) const{
	Result *ret = 0;
	unsigned long hashVal = Hash::calcHashValue(s, s + len, &specialDnaToInt);
	for (unsigned i = 0, movLen = 0; i < len; i ++, movLen += 2){
		hashVal -= (specialDnaToInt(s[i]) << movLen);
		for (unsigned j = 0;; j ++){
			if (dnaString[j] != s[i]) find_p(hashVal, ret);
			if (j < 3) hashVal += (1ULL << movLen);
			else break;
		}
		hashVal -= (3ULL - specialDnaToInt(s[i]) << movLen);
	}
	return ret;
}

void BufferedMatchHash::insert(Dna *seq, int start, int end){
	_elements[_currentHashElement] = HashElement(seq, start, end);
	insert_p(_currentHashElement ++);
}

void BufferedMatchHash::insert(Dna *seq, int start, int end, unsigned long hashValue){
	_elements[_currentHashElement] = HashElement(seq, start, end, hashValue);
	insert_p(_currentHashElement ++);
}

void BufferedMatchHash::remove(Dna *seq, int start, int end){
}

void BufferedMatchHash::remove(Dna *seq, int start, int end, unsigned long hashValue){
}

BufferedMatchHash::HashElement::HashElement() : next(0), seq(0){
}

BufferedMatchHash::HashElement::HashElement(Dna *seq, int start, int end) : next(0), seq(seq -> id()), start(start){
	hashValue = Hash::calcHashValue(seq -> dna() + start, seq -> dna() + end, &specialDnaToInt);
}

BufferedMatchHash::HashElement::HashElement(Dna *seq, int start, int end, unsigned long hashValue) : 
	hashValue(hashValue), next(0), seq(seq -> id()), start(start){
}

BufferedMatchHash::HashElement::~HashElement(){
}

bool BufferedMatchHash::HashElement::operator ==(const BufferedMatchHash::HashElement& element) const{
	return hashValue == element.hashValue;
}


/*
* BinaryHash
*/
BinaryHash::BinaryHash(unsigned long binSize) : _binSize(binSize), _binAnd((1ULL << binSize) - 1), _size(1ULL << _binSize){
	_elements = new HashElement *[_size];
	memset(_elements, 0, sizeof(HashElement *) * _size);
}

BinaryHash::~BinaryHash(){
	delete[] _elements;
}

void BinaryHash::find_p(unsigned int hashVal, MatchHash::Result *&ret) const{
	unsigned long id = hashVal & _binAnd;
	for (HashElement *e = _elements[id]; e; e = e -> next){
		if (e -> hashValue == hashVal){
			Result *res = new Result(e -> seq, e -> start);
			res -> next = ret;
			ret = res;
		}
	}
}

bool BinaryHash::insert_p(MatchHash::HashElement* element){
	unsigned long id = element -> hashValue & _binAnd;
	
	for (HashElement *e = _elements[id]; e; e = e -> next)
		if (e -> seq == element -> seq && e -> start == element -> start) return 0;
	element -> next = _elements[id];
	_elements[id] = element;
	return 1;
}

void BinaryHash::remove_p(Dna* seq, int start, unsigned long hashVal){
	unsigned long id = hashVal & _binAnd;
	
	if (!_elements[id]) return;
	if (_elements[id] -> seq == seq -> id() && _elements[id] -> start == start){
		HashElement *tmp = _elements[id];
		_elements[id] = tmp -> next;
		delete tmp;
		return;
	}
	
	for (HashElement *e = _elements[id]; e -> next; e = e -> next){
		if (e -> seq == seq -> id() && e -> next -> start == start){
			HashElement *tmp = e -> next;
			e -> next = tmp -> next;
			delete tmp;
			break;
		}
	}
}

/*
* BufferedBinaryHash
*/
BufferedBinaryHash::BufferedBinaryHash(unsigned long hashSize, unsigned int binSize) : 
	BufferedMatchHash(hashSize), _binSize(binSize), _binAnd((1ULL << binSize) - 1), _size(1ULL << _binSize){
		_elementBases = new unsigned[_size];
		memset(_elementBases, 0, sizeof(unsigned) * _size);
}

BufferedBinaryHash::~BufferedBinaryHash(){
	delete[] _elements;
}

void BufferedBinaryHash::find_p(unsigned long hashVal, MatchHash::Result *&ret) const{
	unsigned long id = hashVal & _binAnd;
	for (int e = _elementBases[id]; e; e = _elements[e].next){
		if (_elements[e].hashValue == hashVal){
			Result *res = new Result(_elements[e].seq, _elements[e].start);
			res -> next = ret;
			ret = res;
		}
	}
}

void BufferedBinaryHash::insert_p(unsigned elementId){
	unsigned long id = _elements[elementId].hashValue & _binAnd;
	_elements[elementId].next = _elementBases[id];
	_elementBases[id] = elementId;
}

