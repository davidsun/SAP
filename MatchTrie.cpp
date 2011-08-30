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


#include "MatchTrie.h"
#include "MatchStructures.h"

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
* MatchTrie
*/
MatchTrie::MatchTrie(){
	_root = new TrieNode;
}

MatchTrie::~MatchTrie(){
	eraseNode(_root);
}

std::vector <MatchTrie::Result> MatchTrie::find(const char *s, int len, int maxErrorCount){
	std::vector <Result> ret;
	find_p(_root, 0, s, len, maxErrorCount, ret);
	return ret;
}

std::vector <MatchTrie::Result> MatchTrie::segmentFind(const char *s, int len, int segmentSize, int maxErrorCount){
	std::map <std::pair <Dna *, char *>, int> m;
	int count = 1;
	for (int i = 0; i < len; i += segmentSize, count ++){
		std::vector <Result> r;
		int startLoc = i;
		if (i + segmentSize > len){
			if (len < segmentSize) r = find(s + i, len, maxErrorCount);
			else r = find(s + len - segmentSize, segmentSize, maxErrorCount), startLoc = len - segmentSize;
		}  else r = find(s + i, segmentSize, maxErrorCount);
		for (int j = 0; j < r.size(); j ++)
			m[std::make_pair(r[j].seq, r[j].start - startLoc)] ++;
		for (std::map <std::pair <Dna *, char *>, int>::iterator it = m.begin(); it != m.end(); it ++)
			if (it -> second != count) m.erase(it);
	}
	std::vector <Result> ret;
	for (std::map <std::pair <Dna *, char *>, int>::iterator it = m.begin(); it != m.end(); it ++)
		ret.push_back(Result(it -> first.first, it -> first.second));
	return ret;
}


void MatchTrie::insert(Dna *seq, char* start, char* end){
	TrieNode *p = _root;
	for (char *c = start; c < end; c ++){
		int id = dnaToInt(*c);
		if (!(p -> next[id]))
			p -> next[id] = new TrieNode;
		p = p -> next[id];
		if (c - start >= MIN_INSERTING_LENGTH) p -> data.insert(Result(seq, start));
	}
}

void MatchTrie::eraseNode(MatchTrie::TrieNode *&node){
	if (!node) return;
	for (int i = 0; i < 4; i ++) eraseNode(node -> next[i]);
	delete node;
	node = 0;
}

void MatchTrie::remove(Dna* seq, char *start, char *end){
	TrieNode *p = _root;
	for (char *c = start; c < end; c ++){
		int id = dnaToInt(*c);
		if (!(p -> next[id])) return;
		p = p -> next[id];
		if (c - start >= MIN_INSERTING_LENGTH) p -> data.erase(Result(seq, start));
	}
}


void MatchTrie::find_p(MatchTrie::TrieNode *n, int depth, const char *s, int len, int errorCount, std::vector <MatchTrie::Result> &result){
	if (!n || errorCount < 0) return;
	if (depth == len){
		for (std::set <Result>::iterator it = n -> data.begin(); it != n -> data.end(); it ++)
			result.push_back(*it);
	}  else {
		for (int i = 0; i < 4; i ++)
			find_p(n -> next[i], depth + 1, s, len, errorCount - (dnaToInt(s[depth]) != i), result);
	}
}


MatchTrie::Result::Result(const MatchTrie::Result& res) : start(res.start), seq(res.seq){
}

MatchTrie::Result::Result(Dna* seq, char* start) : start(start), seq(seq){
}

bool MatchTrie::Result::operator <(const MatchTrie::Result &res) const{
	if (seq < res.seq) return 1;
	else if (seq > res.seq) return 0;
	else return start < res.start;
}

MatchTrie::TrieNode::TrieNode(){
	for (int i = 0; i < 4; i ++) next[i] = 0;
}
