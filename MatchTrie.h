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




#ifndef MATCHTRIE_H
#define MATCHTRIE_H

#include "MatchStructures.h"

#include <vector>

class MatchTrie{
	public:
		struct Result{
			char *start;
			Dna *seq;
			
			Result(const Result &res);
			Result(Dna *seq, char *start);
			
			bool operator < (const Result &res) const;
		};
		
		MatchTrie();
		~MatchTrie();
		
		std::vector <Result> find(const char *s, int len, int maxErrorCount = 0);
		void insert(Dna *seq, char *start, char *end);
		std::vector <Result> segmentFind(const char *s, int len, int segmentSize, int maxErrorCount = 0);
		void remove(Dna *seq, char *start, char *end);
		
	private:
		struct TrieNode{
			TrieNode();
			
			TrieNode *next[4];
			std::set <Result> data;
		};
		
		void eraseNode(TrieNode *&node);
		void find_p(TrieNode *n, int depth, const char *s, int len, int errorCount, std::vector <Result> &result);
		
		TrieNode *_root;
};

#endif
