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



#ifndef MATCHSTRUCTURES_H
#define MATCHSTRUCTURES_H

#include <map>
#include <vector>
#include <set>
#include <list>
#include <pthread.h>

#include "DynamicArray.h"

#define MIN_INSERTING_LENGTH 9
#define MIN_SEGMENT_SIZE 3

const char dnaString[] = "atgc";

class Dna;
class DnaMatchInfo;
class DnaMatchInfoList;
class MatchTrie;

class Dna{
	public:
		Dna(const Dna &dna);
		Dna(char *dna, int size);
		~Dna();
		
		int id() const;
		char *dna() const;
		unsigned size() const;
		
		char &operator [](int x);
		char operator [](int x) const;
		
	protected:
		int _id;
		char *_dna;
		unsigned _size;
		
	private:
		static int totalDnaCount;
};

class Exon : public Dna{
	public:
		Exon(const Exon &exon);
		Exon(char *name, char *dna, int dnaSize);
		~Exon();
		
		char *name() const;
		
	private:
		char *_name;
};

class MatchExon : public Exon{
	public:
		class Insertion{
			friend class MatchExon;
			
			public:
				Insertion(char *dna, unsigned short size, double score);
				~Insertion();
				
				char *dna() const;
				Insertion *next() const;
				double score() const;
				
			private:
				double _score;
				unsigned short _size;
				Insertion *_next;
				char *_dna;
		};
		
		MatchExon(const Exon &exon);
		MatchExon(char *name, char *exon, int size);
		~MatchExon();
		
		int deleteCount(int loc) const;
		double deleteScore(int loc) const;
		void insert(int loc, char *s, int len, double score);
		Insertion *insertion(int loc);
		void lock(int loc, int len);
		int matchCount(int loc) const;
		int matchCount(int loc, char c) const;
		double matchScore(int loc) const;
		double matchScore(int loc, char c) const;
		char mostProbableDna(int loc) const;
		double totalQ(int loc) const;
		void unlock(int loc, int len);
		void updateDeletionValue(int loc, double score, double quality = 0.0);
		void updateMatchValue(int loc, char c, double score, double quality = 0.0);
		
	private:
		short *_deletionCount;			//Count of deletion happened in every location
		double *_deletionValue;			//Score of deletion happened in every location
		Insertion **_insertion;
		int *_lock;
		short *_matchCount[4];			//Count of matches
		double *_matchValue[4];			//Score of match for ATGC (the sum)
		double *_totalQ;
		
		void init();
};

class ExonList {
	public:
		class iterator{
			friend class ExonList;
			
			public:
				void operator ++();
				void operator ++(int);
				bool isEnd();
				Exon *exon() const;
			
			private:
				iterator(ExonList *parent, const std::map <int, Exon *>::iterator &it);
				
				ExonList *_parent;
				std::map <int, Exon *>::iterator _it;
		};
		
		ExonList();
		ExonList(const std::list <Exon *> &list);
		~ExonList();
		
		iterator begin();
		iterator exonById(int id);
		int readExon(const char *fileName);
		int readMatchExon(const char *fileName);
		iterator removeExon(const iterator &it); 
		int totalExonSize() const;
		
	private:
		std::map <int, Exon *> _infos;
		int _totalExonSize;
};

#endif
