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





#ifndef _HASH_H
#define _HASH_H

#include "MatchStructures.h"

#define HASH_NODE_SIZE 5
/*
* Defines the size of each node of class MatchHash
*/
#define MAX_HASH_SIZE 100000000
/*
* Defines the maximum size of MatchHash
*/

/*
* A virtual class
*/
class Hash{
	public:
		struct Result{
			int start;
			int seq;
			Result *next;
			
			Result(const Result &res);
			Result(int seq, int start);
		};
		
	public:
		Hash();
		~Hash();
		
		static void deleteResult(Result *&res);
		static unsigned long calcHashValue(const char *start, const char *end, unsigned long (*funcDnaToInt)(char c) = 0);
		
		virtual Result* exactFind(const char *s, unsigned len) const = 0;
		virtual Result* oneMismatchFind(const char *s, unsigned len) const = 0;
		virtual void insert(Dna *seq, int start, int end) = 0;
		virtual void insert(Dna *seq, int start, int end, unsigned long hashValue) = 0;
		virtual void remove(Dna *seq, int start, int end) = 0;
		virtual void remove(Dna *seq, int start, int end, unsigned long hashValue) = 0;
};

class MatchHash : public Hash{
	protected:
		struct HashElement{
			friend class MatchHash;
			
			unsigned long hashValue;
			HashElement *next;
			int seq;
			int start;
			
			HashElement(Dna *seq, int start, int end);
			HashElement(Dna *seq, int start, int end, unsigned long hashValue);
			~HashElement();
			
			bool operator == (const HashElement &element) const;
		};
		
	public:
		MatchHash();
		~MatchHash();
		
		Result* exactFind(const char *s, unsigned len) const;
		Result* oneMismatchFind(const char *s, unsigned len) const;
		void insert(Dna *seq, int start, int end);
		void insert(Dna *seq, int start, int end, unsigned long hashValue);
		virtual void remove(Dna *seq, int start, int end);
		virtual void remove(Dna *seq, int start, int end, unsigned long hashValue);
		
	private:
		virtual void find_p(unsigned hashVal, Result *&ret) const = 0;
		virtual bool insert_p(HashElement *element) = 0;
		virtual void remove_p(Dna *seq, int start, unsigned long hashValue) = 0;
};

class BufferedMatchHash : public Hash{
	protected:
		struct HashElement{
			friend class BufferedMatchHash;
			
			unsigned long hashValue;
			unsigned next;
			int seq;
			int start;
			
			HashElement();
			HashElement(Dna *seq, int start, int end);
			HashElement(Dna *seq, int start, int end, unsigned long hashValue);
			~HashElement();
			
			bool operator == (const HashElement &element) const;
		};
		
	public:
		BufferedMatchHash(unsigned hashSize);
		~BufferedMatchHash();
		
		Result* exactFind(const char *s, unsigned len) const;
		Result* oneMismatchFind(const char *s, unsigned len) const;
		void insert(Dna *seq, int start, int end);
		void insert(Dna *seq, int start, int end, unsigned long hashValue);
		void remove(Dna *seq, int start, int end);
		void remove(Dna *seq, int start, int end, unsigned long hashValue);
		
		
	protected:
		HashElement *_elements;
		
	private:
		unsigned _hashSize, _currentHashElement;
		
		virtual void find_p(unsigned long hashVal, Result *&ret) const = 0;
		virtual void insert_p(unsigned elementId) = 0;
};

class BinaryHash : public MatchHash{
	public:
		BinaryHash(unsigned long binSize);
		~BinaryHash();
		
	private:
		HashElement **_elements;
		
		unsigned long _binSize, _binAnd, _size;
		
		void find_p(unsigned hashVal, Result *&ret) const;
		bool insert_p(HashElement *element);
		void remove_p(Dna *seq, int start, unsigned long hashVal);
};

class BufferedBinaryHash : public BufferedMatchHash{
	public:
		BufferedBinaryHash(unsigned long hashSize, unsigned binSize);
		~BufferedBinaryHash();
		
	private:
		unsigned *_elementBases;
		
		unsigned long _binSize, _binAnd, _size;
		
		void find_p(unsigned long hashVal, Result *&ret) const;
		void insert_p(unsigned elementId);
};

#endif
