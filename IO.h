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




#ifndef IO_H
#define IO_H

#include <list>
#include <string>
#include <cstdio>

#define READ_BUFFER_SIZE 131072
#define READ_BUFFER_COUNT 64
#define WRITE_BUFFER_COUNT 16

#include "DynamicArray.h"
#include "MatchStructures.h"

namespace IO{
	class FileReader{
		public:
			FileReader(const char *fileName);
			~FileReader();
			
			bool isOpen() const;
			virtual inline int nextChar();
			virtual inline int readLine(DynamicArray <char> &ret);
			virtual inline std::pair <int, int> readExon(DynamicArray <char> &name, DynamicArray <char> &dna, DynamicArray <char> &quality);
			
		protected:
			FILE *_f;
			
		private:
			std::string _fileName;
	};
	
	class BufferedFileReader : public FileReader{
		public:
			~BufferedFileReader();
			
			static BufferedFileReader *newBufferedFileReader(const char *fileName);
			virtual inline int readLine(DynamicArray <char> &ret);
			virtual inline std::pair <int, int> readExon(DynamicArray <char> &name, DynamicArray <char> &dna, DynamicArray <char> &quality);
			
		private:
			class Buffer{
				public:
					Buffer(FILE *f);
					/*
					Buffer reads data immediately when the object is constructed
					*/
					~Buffer();
					
					inline int size() const;
					inline Buffer *next() const;
					void setNext(Buffer *next);
					inline char operator [](int x) const;
				
				private:
					unsigned _size;
					char _data[READ_BUFFER_SIZE];
					Buffer *_next;
			};
			
			BufferedFileReader(const char* fileName);
			
			bool _isEOF, _isToQuit;
			unsigned _bufferCount;
			unsigned _currentBufferLoc;
			Buffer *_currentBuffer, *_lastBuffer;
			
			pthread_t _readThread;
			pthread_mutex_t _fileReadPMutex;
			pthread_cond_t _fileReadPCondition;
			/*
				when a buffer is erased, a signal is send, to tell the reading thread to read in the new buffer
			*/
			pthread_mutex_t _bufferReadPMutex;
			
			static void *readingProcess(void *arg);
			void nextBuffer_p();
			int nextChar_p();
	};
	
	class FileWriter{
		public:
			FileWriter(const char *fileName);
			~FileWriter();
			
			virtual inline void putString(const DynamicArray <char> &s, unsigned size);
			virtual inline void putString(const char *s, unsigned size);
			
		protected:
			FILE *_f;
			
		private:
			std::string _fileName;
	};
	
	class BufferedFileWriter : public FileWriter{
		public:
			struct Buffer{
				unsigned size;
				char *s;
				Buffer *next;
				
				Buffer(const char *s, unsigned size);
				~Buffer();
			};
			
			~BufferedFileWriter();
			
			static BufferedFileWriter *newBufferedFileWriter(const char *fileName);
			
			virtual inline void putString(const DynamicArray <char> &s, unsigned size);
			virtual inline void putString(const char *s, unsigned size);
			
		private:
			pthread_t _writingThread;
			pthread_cond_t _fileWriteCondition;
			pthread_mutex_t _fileWriteMutex;
			pthread_mutex_t _bufferWriteMutex;
			
			bool _isToQuit;
			int _bufferCount;
			Buffer *_currentBuffer, *_lastBuffer;
			
			BufferedFileWriter(const char *fileName);
			
			static void *writingProcess(void *arg);
			void nextBuffer_p();
	};
	
	int readLine(FileReader &f, DynamicArray <char> &ret);
	
	std::list <Dna *> readDna(const char *fileName);
	
	std::list <Exon *> readExon(const char *fileName);
	
	unsigned long readUnsignedHex(const FileReader &f);
	
	void writeUnsignedHex(FILE *f, unsigned long x);
};

#endif
