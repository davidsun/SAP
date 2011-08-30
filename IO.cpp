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


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <signal.h>

#include "IO.h"
#include "String.h"

#define MIN_VALID_LEN 16

/*
* Static functions
*/
inline static unsigned hexCharToUnsigned(char c){
	if (c <= '9' && c >= '0') return c - '0';
	else return c - 'a' + 10;
}

inline static char unsignedToHexChar(unsigned x){
	if (x < 10) return x + '0';
	else return x - 10 + 'a';
}


/*
* Class FileReader
*/
IO::FileReader::FileReader(const char* fileName) : _fileName(fileName){
	_f = fopen(_fileName.c_str(), "r");
}

IO::FileReader::~FileReader(){
	if (_f) fclose(_f);
}

bool IO::FileReader::isOpen() const{
	return (_f != 0);
}

int IO::FileReader::nextChar(){
	return fgetc(_f);
}

int IO::FileReader::readLine(DynamicArray <char> &ret){
	int size = 0, c = '\n';
	for (; c == '\n'; c = fgetc(_f));
	if (c == EOF) return EOF;
	for (; c != '\n' && c != EOF; c = fgetc(_f)){
		if (size >= ret.size()) ret.expand();
		ret[size ++] = c;
	}
	if (size >= ret.size()) ret.expand();
	ret[size] = 0;
	return size;
}

std::pair <int, int> IO::FileReader::readExon(DynamicArray <char> &name, DynamicArray <char> &dna, DynamicArray <char> &quality){
}

/*
* Class IO::BufferedFileReader::Buffer
*/
IO::BufferedFileReader::Buffer::Buffer(FILE* f) :  _next(0){
	_size = fread(_data, 1, READ_BUFFER_SIZE, f);
}

IO::BufferedFileReader::Buffer::~Buffer(){
}

IO::BufferedFileReader::Buffer* IO::BufferedFileReader::Buffer::next() const{
	return _next;
}

void IO::BufferedFileReader::Buffer::setNext(IO::BufferedFileReader::Buffer* next){
	_next = next;
}

int IO::BufferedFileReader::Buffer::size() const{
	return _size;
}

char IO::BufferedFileReader::Buffer::operator[](int x) const{
	return _data[x];
}


/*
* Class IO::BufferedFileReader
*/
IO::BufferedFileReader* IO::BufferedFileReader::newBufferedFileReader(const char* fileName){
	BufferedFileReader *ret = new BufferedFileReader(fileName);
	if (!ret -> _isEOF){
		pthread_attr_t attr;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
		pthread_create(&ret -> _readThread, &attr, readingProcess, (void *)ret);
	}
	return ret;
}

IO::BufferedFileReader::~BufferedFileReader(){
	Buffer *next;
	for (Buffer *p = _currentBuffer; p; p = next){
		next = p -> next();
		delete p;
	}
	
	_isToQuit = 1;
	pthread_cond_signal(&_fileReadPCondition);
	
	pthread_cond_destroy(&_fileReadPCondition);
	pthread_mutex_destroy(&_fileReadPMutex);
	pthread_mutex_destroy(&_bufferReadPMutex);
}

int IO::BufferedFileReader::readLine(DynamicArray <char> &ret){
	int size = 0;
	pthread_mutex_lock(&_bufferReadPMutex);
	
	int c = nextChar_p();
	if (c == EOF){
		pthread_mutex_unlock(&_bufferReadPMutex);
		return EOF;
	}
	for (; c != '\n' && c != EOF; c = nextChar_p()){
		if (size >= ret.size()) ret.expand();
		ret[size ++] = c;
	}
	pthread_mutex_unlock(&_bufferReadPMutex);
	if (size >= ret.size()) ret.expand();
	ret[size] = 0;
	return size;
}

std::pair <int, int> IO::BufferedFileReader::readExon(DynamicArray <char> &name, DynamicArray <char> &dna, DynamicArray <char> &quality){
	std::pair <int, int> ret(0, 0);
	pthread_mutex_lock(&_bufferReadPMutex);
	
	int c;
	for (c = '\n'; c == '\n'; c = nextChar_p());
	if (c == EOF){
		pthread_mutex_unlock(&_bufferReadPMutex);
		return std::make_pair(EOF, EOF);
	}
	for (; c != '\n' && c != EOF; c = nextChar_p()){
		if (ret.first >= name.size()) name.expand();
		name[ret.first ++] = c;
	}
	
	for (c = '\n'; c == '\n'; c = nextChar_p());
	if (c == EOF){
		pthread_mutex_unlock(&_bufferReadPMutex);
		return std::make_pair(EOF, EOF);
	}
	ret.second = 0;
	for (; c != '\n' && c != EOF; c = nextChar_p()){
		if (ret.second >= dna.size()) dna.expand();
		dna[ret.second ++] = c;
	}
	
	for (c = '\n'; c == '\n'; c = nextChar_p());
	if (c == EOF){
		pthread_mutex_unlock(&_bufferReadPMutex);
		return std::make_pair(EOF, EOF);
	}
	int qLen = 0;
	for (; c != '\n' && c != EOF; c = nextChar_p()){
		if (qLen >= quality.size()) quality.expand();
		quality[qLen ++] = c;
	}
	pthread_mutex_unlock(&_bufferReadPMutex);
	if (ret.first >= name.size()) name.expand();
	name[ret.first] = 0;
	if (ret.second >= dna.size()) dna.expand();
	dna[ret.second] = 0;
	if (qLen >= quality.size()) quality.expand();
	quality[qLen] = 0;
	String::dnaFormat(dna);
	return ret;
}

IO::BufferedFileReader::BufferedFileReader(const char *fileName) : 
	FileReader(fileName), _currentBufferLoc(0), _bufferCount(0), _isEOF(0), _isToQuit(0){
		pthread_cond_init(&_fileReadPCondition, NULL);
		pthread_mutex_init(&_fileReadPMutex, NULL);
		pthread_mutex_init(&_bufferReadPMutex, NULL);
		
		if (isOpen()){
			setvbuf(_f , 0, _IONBF, 0);
			
			_currentBuffer = _lastBuffer = new Buffer(_f);
			_bufferCount ++;
			if (_lastBuffer -> size() < READ_BUFFER_SIZE) _isEOF = 1;
		}
}

void *IO::BufferedFileReader::readingProcess(void *arg){
	BufferedFileReader *reader = (BufferedFileReader *)arg;
	pthread_mutex_lock(&reader -> _fileReadPMutex);
	while (1){
		while (reader -> _bufferCount >= READ_BUFFER_COUNT){
			pthread_cond_wait(&reader -> _fileReadPCondition, &reader -> _fileReadPMutex);
			if (reader -> _isToQuit){
				pthread_exit((void *)0);
				return 0;
			}
		}
		Buffer *nBuffer = new Buffer(reader -> _f);
		reader -> _lastBuffer -> setNext(nBuffer);
		reader -> _lastBuffer = nBuffer;
		reader -> _bufferCount ++;
		if (nBuffer -> size() < READ_BUFFER_SIZE){
			reader -> _isEOF = 1;
			break;
		}
	}
	pthread_mutex_unlock(&reader -> _fileReadPMutex);
	pthread_exit((void *)0);
}

void IO::BufferedFileReader::nextBuffer_p(){
	pthread_mutex_lock(&_fileReadPMutex);
	while (!_isEOF && _currentBuffer -> next() == 0){
		pthread_mutex_unlock(&_fileReadPMutex);
		usleep(30);
		pthread_mutex_lock(&_fileReadPMutex);
	}
	Buffer *last = _currentBuffer;
	_currentBufferLoc = 0;
	_currentBuffer = _currentBuffer -> next();
	
	delete last;
	_bufferCount --;
	pthread_cond_signal(&_fileReadPCondition);
	pthread_mutex_unlock(&_fileReadPMutex);
}

int IO::BufferedFileReader::nextChar_p(){
	while (_currentBuffer && _currentBufferLoc >= _currentBuffer -> size()) nextBuffer_p();
	if (!_currentBuffer) return EOF;
	return (*_currentBuffer)[_currentBufferLoc ++];
}


/*
 * Class FileWriter
 */
IO::FileWriter::FileWriter(const char* fileName){
	_f = fopen(fileName, "w");
}

IO::FileWriter::~FileWriter(){
	fclose(_f);
}

void IO::FileWriter::putString(const char* s, unsigned int size){
	fwrite(s, 1, size, _f);
}

void IO::FileWriter::putString(const DynamicArray <char> &s, unsigned size){
	putString(s.data(), size);
}

/*
 * Class BufferedFileWriter
 */
IO::BufferedFileWriter::Buffer::Buffer(const char *s, unsigned int size) : size(size), next(0){
	this -> s = new char[size];
	memcpy(this -> s, s, size);
}

IO::BufferedFileWriter::Buffer::~Buffer(){
	delete[] s;
}

IO::BufferedFileWriter *IO::BufferedFileWriter::newBufferedFileWriter(const char *fileName){
	BufferedFileWriter *ret = new BufferedFileWriter(fileName);
	
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	pthread_create(&ret -> _writingThread, &attr, writingProcess, (void *)ret);
	return ret;
}


IO::BufferedFileWriter::~BufferedFileWriter(){
	_isToQuit = 1;
	while (_currentBuffer){
		pthread_cond_signal(&_fileWriteCondition);
		usleep(30);
	}
	Buffer *next;
	for (Buffer *p = _currentBuffer; p; p = next){
		next = p -> next;
		delete p;
	}
	
	
	pthread_cond_destroy(&_fileWriteCondition);
	pthread_mutex_destroy(&_bufferWriteMutex);
}

void IO::BufferedFileWriter::putString(const DynamicArray <char> &s, unsigned int size){
	putString(s.data(), size);
}

void IO::BufferedFileWriter::putString(const char *s, unsigned int size){
    Buffer *buf = new Buffer(s, size);
	while (_bufferCount >= WRITE_BUFFER_COUNT) usleep(30);
	pthread_mutex_lock(&_bufferWriteMutex);
	if (_lastBuffer) _lastBuffer -> next = buf;
	else _currentBuffer = buf;
	_lastBuffer = buf;
	
	_bufferCount ++;
	pthread_mutex_unlock(&_bufferWriteMutex);
	pthread_cond_signal(&_fileWriteCondition);
}

IO::BufferedFileWriter::BufferedFileWriter(const char *fileName) : 
	FileWriter(fileName), _isToQuit(0), _bufferCount(0), _currentBuffer(0), _lastBuffer(0){
		setvbuf(_f , 0, _IONBF, 0);
		
		pthread_cond_init(&_fileWriteCondition, NULL);
		pthread_mutex_init(&_bufferWriteMutex, NULL);
}

void *IO::BufferedFileWriter::writingProcess(void *arg){
	BufferedFileWriter *writer = (BufferedFileWriter *)arg;
	pthread_mutex_lock(&writer -> _bufferWriteMutex);
	while (1){
		while (writer -> _bufferCount == 0 && !writer -> _isToQuit)
			pthread_cond_wait(&writer -> _fileWriteCondition, &writer -> _bufferWriteMutex);
		if (writer -> _bufferCount == 0 && writer -> _isToQuit){
			pthread_exit((void *)0);
			return 0;
		}
		bool unlock = 0;
		if (writer -> _currentBuffer != writer -> _lastBuffer) unlock = 1;
		if (unlock) pthread_mutex_unlock(&writer -> _bufferWriteMutex);
		fwrite(writer -> _currentBuffer -> s, 1, writer -> _currentBuffer -> size, writer -> _f);
		if (unlock) pthread_mutex_lock(&writer -> _bufferWriteMutex);
		writer -> nextBuffer_p();
	}
	pthread_mutex_unlock(&writer -> _bufferWriteMutex);
	pthread_exit((void *)0);
}

void IO::BufferedFileWriter::nextBuffer_p(){
	Buffer *last = _currentBuffer;
	_currentBuffer = last -> next;
	delete last;
	if (!_currentBuffer) _lastBuffer = 0;
	
	_bufferCount --;
}


/*
* Namespace IO
*/
namespace IO{
	unsigned long readUnsignedHex(FileReader &f){
		unsigned long ret = 0;
		int c = 0;
		for (; c != EOF && !((c <= '9' && c >= '0') || (c <= 'f' && c >= 'a')); c = f.nextChar());
		if (c == EOF) return 0;
		for (; (c <= '9' && c >= '0') || (c <= 'f' && c >= 'a'); c = f.nextChar())
			ret = (ret << 4U) + hexCharToUnsigned(c);
		return ret;
	}
	
	void writeUnsignedHex(FILE *f, unsigned long x){
		unsigned len = 0;
		char output[16];
		for (; x; x >>= 4) output[len ++] = (x & 15);
		if (len == 0) fputc('0', f);
		else {
			do {
				len --;
				fputc(unsignedToHexChar(output[len]), f);
			}  while (len);
		}
	}
	
	int readLine(FileReader &f, DynamicArray <char> &ret){
		int size = 0, c = '\n';
		for (; c == '\n'; c = f.nextChar());
		if (c == EOF) return -1;
		for (; c != '\n' && c != EOF; c = f.nextChar()){
			if (size >= ret.size()) ret.expand();
			ret[size ++] = c;
		}
		if (size >= ret.size()) ret.expand();
		ret[size] = 0;
		return size;
	}
	
	std::list <Dna *> readDna(const char *fileName){
		std::list <Dna *> ret;
		
		FileReader reader(fileName);
		if (!reader.isOpen()) return ret;
		DynamicArray <char> buffer(1000);
		int size;
		while (size = readLine(reader, buffer)){
			if (buffer[0] != '>'){
				int validLen = String::dnaFormat(buffer);
				if (validLen >= MIN_VALID_LEN) ret.push_back(new Dna(buffer.data(), size));
			}
		}
		return ret;
	}

	std::list <Exon *> readExon(const char *fileName){
		std::list <Exon *> ret;
		
		BufferedFileReader *reader = BufferedFileReader::newBufferedFileReader(fileName);
		if (!reader -> isOpen()) return ret;
		DynamicArray <char> bufferName(100), bufferDna(1000);
		int sizeName, sizeDna;
		while ((sizeName = reader -> readLine(bufferName) >= 0) && (sizeDna = reader -> readLine(bufferDna)) >= 0){
			String::toLower(bufferDna.data());
			ret.push_back(new Exon(bufferName.data(), bufferDna.data(), sizeDna));
		}
		delete reader;
		return ret;
	}
};
