/*    
Copyright (C) 2014 Gavin Wilson, Lincoln Stein

This file is part of RNASequel.

RNASequel is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

RNASequel is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with RNASequel.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _BINARY_IO_H
#define _BINARY_IO_H

#include <string>
#include <stdint.h>
#include <cstdlib>
#include <cstdio>

namespace rnasequel {

class BinaryWrite {
    public:
        BinaryWrite(const std::string &file) {
            open(file);
        }

        BinaryWrite() : _fp(NULL) { }

        bool open(const std::string & file) {
            _fp = fopen(file.c_str(),"wb");
            return _fp != NULL;
        }

	bool is_open() const {
	    return _fp != NULL;
	}

        bool operator! () const{
            return _fp == NULL;
        }

        ~BinaryWrite() {
            close();
        }

        void close(){
            if(_fp != NULL){
                fclose(_fp);
                _fp = NULL;
            }
        }

        template<typename T>
        void write(T v) {
            _write(reinterpret_cast<const char*>(&v), sizeof(T));
        }

        template<typename T>
        void write_vector(const T & v){
            write_vector(v, v.begin());
            /*
            write<size_t>(v.size());
            for(size_t i = 0; i < v.size(); i++){
                write<typename T::value_type>(v[i]);
            }
            */
        }

        template<typename T>
        void write_vector(const T & v, typename T::const_iterator start){
            size_t s = v.end() - start; 
            size_t i = start - v.begin();
            write<size_t>(s);
            const char * ptr = reinterpret_cast<const char *>(&v[i]);
            _write(ptr, sizeof(typename T::value_type) * s);
        }

        void write32(uint32_t v) {
            _write(reinterpret_cast<const char*>(&v), sizeof(uint32_t));
        }

        void write64(uint64_t v) {
            _write(reinterpret_cast<const char*>(&v), sizeof(uint64_t));
        }

        void write8(uint8_t v) {
            _write(reinterpret_cast<const char*>(&v), sizeof(uint8_t));
        }

        void write_size_t(size_t v) {
            _write(reinterpret_cast<const char*>(&v), sizeof(size_t));
        }

        void write_str(const std::string &v) {
            size_t l = v.length();
            write_size_t(l);
            _write(v.data(), l);
        }

        void write_n(const char * v, size_t l){
            _write(v, l);
        }

        long int tell() const {
            return ftell(_fp);
        }

    private:
        BinaryWrite(const BinaryWrite & w);
        BinaryWrite & operator=(const BinaryWrite & w);

        void _write(const char *v, size_t l) {
            fwrite(v, sizeof(char), l, _fp);
        }

        FILE *_fp;
};

class BinaryRead {
    public:
        BinaryRead(const std::string &file, long bsize = 0) {
            open(file, bsize);
        }

        BinaryRead() : _fp(NULL){ }

        bool open(const std::string & file, long bsize = 0) {
            _fp = fopen(file.c_str(),"rb");
            if(bsize > 0) {
                setvbuf( _fp , NULL , _IOFBF , bsize);
            }
            return _fp != NULL;
        }

        bool operator! () const{
            return _fp == NULL;
        }

	bool is_open() const {
	    return _fp != NULL;
	}

        ~BinaryRead() {
            close();
        }

        void close(){
            if(_fp != NULL){
                fclose(_fp);
                _fp = NULL;
            }
        }

        template<typename T>
        T read(){
            T v;
            _read(reinterpret_cast<char*>(&v), sizeof(T));
            return v;
        }

        template<typename T>
        void read(T & v){
            _read(reinterpret_cast<char*>(&v), sizeof(T));
        }

        template<typename T>
        void read_vector(T & v){
            /*
            size_t s = read<size_t>();
            v.resize(s);
            for(size_t i = 0; i < v.size(); i++){
                read<typename T::value_type>(v[i]);
            }*/

            size_t s = read<size_t>();
            v.resize(s);
            s *= sizeof(typename T::value_type);
            char * ptr = reinterpret_cast<char *>(&v[0]);
            _read(ptr, s);
        }

        uint32_t read32() {
	    uint32_t v = 0;
            _read(reinterpret_cast<char*>(&v), sizeof(uint32_t));
	    return v;
        }

        uint64_t read64() {
	    uint64_t v = 0;
            _read(reinterpret_cast<char*>(&v), sizeof(uint64_t));
	    return v;
        }

        size_t read_size_t() {
	    size_t v = 0;
            _read(reinterpret_cast<char*>(&v), sizeof(size_t));
	    return v;
        }

        void read_n(char *v, size_t n) {
            _read(v, n);
        }

        void read8(uint8_t &v) {
            _read(reinterpret_cast<char*>(&v), sizeof(uint8_t));
        }

        void read_str(std::string &v) {
            size_t l = read<size_t>();
            v.resize(l, 0);
            _read(&v[0], l);
        }

        void seek(long int offset) {
            fseek(_fp, offset, SEEK_SET);
        }

        long int tell() const {
            return ftell(_fp);
        }

    private:
        BinaryRead(const BinaryRead & w);
        BinaryRead & operator=(const BinaryRead & w);
        
        void _read(char *v, size_t l) {
            if(fread(v, sizeof(char), l, _fp) != l){
		fputs("Error reading file\n",stderr);
		exit(3);
	    }
        }

        FILE *_fp;
};

} // namespace bwt

#endif
