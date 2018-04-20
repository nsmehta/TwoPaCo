#ifndef BIN_CONV_HPP
#define BIN_CONV_HPP

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>

#define DEBUG 1


extern std::map<int64_t, int64_t> sequence_map;

extern std::map<char, std::string> mp_dna_to_int;

extern std::map<std::string, char> mp_int_to_dna;

std::string dna_to_int(char c);

char int_to_dna(std::string s);

int str_to_bin(std::string s);

std::string bin_to_str(int inp, int chars);

unsigned char dna4_to_ascii(std::string seq);

void sPack(int64_t segment_id, std::string seq, std::ofstream& output);

void pPack(std::vector<int64_t> p_v, std::ofstream& p_output);

std::string ascii_to_dna4(unsigned char inp, int chars);

void unPack(std::ifstream& input);

bool isNext(std::ifstream& input);

std::pair<bool, std::string> getNext(std::ifstream& input);

std::string pUnPack(std::ifstream& input);

std::string sUnPack(std::ifstream& input);

std::string int_to_byte_generic(uint64_t x, int byte_size);

std::string int_to_byte(int x);

uint64_t byte_to_int_generic(std::string s, int byte_size);

int byte_to_int(std::string s);

#endif
