#include "bin_conv.hpp"

std::map<int64_t, int64_t> sequence_map;

std::map<char, std::string> mp_dna_to_int{
    {'A', "00"},
        {'C', "01"},
        {'G', "10"},
        {'T', "11"}
};

std::map<std::string, char> mp_int_to_dna{
    {"00", 'A'},
        {"01", 'C'},
        {"10", 'G'},
        {"11", 'T'}
};


std::string dna_to_int(char c) {
	if (mp_dna_to_int.find(c) != mp_dna_to_int.end()) {
		return mp_dna_to_int[c];
	}
	else {
		throw "invalid char convert fail";
	}
}

char int_to_dna(std::string s) {
	if (mp_int_to_dna.find(s) != mp_int_to_dna.end()) {
		return mp_int_to_dna[s];
	}
	else {
		throw "invalid string convert fail";
	}
}

int str_to_bin(std::string s) {
	int ret = 0;
	for (size_t i = 0; i < s.size(); i++) {
		if (s[i] == '1') ret |= (1<<i);
	}
	return ret;
}

std::string bin_to_str(int inp, int chars = 4) {
	std::string ret;
	for (int i = 0; i < 2*chars; i+=2) {
		std::string here;
		if (inp&(1<<i)) here += "1";
		else here += "0";
		if (inp&(1<<(i+1))) here += "1";
		else here += "0";
		ret += int_to_dna(here);
	}
	return ret;
}

unsigned char dna4_to_ascii(std::string seq) {
	std::string seq_str;
	for (size_t i = 0; i < seq.size(); i++) {
		seq_str += dna_to_int(seq[i]);
	}
	int seq_int = str_to_bin(seq_str);
	unsigned char ret = static_cast<unsigned char>(seq_int);
	return ret;
}


void sPack(int64_t segment_id, std::string seq, std::ofstream& output) {
  	int64_t id = sequence_map.size() + 1;
	sequence_map[segment_id] = id;
	std::pair<std::string, int> ret;
	ret.second = seq.size()%4;
	if (ret.second == 0) ret.second = 4;
	for (size_t i = 0; i < seq.size(); i+=4) {
		ret.first += dna4_to_ascii(seq.substr(i, 4));
	}
	std::string s1 = int_to_byte(ret.first.size());
	std::string s2 = std::to_string(ret.second);
	std::string s3 = ret.first;
	output << '0' << s1 << s2 << s3;

}

void pPack(std::vector<int64_t> p_v, std::ofstream& p_output) {
	int64_t max = 0;
	std::string ret;

	for (size_t i = 0; i < p_v.size(); i++) {
		bool sign = true;
		if (p_v[i] < 0) sign = false;
		p_v[i] = sequence_map[std::abs(p_v[i])];
		if (!sign) p_v[i] *= -1;
		if(max < abs(i)) {
			max = abs(p_v[i]);
		}
	}

	max = std::log2(max);
	max++; // 1 extra bit to store orientation

	max = (max + (8 - (max % 8))) / 8;

	int total_seq = p_v.size();

	std::string total_seq_str  = int_to_byte(total_seq);
	std::string max_str = int_to_byte(max);

	ret += '1' + total_seq_str + max_str;

	for (int64_t s: p_v) {

		if (s < 0) {
			s = (abs(s) << 1); // adding 1 bit for the orientation
		}
		else {
			s = (abs(s) << 1) | 1; // adding 1 bit for the orientation
		}

		std::string s1 = int_to_byte_generic((uint64_t)s, max);

		ret += s1;
	}
	p_output << ret;
}

std::string ascii_to_dna4(unsigned char inp, int chars = 4) {
	int inp_int = static_cast<int>(inp);
	std::string ret = bin_to_str(inp_int, chars);
	return ret;
}

bool isNext(std::ifstream& input) {
	if (input.is_open()) {
		return input.peek() != EOF;
	}
	std::cout << "File not open\n";
	return false;
}

std::pair <bool, std::string> getNext(std::ifstream& input) {
	if (isNext(input)) {
		std::pair <bool, std::string> ret;
		char c;
		input.get(c);
		// S flag
		if (c == '0') {
			ret.first = true;
			ret.second = sUnPack(input);
		} // P flag
		else if(c == '1') {
			ret.first = false;
			ret.second = pUnPack(input);
		}
		return ret;
	}
	else return std::make_pair(false, "");
}

void unPack(std::ifstream& input) {
	if (input.is_open()) {
		while(input.peek() != EOF) {
			char c;
			int s1, s2;
			std::string r, s3;
			input.get(c);
			// S flag
			if (c == '0') {
				sUnPack(input);
			} // P flag
			else if(c == '1') {
				pUnPack(input);
			}
		}
	}
}

std::string pUnPack(std::ifstream& input) {
	std::string ret;
	char c;
	std::string r;
	for (int i = 0; i < 4; i++) {
		input.get(c);
		r += c;
	}
	uint64_t no_of_seq = byte_to_int(r);
	std::string length_str;
	for (int i = 0; i < 4; i++) {
		input.get(c);
		length_str += c;
	}
	uint64_t seq_length = byte_to_int(length_str);

	for (uint64_t i = 0; i < no_of_seq; i++) {
		std::string temp;
		for (uint64_t j = 0; j < seq_length; j++) {
			input.get(c);
			temp += c;
		}
		uint64_t seq_with_rotation = byte_to_int_generic(temp, seq_length);
		char rotation = (seq_with_rotation & 1) ? '+' : '-';
		uint64_t sequence = seq_with_rotation >> 1;
		//std::cout << sequence << rotation;
		ret += std::to_string(sequence);
		ret += rotation;
		if (i < no_of_seq - 1) {
			//std::cout<<",";
			ret += ",";
		}
		else {
			//std::cout<< std::endl;
		}
	}
	if (DEBUG)
		std::cout << "P" << " " << ret << std::endl;
	return ret;
}

std::string sUnPack(std::ifstream& input) {
	std::string ret;
	char c;
	int s1, s2;
	std::string r, s3;

	for (int i = 0; i < 4; i++) {
		input.get(c);
		r += c;
	}
	s1 = byte_to_int(r);
	input.get(c);
	s2 = std::stoi(std::string(1, c));
	for (int i = 0; i < s1; i++) {
		input.get(c);
		s3 += c;
	}

	std::pair<std::string, int> inp = make_pair(s3, s2);
	std::string seq;
	for (size_t i = 0; i < inp.first.size()-1; i++) {
		seq += ascii_to_dna4(inp.first[i]);
	}
	seq += ascii_to_dna4(inp.first.back(), inp.second);
	ret = seq;
	if (DEBUG)
		std::cout << "S " << ret << '\n';
	return ret;
}

// converts 32-bit integer in its binary representation and then in n-character string
std::string int_to_byte_generic(uint64_t x, int byte_size) {
	std::vector<uint8_t> v;
	for (int i = 1; i <= byte_size; i++) {
		uint8_t b = (abs(x) & (((1<<(8*i)) - 1) << 8*(i - 1))) >> (8*(i -1));
		v.push_back(b);
	}

	std::string ret;
	for (uint8_t b_x: v) {
		ret += static_cast<unsigned char>(b_x);
	}
	return ret;
}

// converts 32-bit integer in its binary representation and then in 4-character string
std::string int_to_byte(int x) {
	uint8_t b0 = x & ((1<<8)-1);
	uint8_t b1 = (x & (((1<<8)-1) << 8)) >> 8;
	uint8_t b2 = (x & (((1<<8)-1) << 16)) >> 16;
	uint8_t b3 = (x & (((1<<8)-1) << 24)) >> 24;
	std::vector<uint8_t> v = {b0, b1, b2, b3};
	std::string ret;
	for (uint8_t b_x: v) {
		ret += static_cast<unsigned char>(b_x);
	}
	return ret;
}

uint64_t byte_to_int_generic(std::string s, int byte_size) {
	//std::cout << "Byte size : " << byte_size << '\n';
	uint64_t ret = 0;
	for (int i = 0; i < byte_size; i++) {
		uint64_t b = static_cast<uint8_t>(s[i]);
		ret = ret | (b << (8 * i));
		//std::cout << "val : " << (b << (8*i)) << '\t';
		//std::cout << ret << '\n';
	}
	return ret;
}

int byte_to_int(std::string s) {
	uint32_t b0 = static_cast<uint8_t>(s[0]);
	uint32_t b1 = static_cast<uint8_t>(s[1]);
	uint32_t b2 = static_cast<uint8_t>(s[2]);
	uint32_t b3 = static_cast<uint8_t>(s[3]);
	int ret = (b0 | (b1 << 8) | (b2 << 16) | (b3 << 24));
	return ret;
}

int not_main() {

	std::ofstream output;
	output.open("1.txt", std::ios::out);

	if (output.is_open()) {
		std::vector<std::string> v = {"ACGTAAAATGCATTTTAT", "A", "AC", "ACG", "ACGT", "ACGTT", "ACGTTG", "ACGTTGC", "ACGTTGCA", "ACGTTGCAA", "ACGTTGCAAC", "ACGTTGCAACG", "ACGTTGCAACGT"};
		std::vector<int64_t> cpv = {12, 11, 32, 28, 20, 16, 1, 2, 3, 4, 5, 6, 10};
		for (size_t i = 0; i < v.size(); i++) {
			sPack(cpv[i], v[i], output);
		}
	}

	std::ofstream p_output;

	if (output.is_open()) {
		std::vector<int64_t> p_v;
		p_v = {-12, 11, 32, -28, -20, 16};
		std::string out;
		pPack(p_v, output);

		p_v = {152, 131, 322, -283, -202, 1601, -233};
		//pPack(p_v, output);
		p_v = {1589662, 1345431, 32452, -2583, -2, 1233};
		//pPack(p_v, output);
	}

	if (output.is_open()) {
		std::vector<std::string> v = {"ACGTAAAATGCATTTTAT", "A", "AC", "ACG", "ACGT", "ACGTT", "ACGTTG", "ACGTTGC", "ACGTTGCA", "ACGTTGCAA", "ACGTTGCAAC", "ACGTTGCAACG", "ACGTTGCAACGT"};
		std::vector<int64_t> cpvv = {-12, 11, 32, -28, -20, 16, 1, 2, 3, 4, 5, 6, -10};		
		for (size_t i = 0; i < v.size(); i++) {
			sPack(cpvv[i]+100, v[i], output);
		}

	}

	output.close();


	std::ifstream input;
	input.open("1.txt", std::ios::in);

	unPack(input);

	input.close();
	return 0;
}
