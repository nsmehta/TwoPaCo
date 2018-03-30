#include "bin_conv.hpp"

int main(int argc, char **argv) {
	std::ifstream input;
	input.open(argv[1], std::ios::in);

	std::ofstream output;
	output.open(argv[2], std::ios::out);

	std::pair <bool, std::string> conv;
	while(isNext(input)) {
		conv = getNext(input);
		if (conv.first) output << "S ";
		else output << "P ";
		output << conv.second << '\n';
	}
	input.close();
	output.close();
}
