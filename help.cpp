#include "help.h"
#include <string>

std::string HELP[] = {
	"This program is for getting some data about MS Deconv and Thermo Xtract deconvolution output.\n",
	"\n",
	"The command line arguements are:\n",
	"-f to add file name\n",
	"-l to load file names from a text file\n",
	"-e to specify good EValue border\n",
	"-a to specify mass accuracy\n",
	"-t to add task\n",
	"-h to view this page\n",
	"\n",
	"The following tasks are possible:\n",
	"1 Check for scans with zero mass.\n",
	"2 Check for scans without peaks.\n",
	"3 Check for scans with peaks found by MS Deconv only.\n",
	"4 Check for scans that have the same mass found by both programs.\n",
	"5 Check for scans by length that have correct mass found by each program.\n",
	"6 Find distribution of ratio between Thermo Xtract mass and theoretical mass.\n",
	"7 Check for scans with certain ratio between Thermo Xtract and theoretic mass.\n"
};

void show_help() {
	for (std::string line: HELP) {
		std::cout << line;
	}
}