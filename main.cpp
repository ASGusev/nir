#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <cstring>

#include "scan.h"
#include "zero_mass.h"
#include "no_peaks.h"
#include "test_inclusion.h"
#include "same_mass.h"
#include "scans_finding.h"
#include "mass_parts.h"
#include "thermo_ratio.h"
#include "help.h"

class ArgParser {
private:
	std::vector < std::string > filenames;

	std::vector < int > tasks;

	double evalue_border = 1e-10;
	double eps = 1e-5;

	const char *FILENAME_KEY = "-f";
	const char *LIST_KEY = "-l";
	const char *TASK_KEY = "-t";
	const char *ERROR_BORDER_KEY = "-e";
	const char *EPS_KEY = "-a";
	const char *HELP_KEY = "-h";
public:
	ArgParser(int argc, char **argv) {
		int cur_arg = 1;
		if (argc == 1) {
			tasks.push_back(0);
		}
		while (cur_arg < argc) {
			if (strcmp(FILENAME_KEY, argv[cur_arg]) == 0) {
				filenames.push_back(std::string(argv[++cur_arg]));
			}
			else if (strcmp(LIST_KEY, argv[cur_arg]) == 0) {
				std::string list_name(argv[++cur_arg]);
				std::ifstream file(list_name);
				std::string line;
				while (getline(file, line)) {
					filenames.push_back(line);
				}
				file.close();
			}
			else if (strcmp(TASK_KEY, argv[cur_arg]) == 0) {
				tasks.push_back(std::stoi(argv[++cur_arg]));
			}
			else if (strcmp(ERROR_BORDER_KEY, argv[cur_arg]) == 0) {
				evalue_border = std::stod(argv[++cur_arg]);
			}
			else if (strcmp(EPS_KEY, argv[cur_arg]) == 0) {
				eps = std::stod(argv[++cur_arg]);
			}
			else if (strcmp(HELP_KEY, argv[cur_arg]) == 0) {
				tasks.push_back(0);
			}
			cur_arg++;
		}
	}

	std::vector < std::string > get_filenames() {
		return filenames;
	}

	std::vector < int > get_tasks() {
		return tasks;
	}

	double get_evalue_border() {
		return evalue_border;
	}

	double get_eps() {
		return eps;
	}
};

int main(int argc, char **argv) {
	std::cout.sync_with_stdio(false);
	std::cout.precision(10);

	ArgParser args(argc, argv);
	std::vector < std::string > filenames = args.get_filenames();
	std::vector < int > tasks = args.get_tasks();
	
	for (int t: tasks) {
		switch (t) {
			case 0: {
				show_help();
			}
			case 1: {
				check_zero_mass(filenames, args.get_eps());
				break;
			}
			case 2: {
				check_no_peaks_scans(filenames, args.get_evalue_border());
				break;
			}
			case 3: {
				check_inclusion(filenames, args.get_evalue_border(), args.get_eps());
				break;
			}
			case 4: {
				check_same_mass(filenames, args.get_evalue_border(), args.get_eps());
				break;
			}
			case 5: {
				check_finding(filenames, args.get_eps());
				break;
			}
			case 6: {
				calc_mass_parts(filenames);
				break;
			}
			case 7: {
				check_thermo_scans_for_ratio(filenames, args.get_eps());
				break;
			}
		}
	}

	std::cin.clear();
	std::cin.sync();
	std::cin.get();
	return 0;
}
