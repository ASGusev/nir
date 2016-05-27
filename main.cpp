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

class ArgParser {
private:
	std::vector < std::string > filenames;

	std::vector < int > tasks;

	const char *FILENAME_KEY = "-f";
	const char *LIST_KEY = "-l";
	const char *TASK_KEY = "-t";
public:
	ArgParser(int argc, char **argv) {
		int cur_arg = 1;
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
			cur_arg++;
		}
	}

	std::vector < std::string > get_filenames() {
		return filenames;
	}

	std::vector < int > get_tasks() {
		return tasks;
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
			case 1: {
				check_zero_mass(filenames);
				break;
			}
			case 2: {
				check_no_peaks_scans(filenames);
				break;
			}
			case 3: {
				check_inclusion(filenames);
				break;
			}
			case 4: {
				check_same_mass(filenames);
				break;
			}
			case 5: {
				check_scans_finding(filenames);
				break;
			}
			case 6: {
				calc_mass_parts(filenames);
				break;
			}
			case 7: {
				check_thermo_scans_for_ratio(filenames);
				break;
			}
		}
	}

	std::cin.clear();
	std::cin.sync();
	std::cin.get();
	return 0;
}
