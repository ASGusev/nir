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

const double ERROR_BORDER = 1e-10;

class ZeroPeaksCounter {
private:
	int all_no_peaks, good_no_peaks;
	deconv_program program;
	std::unordered_map < int, Scan > &scan_map;
public:
	ZeroPeaksCounter(deconv_program new_program, std::unordered_map < int, Scan > &map):
		all_no_peaks(0), good_no_peaks(0), program(new_program), scan_map(map) {};

	void operator() (Scan &scan) {
		std::unordered_map < int, Scan>::iterator prediction = scan_map.find(scan.id);
		if (scan.peaks.size() == 0 && prediction != scan_map.end()) {
			if (prediction->second.e_value <= ERROR_BORDER) {
				good_no_peaks++;
				std::cout << program << ' ' << scan.id << ' ' << prediction->second.e_value << std::endl;
			}
			all_no_peaks++;
		}
	}

	int get_all_no_peaks() {
		return all_no_peaks;
	}

	int get_dood_no_peaks() {
		return good_no_peaks;
	}
};

void find_no_peaks_scans(std::unordered_map < int, Scan > &theoretic_scans_map,
	deconv_program program, std::string filename) {
	std::cout << "Looking for scans without peaks in " << filename << std::endl;
	ZeroPeaksCounter zero_peaks(program, theoretic_scans_map);
	go_through_mgf(program, filename, zero_peaks);
	std::cout << "Found " << zero_peaks.get_all_no_peaks() << " scans without peaks including "
		<< zero_peaks.get_dood_no_peaks() << " good scans." << std::endl;
}

template < class Func >
class ScansCollector{
private:
	std::unordered_map < int, Scan > &good_scans;
	Func check_scan;
public:
	ScansCollector(std::unordered_map < int, Scan > &new_map, Func new_checker):
		check_scan(new_checker), good_scans(new_map) {}

	void operator() (Scan &scan) {
		if (check_scan(scan)) {
			good_scans[scan.id] = scan;
		}
	}

	std::unordered_map < int, Scan >& get_map() {
		return good_scans;
	}
};

class EValueTester {
private:
	const double E_VALUE_BORDER;
	std::unordered_map < int, Scan > &theoretic_map;
public:
	EValueTester (std::unordered_map < int, Scan > &new_map, double new_e_value_border):
		E_VALUE_BORDER(new_e_value_border), theoretic_map(new_map) {}

	bool operator()(Scan &scan) {
		std::unordered_map < int, Scan >::iterator prediction = theoretic_map.find(scan.id);
		return prediction != theoretic_map.end() && prediction->second.e_value < E_VALUE_BORDER;
	}
};

class ComparePeaks {
private:
	std::unordered_map < int, Scan > &thermo_scans;
public:
	ComparePeaks(std::unordered_map < int, Scan > &thermo_map): thermo_scans(thermo_map) {}

	void operator() (Scan &scan) {
		std::unordered_map < int, Scan >::iterator thermo_pos = thermo_scans.find(scan.id);
		if (thermo_pos != thermo_scans.end()) {
			Scan &thermo_scan = thermo_pos->second;
			std::sort(scan.peaks.begin(), scan.peaks.end());
			std::sort(thermo_scan.peaks.begin(), thermo_scan.peaks.end());
			int pos = 0;
			bool unique_found = false;
			for (int i = 0; i < scan.peaks.size() && !unique_found; i++) {
				while (pos < thermo_scan.peaks.size() && (thermo_scan.peaks[pos].first <
					scan.peaks[i].first - EPS || (abs(thermo_scan.peaks[pos].first -
						scan.peaks[i].first) < EPS && thermo_scan.peaks[pos].second <
						scan.peaks[i].second - EPS))) {
					++pos;
				}
				if (pos == thermo_scan.peaks.size() ||
					abs(thermo_scan.peaks[pos].first - scan.peaks[i].first) > EPS ||
					abs(thermo_scan.peaks[pos].second - scan.peaks[i].second) > EPS) {
					std::cout << "Scan " << scan.id << " from MS-Align+ contains pikes which are not found by Thermo Xtract.\n";
					unique_found = true;
				}
			}
		}
	}
};

void test_inclusion(std::unordered_map < int, Scan > theoretic_scans_map,
	std::string thermo_filename, std::string ms_align_file_name) {
	const double GOOD_BORDER = 1e-20;
	EValueTester scan_checker(theoretic_scans_map, GOOD_BORDER);
	std::unordered_map < int, Scan > good_thermo_scans;
	ScansCollector < EValueTester > thermo_storage(good_thermo_scans, scan_checker);
	go_through_mgf(Thermo_Xtract, thermo_filename, thermo_storage);
	
	ComparePeaks ms_align_tester(good_thermo_scans);
	go_through_mgf(MS_Align, ms_align_file_name, ms_align_tester);
}

class MassTester{
private:
	std::unordered_map < int, Scan > &sample_scans;
	int same_counter;
	const double MASS_EPS = 1e-1;
public:
	MassTester(std::unordered_map < int, Scan > &ready_scans):
		sample_scans(ready_scans), same_counter(0) {}

	void operator() (Scan &scan) {
		std::unordered_map < int, Scan >::iterator pos = sample_scans.find(scan.id);
		if (pos != sample_scans.end()) {
			Scan &same = pos->second;
			if (abs(scan.mass - same.mass) < MASS_EPS) {
				++same_counter;
			}
		}
	}

	int get_match_number() {
		return same_counter;
	}
};

void find_same_mass(std::unordered_map < int, Scan > theoretic_scans_map,
	std::string thermo_filename, std::string ms_align_file_name) {
	std::cout << "Looking for scans with the same mass in " << thermo_filename <<
		' ' << ms_align_file_name << std::endl;
	EValueTester scan_checker(theoretic_scans_map, ERROR_BORDER);
	std::unordered_map < int, Scan > good_thermo_scans;
	ScansCollector < EValueTester > thermo_scans(good_thermo_scans, scan_checker);
	go_through_mgf(Thermo_Xtract, thermo_filename, thermo_scans);

	MassTester ms_tester(good_thermo_scans);
	go_through_mgf(MS_Align, ms_align_file_name, ms_tester);
	std::cout << ms_tester.get_match_number() << " scans have the same mass in both files.\n";
}

void check_no_peaks_scans(std::vector < std::string > filenames) {
	for (std::string pref: filenames) {
		ScansMapCreator map_creator;
		go_through_tsv(pref + TSV_SUF, map_creator);
		std::unordered_map < int, Scan > theoretic_scans_map = map_creator.get_map();

		find_no_peaks_scans(theoretic_scans_map, MS_Align, pref + MSDECONV_SUF);
		find_no_peaks_scans(theoretic_scans_map, Thermo_Xtract, pref + XTRACT_SUF);
	}
	std::cout << std::endl;
}

void check_inclusion(std::vector < std::string > filenames) {
	for (std::string pref: filenames) {
		ScansMapCreator map_creator;
		go_through_tsv(pref + TSV_SUF, map_creator);
		std::unordered_map < int, Scan > theoretic_scans_map = map_creator.get_map();

		test_inclusion(theoretic_scans_map, pref + XTRACT_SUF, pref + MSDECONV_SUF);
	}
	std::cout << std::endl;
}

void check_same_mass(std::vector < std::string > filenames) {
	for (std::string pref: filenames) {
		ScansMapCreator map_creator;
		go_through_tsv(pref + TSV_SUF, map_creator);
		std::unordered_map < int, Scan > theoretic_scans_map = map_creator.get_map();

		find_same_mass(theoretic_scans_map, pref + XTRACT_SUF, pref + MSDECONV_SUF);
	}
	std::cout << std::endl;
}

class ScanDistributionCounter {
private:
	std::vector < int > &distribution;
public:
	ScanDistributionCounter(std::vector < int > &new_dist): distribution(new_dist) {}

	void operator()(Scan &scan) {
		++distribution[scan.peptide.size()];
	}
};

std::vector < int > find_lengh_distribution(std::string filename) {
	std::vector < int > distribution(MAX_PEPTIDE_LENGTH + 1);
	ScanDistributionCounter counter(distribution);
	go_through_tsv(filename, counter);
	return distribution;
}

template < typename T >
void write_vector(std::vector < T > &vect, std::string filename) {
	std::ofstream fout(filename);
	for (int i = 1; i < vect.size(); ++i) {
		fout << i << ' ' << vect[i] << std::endl;
	}
	fout.close();
}

class MassMatchChecker {
private:
	std::unordered_map < int, Scan > &experimental_scans;
	std::vector < int > &matching_masses;
	const double TOLERABLE_ERROR = 1e-5;
public:
	MassMatchChecker(std::unordered_map < int, Scan > &scans_map, std::vector < int > &new_masses):
		experimental_scans(scans_map), matching_masses(new_masses) {}

	void operator()(Scan &scan) {
		std::unordered_map < int, Scan >::iterator experimental_scan = experimental_scans.find(scan.id);
		if (experimental_scan != experimental_scans.end() && 
			abs(scan.mass - experimental_scan->second.mass) <= scan.mass * TOLERABLE_ERROR) {
				++matching_masses[scan.peptide.size()];
		}
	}
};

void check_mass_calculation(std::string theoretic_filename, std::string experimental_filename,
	deconv_program program, std::vector < int > lenght_distribution, std::string output_filename) {
	
	std::cout << "Counting correct masses in " << experimental_filename << std::endl;
	
	ScansMapCreator scans_map_creator;
	go_through_mgf(program, experimental_filename, scans_map_creator);
	std::unordered_map < int, Scan > output = scans_map_creator.get_map();

	std::vector < int > correct(MAX_PEPTIDE_LENGTH + 1);
	MassMatchChecker checker(output, correct);
	go_through_tsv(theoretic_filename, checker);
	
	std::vector < double > percent;
	for (int j = 0; j <= MAX_PEPTIDE_LENGTH; ++j) {
		if (lenght_distribution[j]) {
			percent.push_back(correct[j] * 100 / lenght_distribution[j]);
		}
		else {
			percent.push_back(0);
		}
	}
	write_vector(percent, output_filename);
}

void check_scans_finding(std::vector < std::string > filename) {
	for (std::string pref: filename) {
		std::cout << "Calculating peptide lengh distribution for " << pref + TSV_SUF << std::endl;
		std::vector < int > dist = find_lengh_distribution(pref + TSV_SUF);
		write_vector(dist, "distribution" + pref + ".txt");

		check_mass_calculation(pref + TSV_SUF, pref + MSDECONV_SUF, MS_Align, dist, "MS" + pref + ".txt");
		check_mass_calculation(pref + TSV_SUF, pref + XTRACT_SUF, Thermo_Xtract, dist, "thermo" + pref + ".txt");

		std::cout << pref + TSV_SUF << " is done." << std::endl;
	}
}

class MassRatioCounter {
private:
	std::vector < int > counters;

	std::unordered_map < int, Scan > &reference_scans;

	static const size_t SECTIONS_NUMBER = 21;
public:
	MassRatioCounter(std::unordered_map < int, Scan > &reference_map):
		counters(SECTIONS_NUMBER), reference_scans(reference_map) {}

	void operator ()(Scan &scan) {
		std::unordered_map < int, Scan >::iterator scan_in_ref = reference_scans.find(scan.id);
		if (scan_in_ref != reference_scans.end()) {
			Scan &reference_scan = scan_in_ref->second;
			double ratio = scan.mass / reference_scan.mass;
			int section = std::min(static_cast <int>(ratio * 10), 20);
			counters[section]++;
		}
	}

	std::vector < int > get_results() {
		return counters;
	}
};

class MassPartCounter {
private:
	const size_t MULTS_NUMBER;

	const std::vector < int > mults;

	std::vector < int > counters;

	std::vector < int > close_counters;

	std::unordered_map < int, Scan > &reference_scans;

	const double TOLERABLE_ERROR = 1e-5;

	const double CLOSE_TO_MATCH = 1e-4;
public:
	MassPartCounter(std::unordered_map < int, Scan > &new_map, const std::vector < int> &new_mults):
		MULTS_NUMBER(new_mults.size()), mults(new_mults), reference_scans(new_map), counters(MULTS_NUMBER),
		close_counters(MULTS_NUMBER) {};

	void operator() (Scan &scan) {
		std::unordered_map < int, Scan >::iterator scan_in_ref = reference_scans.find(scan.id);
		if (scan_in_ref != reference_scans.end()) {
			Scan &reference_scan = scan_in_ref->second;
			for (int i = 0; i < MULTS_NUMBER; ++i) {
				if (abs(scan.mass * mults[i] - reference_scan.mass) <= reference_scan.mass * TOLERABLE_ERROR) {
					std::cout << "Match:\n";
					++counters[i];
				}
				if (abs(scan.mass * mults[i] - reference_scan.mass) <= reference_scan.mass * CLOSE_TO_MATCH) {
					std::cout << "Scan " << scan.id << " has mass " << scan.mass << " and " <<
						reference_scan.mass << " in theory.\n";
					++close_counters[i];
				}
			}
		}
	}

	std::vector < int > get_res(){
		return counters;
	}

	std::vector < int > get_close_res() {
		return close_counters;
	}
};

void calc_mass_parts(std::vector < std::string > filenames) {
	const std::string output_filename_pref = "thermo_mass_ratio_distribution";
	for (std::string pref: filenames) {
		std::cout << "Calculating mass part distribution for " << pref + XTRACT_SUF << std::endl;

		ScansMapCreator theoretic_collector;
		go_through_tsv(pref + TSV_SUF, theoretic_collector);
		std::unordered_map < int, Scan > theoretic_scans = theoretic_collector.get_map();

		MassRatioCounter ratio_counter(theoretic_scans);
		go_through_mgf(Thermo_Xtract, pref + XTRACT_SUF, ratio_counter);
		std::vector < int > distribution = ratio_counter.get_results();
		
		std::ofstream fout(output_filename_pref + pref + ".txt");
		for (int j = 0; j < distribution.size(); ++j) {
			if (j == distribution.size() - 1) {
				fout << j * 10 << "%+ " << distribution[j] << std::endl;
			}
			else {
				fout << j * 10 << "%-" << (j + 1) * 10 << "% " << distribution[j] << std::endl;
			}
		}
		fout.close();

		std::cout << "Done." << std::endl;
	}
}

void check_thermo_scans_for_ratio(std::vector < std::string > filenames) {
	const std::vector < int > ratios = { 3, 4, 5 };
	for (std::string pref: filenames) {
		std::cout << "Looking for integer mass ratios in " << pref + XTRACT_SUF << std::endl;

		ScansMapCreator theoretic_collector;
		go_through_tsv(pref + TSV_SUF, theoretic_collector);
		std::unordered_map < int, Scan > theoretic_scans = theoretic_collector.get_map();

		MassPartCounter ratio_counter(theoretic_scans, ratios);
		go_through_mgf(Thermo_Xtract, pref + XTRACT_SUF, ratio_counter);
		
		std::vector < int > quants = ratio_counter.get_res();
		for (int j = 0; j < ratios.size(); ++j) {
			std::cout << quants[j] << " scan(s) is/are " << ratios[j] << " times less that in theory.\n";
		}

		std::vector < int > near_match = ratio_counter.get_close_res();
		for (int j = 0; j < ratios.size(); ++j) {
			std::cout << near_match[j] << " scan(s) is/are approximately " << ratios[j] << " times less that in theory.\n";
		}
		std::cout << std::endl;
	}
}

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
