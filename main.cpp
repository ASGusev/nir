#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include "scan.h"

const double EPS = 1e-5;

const int GROUPS_OF_SCANS = 2;
const std::string THEORETIC_FILE_NAMES[GROUPS_OF_SCANS] = { "140509QXc1_car_anh_tryp_001.tsv", "140509QXc1_car_anh_tryp_004.tsv" };
const std::string MS_ALIGN_FILE_NAMES[GROUPS_OF_SCANS] = { "140509QXc1_car_anh_tryp_001_msdeconv.mgf", "140509QXc1_car_anh_tryp_004_msdeconv.mgf" };
const std::string THERMO_XTRACT_FILE_NAMES[GROUPS_OF_SCANS] = { "140509QXc1_car_anh_tryp_001_xtract.mgf", "140509QXc1_car_anh_tryp_004_xtract.mgf" };
const double ERROR_BORDER = 1e-10;

void find_zero_mass(std::unordered_map < int, Scan > &theoretic_scans_map,
	deconv_program program, std::string filename) {
	std::cout << "Looking for zero masses in file " << filename << std::endl;
	int counter = 0;
	go_through_mgf(program, filename, [&counter](Scan &scan) {if (abs(scan.mass) < EPS) ++counter;});
	std::cout << counter << " zero mass scans found." << std::endl;
}

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

void check_zero_mass() {
	for (int i = 0; i < GROUPS_OF_SCANS; i++) {
		std::unordered_map < int, Scan > theoretic_scans_map;
		ScansMapCreator map_creator(theoretic_scans_map);
		go_through_tsv(THEORETIC_FILE_NAMES[i], map_creator);

		find_zero_mass(theoretic_scans_map, MS_Align, MS_ALIGN_FILE_NAMES[i]);
		find_zero_mass(theoretic_scans_map, Thermo_Xtract, THERMO_XTRACT_FILE_NAMES[i]);
	}
	std::cout << std::endl;
}

void check_no_peaks_scans() {
	for (int i = 0; i < GROUPS_OF_SCANS; i++) {
		std::unordered_map < int, Scan > theoretic_scans_map;
		ScansMapCreator map_creator(theoretic_scans_map);
		go_through_tsv(THEORETIC_FILE_NAMES[i], map_creator);

		find_no_peaks_scans(theoretic_scans_map, MS_Align, MS_ALIGN_FILE_NAMES[i]);
		find_no_peaks_scans(theoretic_scans_map, Thermo_Xtract, THERMO_XTRACT_FILE_NAMES[i]);
	}
	std::cout << std::endl;
}

void check_inclusion() {
	for (int i = 0; i < GROUPS_OF_SCANS; i++) {
		std::unordered_map < int, Scan > theoretic_scans_map;
		ScansMapCreator map_creator(theoretic_scans_map);
		go_through_tsv(THEORETIC_FILE_NAMES[i], map_creator);
		test_inclusion(theoretic_scans_map, THERMO_XTRACT_FILE_NAMES[i], MS_ALIGN_FILE_NAMES[i]);
	}
	std::cout << std::endl;
}

void check_same_mass() {
	std::unordered_map < int, Scan > theoretic_scans_maps[GROUPS_OF_SCANS];
	for (int i = 0; i < GROUPS_OF_SCANS; i++) {
		ScansMapCreator map_creator(theoretic_scans_maps[i]);
		go_through_tsv(THEORETIC_FILE_NAMES[i], map_creator);
	}

	for (int i = 0; i < GROUPS_OF_SCANS; i++) {
		find_same_mass(theoretic_scans_maps[i], THERMO_XTRACT_FILE_NAMES[i], MS_ALIGN_FILE_NAMES[i]);
	}
	std::cout << std::endl;
}

void make_theoretic_maps(std::unordered_map < int, Scan > *theoretic_scans_maps) {
	for (int i = 0; i < GROUPS_OF_SCANS; i++) {
		ScansMapCreator map_creator(theoretic_scans_maps[i]);
		go_through_tsv(THEORETIC_FILE_NAMES[i], map_creator);
	}
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

void check_lengh_distribution() {
	for (int i = 0; i < GROUPS_OF_SCANS; i++) {
		std::vector < int > dist = find_lengh_distribution(THEORETIC_FILE_NAMES[i]);
		write_vector(dist, "distribution" + std::to_string(i) + ".txt");
	}
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
	
	std::unordered_map < int, Scan > output;
	ScansMapCreator scans_map_creator(output);
	go_through_mgf(program, experimental_filename, scans_map_creator);

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

void check_scans_finding() {
	for (int i = 0; i < GROUPS_OF_SCANS; ++i) {
		std::cout << "Calculating peptide lengh distribution for " << THEORETIC_FILE_NAMES[i] << std::endl;
		std::vector < int > dist = find_lengh_distribution(THEORETIC_FILE_NAMES[i]);
		write_vector(dist, "distribution" + std::to_string(i) + ".txt");

		///check_mass_calculation(THEORETIC_FILE_NAMES[i], MS_ALIGN_FILE_NAMES[i], MS_Align, dist, "MS" + std::to_string(i) + ".txt");
		check_mass_calculation(THEORETIC_FILE_NAMES[i], THERMO_XTRACT_FILE_NAMES[i], Thermo_Xtract, dist, "thermo" + std::to_string(i) + ".txt");

		std::cout << THEORETIC_FILE_NAMES[i] << " is done." << std::endl;
	}
}

int main() {
	std::cout.sync_with_stdio(false);
	
	//check_zero_mass();
	//check_no_peaks_scans();
	//check_inclusion();
	//check_same_mass();
	check_scans_finding();

	std::cin.clear();
	std::cin.sync();
	std::cin.get();
	return 0;
}
