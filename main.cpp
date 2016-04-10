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
const int CHARGE_BORDER = 11;
const double ERROR_BORDER = 1e-9;

class NeutralCounter {
private:
	int counter;
	const double NEUTRAL_MASS = 1.007276;
public:
	NeutralCounter():counter(0) {};

	void operator() (Scan &scan) {
		if (abs(scan.mass - NEUTRAL_MASS) < EPS || abs(scan.mass_2 - NEUTRAL_MASS) < EPS) {
			counter++;
		}
	}

	int get_number() {
		return counter;
	}
};

void find_neutral_mass(std::unordered_map < int, Scan > &theoretic_scans_map,
	std::string program, std::string filename) {
	std::cout << "Looking for neutral masses in " << program << " file " << filename << std::endl;
	NeutralCounter counter;
	go_through_res(theoretic_scans_map, program, filename, counter);
	std::cout << counter.get_number() << " neutral mass scans found." << std::endl;
}

class ZeroPeaksCounter {
private:
	int all_no_peaks, good_no_peaks;
	std::string program_name;
	std::unordered_map < int, Scan > &scan_map;
public:
	ZeroPeaksCounter(std::string program, std::unordered_map < int, Scan > &map):
		all_no_peaks(0), good_no_peaks(0), program_name(program), scan_map(map) {};

	void operator() (Scan &scan) {
		std::unordered_map < int, Scan>::iterator prediction = scan_map.find(scan.id);
		if (scan.peaks.size() == 0 && prediction != scan_map.end()) {
			if (prediction->second.e_value <= ERROR_BORDER) {
				good_no_peaks++;
				std::cout << program_name << ' ' << scan.id << ' ' << prediction->second.e_value << std::endl;
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
	std::string program, std::string filename) {
	std::cerr << "Looking for " << program << " scans without peaks in " << filename << std::endl;
	ZeroPeaksCounter zero_peaks(program, theoretic_scans_map);
	go_through_res(theoretic_scans_map, program, filename, zero_peaks);
	std::cout << "Found " << zero_peaks.get_all_no_peaks() << " scans without peaks including "
		<< zero_peaks.get_dood_no_peaks() << " good scans.\n";
}

class ThermoScansCollector{
private:
	std::unordered_map < int, Scan > &theoretic_map;
	const double E_VALUE_BORDER = 1e-20;
	std::unordered_map < int, Scan > *good_thermo_scans;
public:
	ThermoScansCollector(std::unordered_map < int, Scan > &map): theoretic_map(map) {
		good_thermo_scans = new std::unordered_map < int, Scan >;
	}

	void operator() (Scan &scan) {
		std::unordered_map < int, Scan >::iterator prediction = theoretic_map.find(scan.id);
		if (prediction != theoretic_map.end() && prediction->second.e_value < E_VALUE_BORDER) {
			std::sort(scan.peaks.begin(), scan.peaks.end());
			(*good_thermo_scans)[scan.id] = scan;
		}
	}

	std::unordered_map < int, Scan >* get_map() {
		return good_thermo_scans;
	}
};

class MSAlignTester {
private:
	std::unordered_map < int, Scan > *thermo_scans;
public:
	MSAlignTester(std::unordered_map < int, Scan > *thermo_map) {
		thermo_scans = thermo_map;
	}

	void operator() (Scan &scan) {
		std::unordered_map < int, Scan >::iterator thermo_pos = thermo_scans->find(scan.id);
		if (thermo_pos != thermo_scans->end()) {
			Scan &thermo_scan = thermo_pos->second;
			std::sort(scan.peaks.begin(), scan.peaks.end());
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
	ThermoScansCollector thermo_storage(theoretic_scans_map);
	go_through_res(theoretic_scans_map, "Thermo_Xtract", thermo_filename, thermo_storage);
	std::unordered_map < int, Scan > *good_thermo_scans = thermo_storage.get_map();
	
	MSAlignTester ms_align_tester(good_thermo_scans);
	go_through_res(theoretic_scans_map, "MS-Align", ms_align_file_name, ms_align_tester);
}

int main() {
	std::cout.sync_with_stdio(false);
	std::unordered_map < int, Scan > theoretic_scans_maps[GROUPS_OF_SCANS];
	for (int i = 0; i < GROUPS_OF_SCANS; i++) {
		std::vector < Scan > theoretic_scans;
		read_scans_from_tsv(THEORETIC_FILE_NAMES[i], theoretic_scans);
		for (unsigned int j = 0; j < theoretic_scans.size(); j++) {
			theoretic_scans_maps[i][theoretic_scans[j].id] = theoretic_scans[j];
		}
	}
	
	for (int i = 0; i < GROUPS_OF_SCANS; i++) {
		find_neutral_mass(theoretic_scans_maps[i], "MS-Align", MS_ALIGN_FILE_NAMES[i]);
	}
	for (int i = 0; i < GROUPS_OF_SCANS; i++) {
		find_neutral_mass(theoretic_scans_maps[i], "Thermo Xtract", THERMO_XTRACT_FILE_NAMES[i]);
	}
	std::cout << std::endl;

	for (int i = 0; i < GROUPS_OF_SCANS; i++) {
		find_no_peaks_scans(theoretic_scans_maps[i], "MS-Align", MS_ALIGN_FILE_NAMES[i]);
	}
	for (int i = 0; i < GROUPS_OF_SCANS; i++) {
		find_no_peaks_scans(theoretic_scans_maps[i], "Thermo Xtract", THERMO_XTRACT_FILE_NAMES[i]);
	}
	std::cout << std::endl;
	
	for (int i = 0; i < GROUPS_OF_SCANS; i++) {
		test_inclusion(theoretic_scans_maps[i], THERMO_XTRACT_FILE_NAMES[i], MS_ALIGN_FILE_NAMES[i]);
	}
	
	std::cin.clear();
	std::cin.sync();
	std::cin.get();
	return 0;
}
