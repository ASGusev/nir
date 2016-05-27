#include "same_mass.h"
#include <unordered_map>
#include <iostream>

class MassTester {
private:
	std::unordered_map < int, Scan > &sample_scans;
	int same_counter;
	const double MASS_EPS = 1e-1;
public:
	MassTester(std::unordered_map < int, Scan > &ready_scans) :
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

void check_same_mass(std::vector < std::string > filenames) {
	for (std::string pref : filenames) {
		ScansMapCreator map_creator;
		go_through_tsv(pref + TSV_SUF, map_creator);
		std::unordered_map < int, Scan > theoretic_scans_map = map_creator.get_map();

		find_same_mass(theoretic_scans_map, pref + XTRACT_SUF, pref + MSDECONV_SUF);
	}
	std::cout << std::endl;
}
