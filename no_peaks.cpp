#include "no_peaks.h"
#include <unordered_set>
#include <iostream>

class ZeroPeaksCounter {
private:
	int all_no_peaks, good_no_peaks;
	deconv_program program;
	std::unordered_map < int, Scan > &scan_map;
public:
	ZeroPeaksCounter(deconv_program new_program, std::unordered_map < int, Scan > &map) :
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

void check_no_peaks_scans(std::vector < std::string > filenames) {
	for (std::string pref : filenames) {
		ScansMapCreator map_creator;
		go_through_tsv(pref + TSV_SUF, map_creator);
		std::unordered_map < int, Scan > theoretic_scans_map = map_creator.get_map();

		find_no_peaks_scans(theoretic_scans_map, MS_Align, pref + MSDECONV_SUF);
		find_no_peaks_scans(theoretic_scans_map, Thermo_Xtract, pref + XTRACT_SUF);
	}
	std::cout << std::endl;
}