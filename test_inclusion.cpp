#include "test_inclusion.h"
#include <iostream>
#include <unordered_map>
#include <algorithm>

class ComparePeaks {
private:
	std::unordered_map < int, Scan > &thermo_scans;
public:
	ComparePeaks(std::unordered_map < int, Scan > &thermo_map) : thermo_scans(thermo_map) {}

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

void check_inclusion(std::vector < std::string > filenames) {
	for (std::string pref : filenames) {
		ScansMapCreator map_creator;
		go_through_tsv(pref + TSV_SUF, map_creator);
		std::unordered_map < int, Scan > theoretic_scans_map = map_creator.get_map();

		test_inclusion(theoretic_scans_map, pref + XTRACT_SUF, pref + MSDECONV_SUF);
	}
	std::cout << std::endl;
}
