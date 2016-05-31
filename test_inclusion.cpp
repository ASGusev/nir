#include "test_inclusion.h"
#include <iostream>
#include <unordered_map>
#include <algorithm>

class ComparePeaks {
private:
	const double EPS;

	std::unordered_map < int, Scan > &thermo_scans;
public:
	ComparePeaks(std::unordered_map < int, Scan > &thermo_map, double new_eps):
		thermo_scans(thermo_map), EPS(new_eps) {}

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
					std::cout << "Scan " << scan.id << " from MS-Deconv contains pikes which are not found by Thermo Xtract.\n";
					unique_found = true;
				}
			}
		}
	}
};

void check_inclusion(std::vector < std::string > filenames, double good_border, double eps) {
	for (std::string pref : filenames) {
		ScansMapCreator map_creator;
		go_through_tsv(pref + TSV_SUF, map_creator);
		std::unordered_map < int, Scan > theoretic_scans_map = map_creator.get_map();

		EValueTester scan_checker(theoretic_scans_map, good_border);
		std::unordered_map < int, Scan > good_thermo_scans;
		ScansCollector < EValueTester > thermo_storage(good_thermo_scans, scan_checker);
		go_through_mgf(Thermo_Xtract, pref + XTRACT_SUF, thermo_storage);

		ComparePeaks ms_deconv_tester(good_thermo_scans, eps);
		go_through_mgf(MS_Deconv, pref + MSDECONV_SUF, ms_deconv_tester);
	}
	std::cout << std::endl;
}
