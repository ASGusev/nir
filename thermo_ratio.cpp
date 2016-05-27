#include "thermo_ratio.h"
#include <iostream>

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
	MassPartCounter(std::unordered_map < int, Scan > &new_map, const std::vector < int> &new_mults) :
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

	std::vector < int > get_res() {
		return counters;
	}

	std::vector < int > get_close_res() {
		return close_counters;
	}
};

void check_thermo_scans_for_ratio(std::vector < std::string > filenames) {
	const std::vector < int > ratios = { 3, 4, 5 };
	for (std::string pref : filenames) {
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
