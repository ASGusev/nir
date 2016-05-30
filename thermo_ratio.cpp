#include "thermo_ratio.h"
#include <iostream>

class MassPartCounter {
private:
	const size_t MULTS_NUMBER;

	const std::vector < int > RATIOS;

	std::vector < int > counters;

	std::unordered_map < int, Scan > &reference_scans;

	const double EPS;
public:
	MassPartCounter(std::unordered_map < int, Scan > &new_map, const std::vector < int> &new_rats, double new_eps) :
		MULTS_NUMBER(new_rats.size()), RATIOS(new_rats), reference_scans(new_map), counters(MULTS_NUMBER), EPS(new_eps) {};

	void operator() (Scan &scan) {
		std::unordered_map < int, Scan >::iterator scan_in_ref = reference_scans.find(scan.id);
		if (scan_in_ref != reference_scans.end()) {
			Scan &reference_scan = scan_in_ref->second;
			for (int i = 0; i < MULTS_NUMBER; ++i) {
				if (abs(scan.mass * RATIOS[i] - reference_scan.mass) <= reference_scan.mass * EPS) {
					++counters[i];
				}
			}
		}
	}

	std::vector < int > get_res() {
		return counters;
	}
};

void check_thermo_scans_for_ratio(std::vector < std::string > filenames, double eps) {
	const std::vector < int > RATIOS = { 3, 4, 5 };
	for (std::string pref : filenames) {
		std::cout << "Looking for integer mass ratios in " << pref + XTRACT_SUF << std::endl;

		ScansMapCreator theoretic_collector;
		go_through_tsv(pref + TSV_SUF, theoretic_collector);
		std::unordered_map < int, Scan > theoretic_scans = theoretic_collector.get_map();

		MassPartCounter ratio_counter(theoretic_scans, RATIOS, eps);
		go_through_mgf(Thermo_Xtract, pref + XTRACT_SUF, ratio_counter);

		std::vector < int > quants = ratio_counter.get_res();
		for (int j = 0; j < RATIOS.size(); ++j) {
			std::cout << quants[j] << " scan(s) is/are " << RATIOS[j] << " times less that in theory.\n";
		}

		std::cout << std::endl;
	}
}
