#include "mass_parts.h"
#include <algorithm>
#include <iostream>

class MassRatioCounter {
private:
	std::vector < int > counters;

	std::unordered_map < int, Scan > &reference_scans;

	static const size_t SECTIONS_NUMBER = 21;
public:
	MassRatioCounter(std::unordered_map < int, Scan > &reference_map) :
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

void calc_mass_parts(std::vector < std::string > filenames) {
	const std::string output_filename_pref = "thermo_mass_ratio_distribution";
	for (std::string pref : filenames) {
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
