#include "scans_finding.h"
#include <iostream>

class ScanDistributionCounter {
private:
	std::vector < int > distribution;
public:
	ScanDistributionCounter(): distribution(MAX_PEPTIDE_LENGTH + 1) {}

	void operator()(Scan &scan) {
		++distribution[scan.peptide.size()];
	}

	std::vector < int > get_distribution() {
		return distribution;
	}
};

class MassMatchChecker {
private:
	std::unordered_map < int, Scan > &experimental_scans;

	std::vector < int > matching_masses;
	
	const double EPS;
public:
	MassMatchChecker(std::unordered_map < int, Scan > &scans_map, double new_eps):
		experimental_scans(scans_map), matching_masses(MAX_PEPTIDE_LENGTH + 1), EPS(new_eps) {}

	void operator()(Scan &scan) {
		std::unordered_map < int, Scan >::iterator experimental_scan = experimental_scans.find(scan.id);
		if (experimental_scan != experimental_scans.end() &&
			abs(scan.mass - experimental_scan->second.mass) <= scan.mass * EPS) {
			++matching_masses[scan.peptide.size()];
		}
	}

	std::vector < int > get_results() {
		return matching_masses;
	}
};

std::vector < int > find_lengh_distribution(std::string filename) {
	ScanDistributionCounter counter;
	go_through_tsv(filename, counter);
	std::vector < int > distribution = counter.get_distribution();
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

void check_mass_calculation(std::string theoretic_filename, std::string experimental_filename,
	deconv_program program, std::vector < int > lenght_distribution, std::string output_filename, double eps) {

	std::cout << "Counting correct masses in " << experimental_filename << std::endl;

	ScansMapCreator scans_map_creator;
	go_through_mgf(program, experimental_filename, scans_map_creator);
	std::unordered_map < int, Scan > output = scans_map_creator.get_map();

	MassMatchChecker checker(output, eps);
	go_through_tsv(theoretic_filename, checker);
	std::vector < int > correct = checker.get_results();

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

void check_finding(std::vector < std::string > filename, double eps) {
	for (std::string pref: filename) {
		std::cout << "Calculating peptide lengh distribution for " << pref + TSV_SUF << std::endl;
		ScanDistributionCounter counter;
		go_through_tsv(pref + TSV_SUF, counter);
		std::vector < int > dist = counter.get_distribution();
		write_vector(dist, "distribution" + pref + ".txt");

		check_mass_calculation(pref + TSV_SUF, pref + MSDECONV_SUF, MS_Deconv, dist, "MS" + pref + ".txt", eps);
		check_mass_calculation(pref + TSV_SUF, pref + XTRACT_SUF, Thermo_Xtract, dist, "thermo" + pref + ".txt", eps);

		std::cout << pref + TSV_SUF << " is done." << std::endl;
	}
}
