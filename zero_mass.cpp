#include "zero_mass.h"
#include "scan.h"

void find_zero_mass(std::unordered_map < int, Scan > &theoretic_scans_map,
	deconv_program program, std::string filename) {
	std::cout << "Looking for zero masses in file " << filename << std::endl;
	int counter = 0;
	go_through_mgf(program, filename, [&counter](Scan &scan) {if (abs(scan.mass) < EPS) ++counter;});
	std::cout << counter << " zero mass scans found." << std::endl;
}

void check_zero_mass(std::vector < std::string > filenames) {
	for (std::string pref : filenames) {
		ScansMapCreator map_creator;
		go_through_tsv(pref + TSV_SUF, map_creator);
		std::unordered_map < int, Scan > theoretic_scans_map = map_creator.get_map();

		find_zero_mass(theoretic_scans_map, MS_Align, pref + MSDECONV_SUF);
		find_zero_mass(theoretic_scans_map, Thermo_Xtract, pref + XTRACT_SUF);
	}
	std::cout << std::endl;
}
