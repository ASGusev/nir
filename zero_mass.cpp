#include "zero_mass.h"
#include <iostream>

void find_zero_mass(deconv_program program, std::string filename, double accuracy) {
	std::cout << "Looking for zero masses in file " << filename << std::endl;
	int counter = 0;
	go_through_mgf(program, filename, [&counter, accuracy](Scan &scan) {if (abs(scan.mass) < accuracy) ++counter;});
	std::cout << counter << " zero mass scans found." << std::endl;
}

void check_zero_mass(std::vector < std::string > filenames, double accuracy) {
	for (std::string pref : filenames) {
		ScansMapCreator map_creator;
		go_through_tsv(pref + TSV_SUF, map_creator);
		std::unordered_map < int, Scan > theoretic_scans_map = map_creator.get_map();

		find_zero_mass(MS_Align, pref + MSDECONV_SUF, accuracy);
		find_zero_mass(Thermo_Xtract, pref + XTRACT_SUF, accuracy);
	}
	std::cout << std::endl;
}
