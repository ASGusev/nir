#include "scan.h"
#include <fstream>
#include <sstream>

void reset_scan(Scan &scan) {
	scan.id = 0;
	scan.mass = 0;
	scan.isotope_error = 0;
	scan.precursor_error = 0;
	scan.charge = 0;
	scan.spec_e_value = 0;
	scan.e_value = 0;
	scan.peaks.clear();
	scan.mass_2 = 0;
}

void Scan::operator=(Scan &other) {
	id = other.id;
	mass = other.mass;
	isotope_error = other.isotope_error;
	precursor_error = other.precursor_error;
	charge = other.charge;
	peptide = other.peptide;
	protein = other.protein;
	de_novo_scope = other.de_novo_scope;
	msgf_scope = other.msgf_scope;
	spec_e_value = other.spec_e_value;
	e_value = other.e_value;

	peaks = other.peaks;

	mass_2 = other.mass_2;
}

Scan parse_tsv_line(std::string &line) {
	const std::string ID_PREF = "scan=";

	Scan ans;
	std::stringstream line_stream(line);
	std::string temp;
	line_stream >> temp;
	line_stream >> temp;
	line_stream >> ans.id;
	line_stream >> temp;
	line_stream >> ans.mass;
	line_stream >> ans.isotope_error;
	line_stream >> ans.precursor_error;
	line_stream >> ans.charge;
	ans.mass = (ans.mass - DA) * ans.charge;
	line_stream >> ans.peptide;
	line_stream >> ans.protein;
	line_stream >> ans.de_novo_scope;
	line_stream >> ans.msgf_scope;
	line_stream >> ans.spec_e_value;
	line_stream >> ans.e_value;

	return ans;
}

int msalign_id(std::string title) {
	return std::stoi(title.substr(ID_PREF.size() + 5, title.size() - (ID_PREF.size() + 5)));
}

int thermo_xtract_id(std::string title) {
	int eq_pos = title.find_last_of("=");
	return std::stoi(title.substr(eq_pos + 1, title.size() - eq_pos - 2));
}

void ScansMapCreator::operator() (Scan &scan) {
	scans_map[scan.id] = scan;
}

std::unordered_map < int, Scan > ScansMapCreator::get_map() {
	return scans_map;
}
