#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>

struct Scan {
	int id;
	double mass;
	int isotope_error;
	double precursor_error;
	int charge;
	std::string peptide;
	std::string protein;
	int de_novo_scope;
	int msgf_scope;
	double spec_e_value;
	double e_value;

	std::vector < std::pair < double, double > > peaks;

	double mass_2;
	
	void operator= (Scan &other);
};

enum deconv_program {MS_Align, Thermo_Xtract};
const std::string program_name[] = { "MS Align", "Thermo Xtract" };

const std::string MGF_SCAN_BEGINNING = "BEGIN IONS";
const std::string MGF_SCAN_ENDING = "END IONS";
const std::string ID_PREF = "TITLE=";
const std::string PEPMASS_PREF = "PEPMASS=";
const std::string CHARGE_PREF = "CHARGE=";
const double DA = 1.007276;
const int MAX_PEPTIDE_LENGTH = 70;
const double EPS = 1e-5;
const std::string TSV_SUF = ".tsv";
const std::string XTRACT_SUF = "_xtract.mgf";
const std::string MSDECONV_SUF = "_msdeconv.mgf";

void reset_scan(Scan &scan);

Scan parse_tsv_line(std::string &line);

int thermo_xtract_id(std::string);

int msalign_id(std::string);

template < class Func >
void go_through_mgf(deconv_program format, std::string filename, Func &action) {
	std::ifstream file(filename);
	Scan cur_scan;
	while (!file.eof()) {
		std::string line;
		getline(file, line);
		if (line == MGF_SCAN_BEGINNING) {
			reset_scan(cur_scan);
		}
		else if (line == MGF_SCAN_ENDING) {
			action(cur_scan);
		}
		else if (line.compare(0, ID_PREF.size(), ID_PREF) == 0) {
			if (format == MS_Align) {
				cur_scan.id = msalign_id(line);
			}
			else {
				cur_scan.id = thermo_xtract_id(line);
			}
		}
		else if (line.compare(0, PEPMASS_PREF.size(), PEPMASS_PREF) == 0) {
			cur_scan.mass = std::stod(line.substr(PEPMASS_PREF.size(), line.size() - PEPMASS_PREF.size())) - DA;
			if (format == Thermo_Xtract) {
				int space_pos = line.find(' ');
				if (space_pos != std::string::npos) {
					cur_scan.mass_2 = std::stod(line.substr(space_pos, line.size() - space_pos));
				}
			}
		}
		else if (line.compare(0, CHARGE_PREF.size(), CHARGE_PREF) == 0) {
			int charge_beginning = CHARGE_PREF.size();
			int charge_len = line.find_first_not_of("1234567890", charge_beginning) - charge_beginning;
			cur_scan.charge = std::stoi(line.substr(charge_beginning, charge_len));
		}
		else if (line[0] >= '0' && line[0] <= '9') {
			int space_pos = line.find(' ');
			std::pair < double, double > peak;
			peak.first = std::stod(line.substr(0, space_pos));
			peak.second = std::stod(line.substr(space_pos, line.size() - space_pos));
			cur_scan.peaks.push_back(peak);
		}
	}
	file.close();
}

template < class Func >
void go_through_tsv(std::string filename, Func &f) {
	std::ifstream file(filename);
	std::string line;
	getline(file, line);
	while (getline(file, line)) {
		Scan cur_scan = parse_tsv_line(line);
		f(cur_scan);
	}
	file.close();
}

class ScansMapCreator {
private:
	std::unordered_map < int, Scan > *scans_map;
public:
	ScansMapCreator();

	void operator() (Scan &scan);

	std::unordered_map < int, Scan > get_map();

	~ScansMapCreator();
};