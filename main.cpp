#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <algorithm>

const double EPS = 1e-5;

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

	int pikes_number;
};

void reset_scan(Scan &scan) {
	scan.id = 0;
	scan.mass = 0;
	scan.isotope_error = 0;
	scan.precursor_error = 0;
	scan.charge = 0;
	scan.spec_e_value = 0;
	scan.e_value = 0;
	scan.pikes_number = 0;
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
	//ans.mass *= ans.charge;
	for (int i = 0; i < 5; i++) {
		line_stream >> temp;
	}
	line_stream >> ans.e_value;

	return ans;
}

void read_scans_from_tsv(std::string filename, std::vector < Scan > &scans) {
	std::ifstream file(filename);
	std::string line;
	while (getline(file, line)) {
		if (line[0] != '#') {
			scans.push_back(parse_tsv_line(line));
		}
	}
	file.close();
}

const std::string THEORETIC_FILE_NAME = "140509QXc1_car_anh_tryp_001.tsv";
const std::string MS_ALIGN_FILE_NAME = "140509QXc1_car_anh_tryp_001_msdeconv.mgf";
const std::string THERMO_XTRACT_FILE_NAME = "140509QXc1_car_anh_tryp_001_xtract.mgf";
const std::string MGF_SCAN_BEGINNING = "BEGIN IONS";
const std::string MGF_SCAN_ENDING = "END IONS";
const std::string ID_PREF = "TITLE=";
const std::string PEPMASS_PREF = "PEPMASS=";
const std::string CHARGE_PREF = "CHARGE=";
const int CHARGE_BORDER = 11;

int get_int(std::string s) {
	std::stringstream stream(s);
	int ans;
	stream >> ans;
	return ans;
}

double get_double(std::string s) {
	std::stringstream stream(s);
	double ans;
	stream >> ans;
	return ans;
}

int msalign_id(std::string title) {
	return std::stoi(title.substr(ID_PREF.size() + 5, title.size() - (ID_PREF.size() + 5)));
}

int thermo_xtract_id(std::string title) {
	int eq_pos = title.find_last_of("=");
	return std::stoi(title.substr(eq_pos + 1, title.size() - eq_pos - 2));
}

bool check_scan(Scan &scan, std::unordered_map < int, Scan > &theoretic_scans_map) {
	std::unordered_map < int, Scan > ::iterator in_theoretic = theoretic_scans_map.find(scan.id);
	return (in_theoretic != theoretic_scans_map.end()) &&
		(scan.mass < EPS || scan.charge >= CHARGE_BORDER || in_theoretic->second.charge >= CHARGE_BORDER
			|| in_theoretic->second.mass == 0);
}

void go_through_res(std::unordered_map < int, Scan > &theoretic_scans_map, std::string program, std::string filename) {
	std::ifstream file(filename);
	Scan cur_scan;
	while (!file.eof()) {
		std::string line;
		getline(file, line);
		if (line == MGF_SCAN_BEGINNING) {
			reset_scan(cur_scan);
		}
		else if (line == MGF_SCAN_ENDING) {
			if (check_scan(cur_scan, theoretic_scans_map)) {
				std::unordered_map < int, Scan > ::iterator in_theoretic = theoretic_scans_map.find(cur_scan.id);
				std::cout << program << ' ' << cur_scan.id << ' ' << in_theoretic->second.e_value << ' ' << cur_scan.mass << ' ' << cur_scan.charge << std::endl;
			}
		}
		else if (line.compare(0, ID_PREF.size(), ID_PREF) == 0) {
			if (program[0] == 'M') {
				cur_scan.id = msalign_id(line);
			}
			else {
				cur_scan.id = thermo_xtract_id(line);
			}
		}
		else if (line.compare(0, PEPMASS_PREF.size(), PEPMASS_PREF) == 0) {
			cur_scan.mass = std::stod(line.substr(PEPMASS_PREF.size(), line.size() - PEPMASS_PREF.size()));
		}
		else if (line.compare(0, CHARGE_PREF.size(), CHARGE_PREF) == 0) {
			int charge_beginning = CHARGE_PREF.size();
			int charge_len = line.find_first_not_of("1234567890", charge_beginning) - charge_beginning;
			cur_scan.charge = std::stoi(line.substr(charge_beginning, charge_len));
		}
		else if (line[0] >= '0' && line[0] <= '9') {
			cur_scan.pikes_number++;
		}
	}
	file.close();
}

int main() {
	std::cout.sync_with_stdio(false);
	std::vector < Scan > theoretic_scans;
	read_scans_from_tsv(THEORETIC_FILE_NAME, theoretic_scans);
	std::unordered_map < int, Scan > theoretic_scans_map;
	for (unsigned int i = 0; i < theoretic_scans.size(); i++) {
		theoretic_scans_map[theoretic_scans[i].id] = theoretic_scans[i];
	}
	go_through_res(theoretic_scans_map, "Thermo_Xtraxt", THERMO_XTRACT_FILE_NAME);
	go_through_res(theoretic_scans_map, "MS_Align", MS_ALIGN_FILE_NAME);

	std::cin.clear();
	std::cin.sync();
	std::cin.get();
	return 0;
}
