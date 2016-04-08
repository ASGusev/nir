#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <cmath>

const double EPS = 1e-9;

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

const int GROUPS_OF_SCANS = 2;
const std::string THEORETIC_FILE_NAMES[GROUPS_OF_SCANS] = { "140509QXc1_car_anh_tryp_001.tsv", "140509QXc1_car_anh_tryp_004.tsv" };
const std::string MS_ALIGN_FILE_NAMES[GROUPS_OF_SCANS] = { "140509QXc1_car_anh_tryp_001_msdeconv.mgf", "140509QXc1_car_anh_tryp_004_msdeconv.mgf" };
const std::string THERMO_XTRACT_FILE_NAMES[GROUPS_OF_SCANS] = { "140509QXc1_car_anh_tryp_001_xtract.mgf", "140509QXc1_car_anh_tryp_004_xtract.mgf" };
const std::string MGF_SCAN_BEGINNING = "BEGIN IONS";
const std::string MGF_SCAN_ENDING = "END IONS";
const std::string ID_PREF = "TITLE=";
const std::string PEPMASS_PREF = "PEPMASS=";
const std::string CHARGE_PREF = "CHARGE=";
const int CHARGE_BORDER = 11;

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

const double NEUTRAL_MASS = 1.007276;

class NeutralCounter {
private:
	int counter;
public:
	NeutralCounter():counter(0) {};

	void operator() (Scan &scan) {
		if (abs(scan.mass - NEUTRAL_MASS) < EPS) {
			counter++;
		}
	}

	int get_number() {
		return counter;
	}
};

class ZeroPeaksCounter {
private:
	int counter;
	std::string program_name;
	std::unordered_map < int, Scan > &scan_map;
public:
	ZeroPeaksCounter(std::string program, std::unordered_map < int, Scan > &map):counter(0), program_name(program), scan_map(map) {};

	void operator() (Scan &scan) {
		if (scan.pikes_number == 0 && scan_map.find(scan.id) != scan_map.end()) {
			std::cout << program_name << ' ' << scan.id << ' ' << scan_map.find(scan.id)->second.e_value << std::endl;
			counter++;
		}
	}

	int get_number() {
		return counter;
	}
};

template < class Func >
void go_through_res(std::unordered_map < int, Scan > &theoretic_scans_map, std::string program, std::string filename, Func &action) {
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
	std::unordered_map < int, Scan > theoretic_scans_maps[GROUPS_OF_SCANS];
	for (int i = 0; i < GROUPS_OF_SCANS; i++) {
		std::vector < Scan > theoretic_scans;
		read_scans_from_tsv(THEORETIC_FILE_NAMES[i], theoretic_scans);
		for (unsigned int j = 0; j < theoretic_scans.size(); j++) {
			theoretic_scans_maps[i][theoretic_scans[j].id] = theoretic_scans[j];
		}
	}
	/*
	for (int i = 0; i < GROUPS_OF_SCANS; i++) {
		std::cerr << "Looking for scans with no peaks in part " << i + 1 << " of MS-Align+ results.\n";
		ZeroPeaksCounter ms_align_zero_peaks("MS_Align", theoretic_scans_maps[i]);
		go_through_res(theoretic_scans_maps[i], "MS_Align", MS_ALIGN_FILE_NAMES[i], ms_align_zero_peaks);
	}
	for (int i = 0; i < GROUPS_OF_SCANS; i++) {
		std::cout << "Looking for scans with no peaks in part " << i + 1 << " of Thermo Xtract results.\n";
		ZeroPeaksCounter thermo_zero_peaks("Thermo_Xtract", theoretic_scans_maps[i]);
		go_through_res(theoretic_scans_maps[i], "Thermo_Xtraxt", THERMO_XTRACT_FILE_NAMES[i], thermo_zero_peaks);
	}
	*/
	for (int i = 0; i < GROUPS_OF_SCANS; i++) {
		NeutralCounter thermo_counter;
		go_through_res(theoretic_scans_maps[i], "Thermo_Xtraxt", THERMO_XTRACT_FILE_NAMES[i], thermo_counter);
		std::cout << thermo_counter.get_number() << " neutral mass scans found in part " << i + 1 << " of Thermo Xtract output." << std::endl;
	}
	for (int i = 0; i < GROUPS_OF_SCANS; i++) {
		NeutralCounter ms_align_counter;
		go_through_res(theoretic_scans_maps[i], "MS_Align", MS_ALIGN_FILE_NAMES[i], ms_align_counter);
		std::cout << ms_align_counter.get_number() << " neutral mass scans found in part " << i + 1 << " of MS-ALign+ output." << std::endl;
	}
	
	std::cin.clear();
	std::cin.sync();
	std::cin.get();
	return 0;
}
