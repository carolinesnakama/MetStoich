#include "stdafx.h"
#include "network.h"
#include <sstream>
#include <ctype.h>

using namespace std;
using namespace network;

/*
network.cpp
Purpose: reads a .dat or .txt file (similar to METATOOL) with a metabolic network.
It builds the stoichiometry matrix and classifies metabolites into external and internal, 
and reactions into reversible and irreversible.

@author Caroline S. M. Nakama
@version 1.0 04/17

*/

// Functions for Network

Network::Network(Eigen::MatrixXd S, bool intern){
	set_st_matrix(S, intern);
}

Network::Network(Eigen::MatrixXd S, bool intern, deque<bool> r){
	set_st_matrix(S, intern);
	try{
		set_rev_vector(r);
	}
	catch (length_error e){
		throw length_error(e.what());
	}
}

void Network::set_st_matrix(Eigen::MatrixXd S, bool intern){
	if (intern)
		int_st_matrix = S;
	else
		st_matrix = S;
}

void Network::set_rev_vector(deque<bool> r){
	if (!int_st_matrix.isZero() || !st_matrix.isZero()){
		if (int_st_matrix.cols() == r.size() || st_matrix.cols() == r.size())
			rev_vector = r;
		else
			throw length_error("Reversibility vector is not consistent with the stoichiometry matrix. Please check.");
	}
	else
		throw runtime_error("Stoichiometry matrix has not been defined.");
}

void Network::set_v_exp(Eigen::VectorXd v){
	v_exp = v;
}

Eigen::VectorXi Network::get_met_carb(){
	throw domain_error("This information is not available.");
}

void Network::print_metabolites(ofstream &file){
	throw domain_error("A detailed metabolic network was not provided. There is not a list of metabolites.");
}

void Network::print_reactions(ofstream &file){
	throw domain_error("A detailed metabolic network was not provided. A list of the reactions cannot be printed.");
}

void Network::clear_rev_vector(){
	rev_vector.clear();
}

void Network::check_st_matrix(){
	if (!st_matrix.isZero()){
		gen_ext_met_vector();
		if (ext_met_vector.empty()){
			throw "This stoichiometry matrix has no rows corresponding to external metabolites. Please check if the stoichiometry matrix provided is complete.";
		}
	} 
}

void Network::gen_ext_met_vector(){
	ext_met_vector.clear();
	for (size_t i = 0; i < st_matrix.rows(); i++){
		if (!(st_matrix.row(i).array() > 0.0).any() || !(st_matrix.row(i).array() < 0.0).any())
			ext_met_vector.push_back(i);
	}
}

void Network::gen_int_st_matrix(){
	if (int_st_matrix.isZero()){
		if (ext_met_vector.empty())
			gen_ext_met_vector();
		int_st_matrix.resize(st_matrix.rows() - ext_met_vector.size(), st_matrix.cols());
		ext_st_matrix.resize(ext_met_vector.size(), st_matrix.cols());
		size_t i = 0; size_t j = 0;
		while (i < st_matrix.rows()){
			if (i == ext_met_vector[j]){
				ext_st_matrix.row(j) = st_matrix.row(i);
				j++;
			}
			else
				int_st_matrix.row(i - j) = st_matrix.row(i);
			i++;
		}
	}
}

//void Network::gen_border_reac(){
//	border_reac.clear();
//	if (!int_st_matrix.isZero()){
//		for (size_t i = 0; i < int_st_matrix.cols(); i++){
//			if (!(int_st_matrix.col(i).array() > 0.0).any() || !(int_st_matrix.col(i).array() < 0.0).any())
//				border_reac.push_back(i);
//		}
//	}
//}

//std::vector<std::size_t> Network::get_border_reac(){
//	if (border_reac.empty())
//		gen_border_reac();
//	return border_reac;
//}

int Network::get_num_ext_met(){
	if (ext_met_vector.empty())
		gen_ext_met_vector();
	return ext_met_vector.size();
}

int Network::get_num_int_met(){
	if (int_st_matrix.isZero())
		gen_int_st_matrix();
	return int_st_matrix.rows();
}

int Network::get_num_total_reaction(){
	if (st_matrix.isZero())
		return int_st_matrix.cols();
	else
		return st_matrix.cols();
}

int Network::get_num_rev_reaction(){
	int i = 0;
	for (auto it = rev_vector.begin(); it != rev_vector.end(); it++){
		if (*it)
			i++;
	}
	return i;
}

int Network::get_num_irrev_reaction(){
	int i = 0;
	for (auto it = rev_vector.begin(); it != rev_vector.end(); it++){
		if (!*it)
			i++;
	}
	return i;
}

std::deque<bool> Network::get_rev_vector(){
	if (!rev_vector.empty())
		return rev_vector;
	else
		throw runtime_error("A vector indicating reaction reversibility has not been defined.");
}

Eigen::MatrixXd Network::get_st_matrix(){
	if (!st_matrix.isZero())
		return st_matrix;
	else
		throw length_error("A complete stoichiometry matrix was not defined.");
}

Eigen::MatrixXd Network::get_int_st_matrix(){
	if (int_st_matrix.isZero()){
		if (!st_matrix.isZero())
			gen_int_st_matrix();
		else
			throw length_error("No stoichiometry matrix was defined.");
	}
//	cout << int_st_matrix << endl;            //print para teste
	return int_st_matrix;
}

Eigen::MatrixXd Network::get_ext_st_matrix(){
	if (ext_st_matrix.isZero()){
		if (!st_matrix.isZero())
			gen_int_st_matrix();
		else
			throw length_error("No stoichiometry matrix was defined.");
	}
//	cout << ext_st_matrix << endl;           //print para teste
	return ext_st_matrix;
}

vector<int> Network::get_ids_upt_reactions(){        //TODO
	vector<int> ids;
	return ids;
}

vector<int> Network::get_ids_prod_reactions(){      //TODO
	vector<int> ids;
	return ids;
}

std::vector<std::string> Network::get_list_reactions(){
	vector<string> l;
	l.push_back("there are no names for reastions");
	return l;
}

void Network::clear(){
	delete this;
}

// Functions to Detailed Network

DetailedNetwork::DetailedNetwork(string arquivo){
	ext_met = 0;
	int_met = 0;
	km_met = 0;
	rev_reac = 0;
	mes_reac = 0;
	km_reac = 0;
	irrev_reac = 0;
	upt_reac = 0;
	exp_met = 0;
	redox = false;
	OxNum.resize(6);
	OxNum << 4, 1, -2, -3, 5, 6;      //C, H, O, N, P, S
	try{
		read_file(arquivo);
		gen_st_matrix();
		if (carbon)
			gen_carb_vector();
		gen_int_st_matrix();
		if (mes_reac > 0)
			gen_mes_st_matrix();
	}
	catch (invalid_argument e){
		cout << e.what();
	}
}

void DetailedNetwork::check_st_matrix(){
	ofstream log("log.txt");
	bool passed = true;
	bool ext_in = false;
	bool ext_out = false;
	bool rmet_int = false;
	vector<int> rrow;
	vector<int> in, out;
	log << "Internal metabolite verification" << endl;
	for (int i = 0; i < total_met; i++){
		in.push_back((st_matrix.row(i).array() > 0.0).count());
		out.push_back((st_matrix.row(i).array() < 0.0).count());
		if (!metabolites[i].is_ext()){
			//if (i >= ext_met){
			if (in[i] + out[i] > 1)
				log << metabolites[i].get_name() << " passed: Internal metabolite participates in at least two reactions" << endl;
			else {
				log << metabolites[i].get_name() << " failed: Internal metabolite must participate in ar least two reactions" << endl;
				passed = false;
			}
		}
		else{
			if (!ext_in)
				ext_in = in[i];
			if (!ext_out)
				ext_out = out[i];
		}
		if (in[i] == 0 && out[i] == 0){
			log << metabolites[i].get_name() << " does not participate in any reaction. Metabolite removed." << endl;
			rrow.push_back(i);
			if (metabolites[i].is_ext())
				ext_met -= 1;
			else{
				int_met -= 1;
				rmet_int = true;
			}
		}
	}
	log << endl << "Network overall verifcation" << endl;
	if (ext_in && ext_out)
		log << "Network passed: Network has inlet(s) and outlet(s)" << endl;
	else{
		log << "Network failed: Network must have at least one inlet e one outlet" << endl;
		passed = false;
	}
	if (!rrow.empty()){
		int size = rrow.size();
		for (int i = size - 1; i >= 0; i--){
			metabolites.erase(metabolites.begin() + rrow[i]);
			remove_row(st_matrix, rrow[i]);
			if (carbon)
				remove_el(met_carb, rrow[i]);
		}
		total_met = ext_met + int_met;
		if (rmet_int)
			gen_int_st_matrix();
	}
	log << endl << "Carbon balance verification" << endl;
	if (!check_carb_bal(log)){
		passed = false;
		log << "Carbon atoms are not balanced" << endl;
	}
	if (!passed)
		throw "There is a problem with this metabolic network. Please check.";
	if (exp)
		check_exp_bal(log);
}

void DetailedNetwork::check_exp_bal(ofstream &log){
	if (carbon){
		log << endl << "Carbon balance of experimental data: ";
		log << balance(true) << endl;
	}
	if (redox){
		log << endl << "Redox balance of experimental data: ";
		log << balance(false) << endl;
	}
}

double DetailedNetwork::balance(bool carb){
	int posmet;
	Eigen::VectorXd c = Eigen::VectorXd::Zero(mes_reac);
	double c_bal1 = 0.0; double c_bal2 = 0.0;
	for (int i = 0; i < mes_reac; i++){
		for (auto it = reactions[total_reac - mes_reac + i].get_metabolite_begin(); it != reactions[total_reac - mes_reac + i].get_metabolite_end(); it++){
			posmet = find_met_position(it->first);
			if (metabolites[posmet].is_exp()){
				if (carb)
					c(i) = carbon_balance(posmet, reactions[total_reac - mes_reac + i].get_flux());
				else
					c(i) = redox_balance(posmet, reactions[total_reac - mes_reac + i].get_flux());
			}
		}
		if (i < upt_reac)
			c_bal1 += c(i);
		else
			c_bal2 += c(i);
	}
	cout << c_bal1 << "  " << c_bal2 << endl;                         //print para teste
	return abs(c_bal2 / c_bal1);
}

double DetailedNetwork::carbon_balance(int pos, double flux){
	return metabolites[pos].get_carbon() * flux;
}

double DetailedNetwork::redox_balance(int pos, double flux){
	cout << OxNum.transpose() * metabolites[pos].get_composition() * flux << endl;
	return OxNum.transpose() * metabolites[pos].get_composition() * flux;
}

bool DetailedNetwork::check_carb_bal(ofstream &log){
	bool passed = true;
	if (carbon){               
		Eigen::VectorXd result;
		result = st_matrix.transpose() * met_carb.cast<double>();
		for (int i = 0; i < total_reac; i++){
			if (result(i) != 0){
				log << reactions[i].get_name() << " failed: Equation is not balanced with respect to C." << endl;
				passed = false;
			}
			else
				log << reactions[i].get_name() << " passed: Equation is balanced with respect to C." << endl;
		}
	}
	else{
		log << "Number of C atoms in each metabolite was not provided. Balance of carbon atoms was not performed." << endl;
	}
	return passed;
}

void DetailedNetwork::print_metabolites(ofstream &file){     //TODO: Adaptar para modelo cinetico
	for (int i = 0; i < total_met; i++)
		file << metabolites[i].get_name() << endl;
	if (kmodel){
		for (auto it = elst_enz_list.begin(); it != elst_enz_list.end(); it++){
			file << *it << endl;
		}
	}
}

void DetailedNetwork::print_reactions(ofstream &file){    //TODO: Consertar e adaptar para modelo cinetico
	//if (!kmodel){
		for (int i = 0; i < total_reac; i++){
			file << reactions[i].get_name() << " " << reactions[i].get_num_elem_steps() << endl; //<< ": ";
			//reactions[i].print_reaction(file);
			//if (i == 59){
			//	for (auto it = reactions[i].get_metabolite_begin(); it != reactions[i].get_metabolite_end(); it++){
			//		file << it->first << " " << it->second << endl;
			//	}
			//	file << endl;
			//}
		}
	//}
		file << endl << endl;
	//else{
		for (auto it = elst_reac_list.begin(); it != elst_reac_list.end(); it++){
			file << *it << endl;
		}
	//}
}

void DetailedNetwork::read_file(string arquivo){
	ifstream file(arquivo);
	string line, word;
	vector<string> vetor, met;
	vector<int> carb;
	int size = 0;
		
	if (file.is_open())
	{
		while (getline(file, line))
		{
			if (line.empty() || line[0] == '#' || line[0] == ' ')   // uma linha que contem informacao nao pode comecar com espaco
				continue;
			if (line[0] == '-')                                     //Nao pode ter espaco entre o titulo e o corpo
			{
				if (line.find("ENZREV") != string::npos)
				{
					kmodel = (line.find("[K]") != string::npos);
					getline(file, line);
					gen_vector(line, vetor);
					rev_reac = ((kmodel) ? vetor.size() / 2 : vetor.size());
					pop_enz_vector(vetor, true, false, kmodel);
				}	
				else if (line.find("ENZIRREV") != string::npos)
				{
					kmodel = (line.find("[K]") != string::npos);
					getline(file, line);
					gen_vector(line, vetor);
					irrev_reac = ((kmodel) ? vetor.size() / 2 : vetor.size());
					pop_enz_vector(vetor, false, false, kmodel);
				}
				else if (line.find("ENZMEAS") != string::npos)
				{
					bool k = (line.find("[K]") != string::npos);
					getline(file, line);
					gen_vector(line, vetor);
					mes_reac = ((k) ? vetor.size() / 2 : vetor.size());
					pop_enz_vector(vetor, false, true, k);
					getline(file, line);
					exp = !line.empty();
					if (exp){
						gen_vector(line, vetor);
						getline(file, line);
						// if (!line.empty()){
						// 	vector<string> vetor2;
						// 	gen_vector(line, vetor2);
						// 	class_enz_mes(vetor, vetor2);
						// }
						// else
						class_enz_mes(vetor);
					}
					irrev_reac += mes_reac;
				}
				else if (line.find("METINT") != string::npos)
				{
					carbon = (line.find("[C]") != string::npos);
					getline(file, line);
					gen_vector(line, vetor);
					get_carbon(vetor, met, carb, carbon);
					int_met = met.size() - ext_met;
					pop_met_vector(met, carb, false);
				}
				else if (line.find("METEXT") != string::npos)
				{
					carbon = (line.find("[C]") != string::npos);
					getline(file, line);
					gen_vector(line, vetor);
					get_carbon(vetor, met, carb, carbon);
					ext_met = met.size() - int_met;
					pop_met_vector(met, carb, true);
					getline(file, line);
					gen_vector(line, vetor);
					class_ext_met(vetor);
				}
				else if (line.find("CAT") != string::npos)                  //tem que vir depois das declaracoes das reacoes reversiveis, irreversiveis e medidas
				{    
					total_reac = rev_reac + irrev_reac;
					total_met = int_met + ext_met;
					vector<string> list;
					for (int i = 0; i < total_reac; i++){
						getline(file, line);
						gen_vector(line, list);
						add_reaction(list);
					}					
				}
				else if (line.find("COMP") != string::npos)
				{
					redox = true;
					for (int i = 0; i < exp_met; i++){
						getline(file, line);
						if (!line.empty()){
							gen_vector(line, vetor);
							add_met_composition(vetor);
						}
						else{
							throw length_error("Composition for all metabolites with experimentally measured consumption or production rates must be provided.");
						}
					}
				}
			}		
		}
	}
	else
	{
		throw invalid_argument("File could not be opened.");
	}
}

// void DetailedNetwork::class_enz_mes(vector<std::string> &vetor, vector<string> &vetor2){
void DetailedNetwork::class_enz_mes(vector<std::string> &vetor){
	int i = 0;
	// bool c = (vetor.size() == vetor2.size());
	for (auto it = vetor.begin(); it != vetor.end(); it++){
		//if ((*it)[0] == '-')
		//	upt_reac++;
		reactions[rev_reac + irrev_reac + i].add_flux(stod(*it)+0.0);
		// if (c)
		// 	reactions[rev_reac + irrev_reac + i].add_var(stod(vetor2[i]) + 0.0);
		i++;
	} 
	//cout << "upt_reac = " << upt_reac << endl;       //print para teste
}

void DetailedNetwork::class_ext_met(vector<std::string>& vetor){
	int i = 0;
	for (auto it = vetor.begin(); it != vetor.end(); it++){
		if (stoi(*it) == 1){
			exp_met++;
			metabolites[int_met + i].set_as_exp(true);
		}
		i++;
	}
	cout << "exp_met = " << exp_met << endl;       //print para teste
}

void DetailedNetwork::get_carbon(vector<string>& vetor, vector<string>& met, vector<int>& carb, bool carbon){
	if (carbon){
		size_t size = vetor.size();
		for (size_t i = 0; i < size; i = i + 2){
			met.push_back(vetor[i]);
			size_t pos = vetor[i + 1].find_first_of("]");
			carb.push_back(stoi(vetor[i + 1].substr(1, pos)));
		}
	}
	else
		met.insert(met.end(), vetor.begin(), vetor.end());
}

void DetailedNetwork::pop_enz_vector(vector<string>& vetor, bool rev, bool mes, bool k){
	size_t size = vetor.size();
	bool k2;
	for (size_t i = 0; i < size; i++){
		if (k){
			k2 = (vetor[i + 1].find("[1]") != string::npos);
			reactions.push_back(Reaction(vetor[i], rev, k2, mes));
			i++;
		}
		else
			reactions.push_back(Reaction(vetor[i], rev, false, mes));
	}
}

void DetailedNetwork::pop_met_vector(vector<string>& vetor, vector<int>& carb, bool ext){
	bool carbon = carb.empty();
	int sizei = ((ext) ? int_met : ext_met);
	int sizef = ext_met + int_met;
	for (int i = sizei; i < sizef; i++){
		metabolites.push_back(Metabolite(vetor[i], ext));
		if (!carbon)
			metabolites[i].add_carbon(carb[i]);
	}
}

void DetailedNetwork::gen_vector(string& line, vector<string>& list){
	list.clear();
	istringstream iss(line);
	for (string s; iss >> s;){
		list.push_back(s);
	//	cout << s << endl;                //print para teste
	}
}

int DetailedNetwork::find_reaction_position(string name){
	int i;
	for (i = 0; i < total_reac; i++){
		if (reactions[i].get_name() == name)
			break;
	}
	return i;
}

int DetailedNetwork::find_met_position(string name){
	int i;
	for (i = 0; i < total_met; i++){
		if (metabolites[i].get_name() == name)
			break;
	}
	return i;
}

void DetailedNetwork::add_reaction(vector<string> list){
	int posmet;
	double coef;

	int list_end = list.size() - 1;
	int posreac = find_reaction_position(list[0]);
	int poseq = find(list.begin(), list.end(), "=") - list.begin();

	for (int i = 2; i < list_end; i++){
		if (list[i] == "+" || list[i] == "=")
			continue;
		//if (isdigit(list[i][0]) || list[i][0] == '.'){
		//if (all_of(list[i].begin(), list[i].end(), ::isdigit) || list[i][0] == '.'){
		if (list[i].find_first_not_of("0123456789.") == string::npos){
			coef = stod(list[i]);
			++i;
		}
		else{
			coef = 1;
		}
		posmet = find_met_position(list[i]);
		if (i < poseq){
			reactions[posreac].add_metabolite(list[i], -coef, metabolites[posmet].get_carbon());
//			reactions[posreac].add_metabolite(posmet, -coef);
//			st_matrix(posmet, posreac) = -coef;
		}
		else{
//			st_matrix(posmet, posreac) = coef;
//			reactions[posreac].add_metabolite(posmet, coef);
			reactions[posreac].add_metabolite(list[i], coef, metabolites[posmet].get_carbon());
		}
	}
	if (kmodel)
		reactions[posreac].decomp_elem_steps();
}

void DetailedNetwork::add_met_composition(vector<string> &list){
	Eigen::VectorXi c = Eigen::VectorXi::Zero(6);
	map<string, int> order = { { "C", 0 }, { "H", 1 }, { "O", 2 }, { "N", 3 }, { "P", 4 }, { "S", 5 } };
	int posmet = find_met_position(list[0]);
	for (auto it = list.begin() + 2; it != list.end() - 1; it++){
		if (*it == "+")
			continue;
		if (it->find_first_not_of("0123456789") == string::npos){
			c(order[*(it + 1)]) = stoi(*it);
			it++;
		}
		else
			c(order[*(it)]) = 1;
	}
	//cout << c << endl << endl;                               // print para teste
	metabolites[posmet].add_composition(c); //TODO
}

void DetailedNetwork::gen_st_matrix(){
	int size, posmet;
	st_matrix.resize(total_met, total_reac);
	st_matrix << Eigen::MatrixXd::Zero(total_met, total_reac);
	for (int i = 0; i < total_reac; i++){
		for (auto it = reactions[i].get_metabolite_begin(); it != reactions[i].get_metabolite_end(); it++){
			posmet = find_met_position(it->first);
			st_matrix(posmet, i) = it->second;
		}
		//cout << reactions[i].get_name() << endl;                                     //print para teste
	}
}

void DetailedNetwork::gen_carb_vector(){
	met_carb.resize(total_met);
	for (int i = 0; i < total_met; i++){
		met_carb(i) = metabolites[i].get_carbon();
	}
}

void DetailedNetwork::gen_int_st_matrix(){
	if (metabolites[0].is_ext()){
		int_st_matrix = st_matrix.bottomRows(int_met);
		ext_st_matrix = st_matrix.topRows(ext_met);
	}
	else{
		int_st_matrix = st_matrix.topRows(int_met);
		ext_st_matrix = st_matrix.bottomRows(ext_met);
	}
}

void DetailedNetwork::gen_mes_st_matrix(){
	//vector<int> not_mes;
	mes_st_matrix = ext_st_matrix.bottomRows(exp_met);
	Eigen::MatrixXd m;
	m = mes_st_matrix.rightCols(mes_reac);
	for (size_t i = 0; i < mes_reac; i++){
		if ((m.col(i).array() < 0.0).any()){
			upt_reac++;
		}
	}
	//cout << "upt_reac = " << upt_reac << endl;
}

void DetailedNetwork::gen_elst_st_matrix(){
	calc_num_elst();
	elst_st_matrix = Eigen::MatrixXd::Zero(km_met, km_reac);
	elst_subs_st_matrix = Eigen::MatrixXd::Zero(km_met, km_reac);
	int count_met = total_met;
	int count_reac = 0;
	int reac_elst = 0;
	int posmet;
	vector<string> vetor;
	for (auto it = reactions.begin(); it != reactions.end(); it++){
		if (!it->elem_steps_decomp()){
			for (auto itt = it->get_metabolite_begin(); itt != it->get_metabolite_end(); itt++){
				posmet = find_met_position(itt->first);
				elst_st_matrix(posmet, count_reac) = itt->second;
			}
			count_reac++;
			elst_reac_list.push_back(it->get_name());
		}
		else{
			reac_elst = it->get_num_elem_steps();
			for (auto itt = it->get_elem_steps_enz_begin(true); itt != it->get_elem_steps_enz_end(true); itt++){
				elst_st_matrix(count_met + itt->first, count_reac + itt->second) = -1;
				elst_subs_st_matrix(count_met + itt->first, count_reac + itt->second) = -1;
			}
			for (auto itt = it->get_elem_steps_enz_begin(false); itt != it->get_elem_steps_enz_end(false); itt++){
				elst_st_matrix(count_met + itt->first, count_reac + itt->second) = 1;
			}
			for (auto itt = it->get_elem_steps_met_begin(true); itt != it->get_elem_steps_met_end(true); itt++){
				posmet = find_met_position(itt->first);
				elst_st_matrix(posmet, count_reac + itt->second) = -1;
				elst_subs_st_matrix(posmet, count_reac + itt->second) = -1;
			}
			for (auto itt = it->get_elem_steps_met_begin(false); itt != it->get_elem_steps_met_end(false); itt++){
				posmet = find_met_position(itt->first);
				elst_st_matrix(posmet, count_reac + itt->second) = 1;
			}
			count_met = count_met + reac_elst / 2;
			vetor = it->get_enz_list();
			elst_enz_list.insert(elst_enz_list.end(), vetor.begin(), vetor.end());
			count_reac = count_reac + reac_elst;
			vetor = it->get_elst_list();
			elst_reac_list.insert(elst_reac_list.end(), vetor.begin(), vetor.end());
		}
	}
}

void DetailedNetwork::gen_rev_st_matrix(){
	rev_st_matrix = Eigen::MatrixXd::Zero(total_met, total_reac + rev_reac);
	int j = 0;
	for (int i = 0; i < rev_reac; i++){
		rev_st_matrix.col(j) = st_matrix.col(i);
		rev_st_matrix.col(j + 1) = st_matrix.col(i) * -1;
		j = j + 2;
	}
	rev_st_matrix.rightCols(total_reac - rev_reac) = st_matrix.rightCols(total_reac - rev_reac);
}

void DetailedNetwork::calc_num_elst(){
	int ne;
	km_met = total_met;
	for (auto it = reactions.begin(); it != reactions.end(); it++){
		if (it->elem_steps_decomp()){
			ne = it->get_num_elem_steps();
			km_reac = km_reac + ne;
			km_met = km_met + ne / 2;
		}
		else
			km_reac++;
	}
}

deque<bool> DetailedNetwork::gen_rev_vector(){
	deque<bool> r;
	for (int i = 0; i < total_reac; i++)
		r.push_back(reactions[i].is_rev());
	return r;
}

deque<bool> DetailedNetwork::get_rev_vector(){
	return gen_rev_vector();
}

vector<int> DetailedNetwork::get_ids_upt_reactions(){
	vector<int> ids;
	for (int i = 0; i < upt_reac; i++){
		ids.push_back(total_reac - mes_reac + i);
		//cout << total_reac - mes_reac + i << endl;        //print para teste
	}
	return ids;
}

vector<int> DetailedNetwork::get_ids_prod_reactions(){
	vector<int> ids;
	for (int i = 0; i < upt_reac; i++){
		ids.push_back(total_reac - mes_reac + upt_reac + i);
		//cout << total_reac - mes_reac + upt_reac + i << endl;     //print para teste
	}
	return ids;
}

std::vector<std::string> DetailedNetwork::get_list_reactions(){
	vector<string> reac_names;
	for (auto it = reactions.begin(); it != reactions.end(); it++){
		reac_names.push_back(it->get_name());
	}
	return reac_names;
}

void DetailedNetwork::set_v_exp(Eigen::VectorXd v){
	if (v.size() == mes_reac){
		for (size_t i = 0; i < mes_reac; i++){
			reactions[rev_reac + irrev_reac + i].add_flux(v(i));
		}
	}
	else
		throw length_error("Size of the vector must match the number of measured reactions declared");
}

Eigen::VectorXd DetailedNetwork::get_v_exp(){
	if (exp){
		Eigen::VectorXd v(mes_reac);
		for (size_t i = 0; i < mes_reac; i++){
			v(i) = reactions[rev_reac + irrev_reac - mes_reac + i].get_flux();
		}
		//cout << v << endl;                       //print para teste
		return v;
	}
	else{
		return Eigen::VectorXd::Zero(mes_reac);
	}
}

Eigen::VectorXd DetailedNetwork::get_v_exp_rec(){
	if (exp){
		Eigen::MatrixXd M = get_mes_st_matrix().rightCols(mes_reac);
		Eigen::VectorXd v = M * get_v_exp();
		Eigen::VectorXd w(mes_reac);
		for (size_t i = 0; i < mes_reac; i++){
			w(i) = reactions[rev_reac + irrev_reac - mes_reac + i].get_var();
		}
		Eigen::MatrixXd W = (M.cwiseAbs() * w).asDiagonal();
		Eigen::MatrixXd A = Eigen::MatrixXd::Zero(1, v.size());
		// Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2, v.size());
		for (size_t i = 0; i < v.size(); i++){
			// A(1, i) = OxNum.transpose() * metabolites[total_met - exp_met + i].get_composition();
			A(0, i) = metabolites[total_met - exp_met + i].get_carbon();
		}
		Eigen::VectorXd v_rec;
		v_rec = v - W * A.transpose() * (A * W * A.transpose()).colPivHouseholderQr().inverse() * A * v;
		//cout << endl << get_v_exp() << endl;
		//cout << endl << M << endl;
		//cout << endl << v << endl;
		//cout << endl << A << endl;
		//cout << endl << W << endl;
		//cout << endl << v_rec << endl;
		//cout << endl << v_rec.transpose() * M << endl;
		return (v_rec.transpose() * M).transpose();
	}
	else{
		return Eigen::VectorXd::Zero(mes_reac);
	}
}




