#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <deque>
#include <list>

#ifndef REACTION_H
#define REACTION_H

class Reaction{
private:
	bool rev, mes, kin;
	int subs, prod, fc, cs, cp;
	double flux, var;
	std::string name;
	std::vector<std::pair<std::string, double>> metabolites;
	std::vector<std::string> enz_comp;
	//std::vector<std::unique_ptr<Reaction>> elem_steps;
	//std::vector<std::pair<int, double>> elem_steps;
	//std::vector<std::vector<std::pair<std::string, int>>> elem_steps;
	//std::vector<std::vector<int>> elem_steps;     //vetor de 4 ints: id da etapa elementar, met (0) ou enz_comp (1), id do met ou enz_comp, subs ou prod  
	std::vector<std::string> elem_steps;
	std::vector<std::pair<std::string, int>> elem_steps_met_subs;
	std::vector<std::pair<std::string, int>> elem_steps_met_prod;
	std::vector<std::pair<int, int>> elem_steps_enz_subs;
	std::vector<std::pair<int, int>> elem_steps_enz_prod;
	//std::vector<std::pair<int, int>> elem_steps_met_subs;
	//std::vector<std::pair<int, int>> elem_steps_met_prod;
	std::deque<bool> form_complex;
	void list_enz_complexes();
public:
	Reaction(){};
	Reaction(std::string nome, bool r, bool k, bool m);
	void add_metabolite(std::string met, double coef, bool c = false);
	void add_flux(double v){ flux = v; };
	void add_var(double w){ var = w; };
	void decomp_elem_steps();
	//void print_reaction(std::ofstream &file);
	std::string get_name(){ return name; };
	bool is_rev(){ return rev; };
	bool measured() { return mes; };
	bool elem_steps_decomp() { return kin; };
	double get_flux() { return flux; };
	double get_var() { if (var == 0) return flux * 0.025; else return var; };
	double get_coefficient(std::string met);
	int get_num_elem_steps(){ if (kin) return (cs + cp + 1) * 2; else return 1; };
	std::vector<std::string> get_elst_list(){ return elem_steps; };
	std::vector<std::string> get_enz_list(){ return enz_comp; };
	std::vector<std::pair<std::string, double>>::const_iterator get_metabolite_begin() const { return metabolites.begin(); };
	std::vector<std::pair<std::string, double>>::const_iterator get_metabolite_end() const { return metabolites.end(); };
	std::vector<std::pair<std::string, int>>::const_iterator get_elem_steps_met_begin(bool subs) const { if (subs) return elem_steps_met_subs.begin(); else return elem_steps_met_prod.begin(); };
	std::vector<std::pair<std::string, int>>::const_iterator get_elem_steps_met_end(bool subs) const { if (subs) return elem_steps_met_subs.end(); else return elem_steps_met_prod.end(); };
	std::vector<std::pair<int, int>>::const_iterator get_elem_steps_enz_begin(bool subs) const { if (subs) return elem_steps_enz_subs.begin(); else return elem_steps_enz_prod.begin(); };
	std::vector<std::pair<int, int>>::const_iterator get_elem_steps_enz_end(bool subs) const { if (subs) return elem_steps_enz_subs.end(); else return elem_steps_enz_prod.end(); };
};


#endif