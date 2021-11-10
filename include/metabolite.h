#include <string>
#include <array>
#include <iostream>
#include <Eigen/Dense>

#ifndef METABOLITE_H
#define METABOLITE_H

class Metabolite{
private:
	bool ext, exp;
	std::string name;
	Eigen::VectorXi composition;    //C, H, O, N, P, S 
public:
	Metabolite(std::string nome, bool e) : name(nome), ext(e), exp(false), composition(Eigen::VectorXi::Zero(6)){};
	void add_carbon(int c){ composition[0] = c; };
	void add_composition(Eigen::VectorXi c){ composition = c; };
	void set_as_exp(bool e){ exp = e; };
	std::string get_name() const { return name; };
	int get_carbon() const { return composition(0); };
	bool is_ext() const { return ext; };
	bool is_exp() const { return exp; };
	Eigen::VectorXi get_composition() const { return composition; };
};

#endif