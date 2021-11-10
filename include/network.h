#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <string>
#include <Eigen/Dense>
#include <map>
#include "metabolite.h"
#include "reaction.h"
#include "auxmatrixop.h"

#ifndef NETWORK_H
#define NETWORK_H

namespace network{
	class Network{
	private:
		std::deque<bool> rev_vector;
		std::vector<std::size_t> ext_met_vector;
		Eigen::VectorXd v_exp;
		virtual void gen_int_st_matrix();
		void gen_ext_met_vector();
		int mes_reac;
	protected:
		Eigen::MatrixXd st_matrix, int_st_matrix, ext_st_matrix, mes_st_matrix, elst_st_matrix, elst_subs_st_matrix, rev_st_matrix;
		//std::vector<std::size_t> border_reac;
		void gen_border_reac();
	public:
		Network(){};
		~Network(){};
		Network(Eigen::MatrixXd S, bool intern = false);
		Network(Eigen::MatrixXd S, bool intern, std::deque<bool> r);
		//std::vector<std::size_t> get_border_reac();
		void clear();
		virtual void set_v_exp(Eigen::VectorXd v);
		virtual std::deque<bool> get_rev_vector();
		virtual Eigen::VectorXi get_met_carb();
		virtual void set_st_matrix(Eigen::MatrixXd S, bool intern = false);
		virtual void set_rev_vector(std::deque<bool> r);
		virtual void print_metabolites(std::ofstream &file);
		virtual void print_reactions(std::ofstream &file);
		virtual void check_st_matrix();
		virtual void clear_rev_vector();
		virtual int get_num_ext_met();
		virtual int get_num_int_met();
		virtual int get_num_total_reaction();
		virtual int get_num_rev_reaction();
		virtual int get_num_irrev_reaction();
		virtual int get_num_mes_reaction() const { return mes_reac; };
		virtual std::vector<int> get_ids_upt_reactions();   //TODO
		virtual std::vector<int> get_ids_prod_reactions();   //TODO
		virtual Eigen::MatrixXd get_st_matrix();
		virtual Eigen::MatrixXd get_int_st_matrix();
		virtual Eigen::MatrixXd get_ext_st_matrix();
		virtual Eigen::MatrixXd get_mes_st_matrix(){ return mes_st_matrix; };
		virtual Eigen::MatrixXd get_elst_st_matrix(){ return Eigen::MatrixXd::Zero(1,1); };
		virtual Eigen::MatrixXd get_elst_subs_st_matrix(){ return Eigen::MatrixXd::Zero(1, 1); };
		virtual Eigen::MatrixXd get_rev_st_matrix() { return Eigen::MatrixXd::Zero(1, 1); };
		virtual std::vector<std::string> get_list_reactions();
		virtual Eigen::VectorXd get_v_exp() { return v_exp; };
		virtual Eigen::VectorXd get_v_exp_rec() { return v_exp; };
	};

	class DetailedNetwork : public Network{
	private:
		Eigen::VectorXi OxNum;    //C, H, O, N, P, S 
		//atributos
		//Eigen::MatrixXd elst_st_matrix, elst_subs_st_matrix;
		Eigen::VectorXi met_carb;
		bool carbon, exp, redox, kmodel;
		int total_reac, rev_reac, irrev_reac, mes_reac, upt_reac, km_reac;
		int total_met, int_met, ext_met, exp_met, km_met;
		std::vector<Metabolite> metabolites;
		std::vector<Reaction> reactions;
		std::vector<std::string> elst_enz_list;
		std::vector<std::string> elst_reac_list;
		//funcoes
		void read_file(std::string file);
		void add_reaction(std::vector<std::string> list);
		void add_met_composition(std::vector<std::string> &list);
		void gen_st_matrix();
		std::deque<bool> gen_rev_vector();
		virtual void gen_int_st_matrix();
		void gen_elst_st_matrix();
		void gen_mes_st_matrix();
		void gen_rev_st_matrix();
		void calc_num_elst();
		bool check_carb_bal(std::ofstream& log);
		void check_exp_bal(std::ofstream& log);
		void gen_carb_vector();
		void gen_vector(std::string& line, std::vector<std::string>& reaction);
		void pop_enz_vector(std::vector<std::string>& vetor, bool rev, bool mes, bool k);
		void pop_met_vector(std::vector<std::string>& vetor, std::vector<int>& carb, bool ext);
		void get_carbon(std::vector<std::string>& vetor, std::vector<std::string>& met, std::vector<int>& carb, bool carbon);
		void class_ext_met(std::vector<std::string>& vetor);
		// void class_enz_mes(std::vector<std::string>& vetor, std::vector<std::string>& vetor2 = std::vector<std::string>());
		void class_enz_mes(std::vector<std::string>& vetor);
		int find_reaction_position(std::string name);
		int find_met_position(std::string name);
		double balance(bool carb);
		double carbon_balance(int pos, double flux);
		double redox_balance(int pos, double flux);
	public:
		DetailedNetwork() = delete;
		~DetailedNetwork(){};
		DetailedNetwork(std::string file);
		virtual std::deque<bool> get_rev_vector();
		virtual void set_v_exp(Eigen::VectorXd v);
		virtual void set_st_matrix(Eigen::MatrixXd S, bool intern = false){}; //COLOCAR UMA MENSAGEM DE ERRO
		virtual void set_rev_vector(std::deque<bool> r){};   //COLOCAR UMA MENSAGEM DE ERRO
		virtual void clear_rev_vector(){};    //COLOCAR UMA MENSAGEM DE ERRO
		virtual Eigen::VectorXi get_met_carb() const { return met_carb; };
		virtual void print_metabolites(std::ofstream &file);
		virtual void print_reactions(std::ofstream &file);
		virtual void check_st_matrix();
		virtual Eigen::MatrixXd get_st_matrix() const { return st_matrix; };
		virtual Eigen::MatrixXd get_int_st_matrix() const { return int_st_matrix; };
		virtual Eigen::MatrixXd get_ext_st_matrix() const { return ext_st_matrix; };
		virtual Eigen::MatrixXd get_mes_st_matrix() const { return mes_st_matrix; };
		virtual Eigen::MatrixXd get_elst_st_matrix() { if (elst_st_matrix.isZero()) gen_elst_st_matrix(); return elst_st_matrix; };
		virtual Eigen::MatrixXd get_elst_subs_st_matrix() { /*if (elst_st_matrix.isZero()) gen_elst_st_matrix();*/ return elst_subs_st_matrix; };
		virtual Eigen::MatrixXd get_rev_st_matrix() { if (rev_st_matrix.isZero()) gen_rev_st_matrix(); return rev_st_matrix; };
		virtual int get_num_ext_met() const { return ext_met; };
		virtual int get_num_int_met() const { return int_met; };
		virtual int get_num_total_reaction() const { return total_reac; };
		virtual int get_num_rev_reaction() const { return rev_reac; };
		virtual int get_num_irrev_reaction() const { return irrev_reac; };
		virtual int get_num_mes_reaction() const { return mes_reac; };
		virtual std::vector<int> get_ids_upt_reactions();
		virtual std::vector<int> get_ids_prod_reactions();
		virtual std::vector<std::string> get_list_reactions();
		virtual Eigen::VectorXd get_v_exp();
		virtual Eigen::VectorXd get_v_exp_rec();
	};
}
#endif