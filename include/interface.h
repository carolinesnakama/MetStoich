#include <Eigen/Dense>
#include <vector>
#include <deque>
#include <string>
#include <memory>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <fstream>
#include "stoichiometry.h"
#include "network.h"
#include "expdata.h"
#include "mfa.h"
#include "efmanalysis.h"
#include "efv.h"
#include "printvector.h"

#ifndef INTERFACE_H
#define INTERFACE_H

class Interface{
private:
	std::map<std::string, std::unique_ptr<stoichiometry::StoichAnalysis>> stoich;
	std::unique_ptr<mfa::MFA> M;
	std::unique_ptr<network::Network> rede;
	std::unique_ptr<expdata::ExpData> F;
	std::unique_ptr<efmanalysis::EFMAnalysis> EA;
	std::unique_ptr<efv::EFV> EV;
	bool detailed = true;
	bool calcPosConsRel();
	bool calcEFM();
	bool calcEFMNS();
	bool calcMFA();
	bool calcEFV();
	void check_network();
	void printVector(std::ofstream &file, const std::vector<std::size_t> &vector);
	void printVectorVector(std::ofstream &file, const std::vector<std::vector<std::size_t>> &vector);
public:
	Interface(){};
	~Interface(){};
	Interface(Eigen::MatrixXd S, bool intern = false);                  
	Interface(Eigen::MatrixXd S, bool intern, std::deque<bool> r);      
	Interface(std::string file);                                          
	void addStoichiometryMatrix(Eigen::MatrixXd S, bool intern = false);
	void addReversibilityVector(std::deque<bool> r);
	void addFluxVector(Eigen::VectorXd v);
	void updateNetwork(std::string file);                               
	void clearCalculation(std::string calc);
	void clearCalculations();
	Eigen::MatrixXd calculate(std::string name);
	Eigen::MatrixXd getStoichResult(std::string name);
	Eigen::MatrixXd getStoichMatrix(){ return rede->get_st_matrix(); };
	Eigen::MatrixXd getElstStoichMatrix(){ return rede->get_elst_st_matrix(); };
	Eigen::MatrixXd getElstSubsStoichMatrix(){ return rede->get_elst_subs_st_matrix(); };
	Eigen::MatrixXd getRevStoichMatrix() { return rede->get_rev_st_matrix(); };
	void generateListMetabolites(std::ofstream &file){ rede->print_metabolites(file); };
	void generateListReactions(std::ofstream &file){ rede->print_reactions(file); };
	void generateOutputFile(std::string file, Eigen::MatrixXd EFM = Eigen::MatrixXd::Zero(1,1));
	void analyseEFM();
};



#endif