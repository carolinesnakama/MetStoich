#include "stdafx.h"
#include "interface.h"

using namespace std;
using namespace stoichiometry;
using namespace network;
using namespace efmanalysis;
using namespace mfa;
using namespace efv;

/**
interface.cpp
Purpose: connects all classes among each other.
Upon setting a network, it is possible to perform different calculations.

@author Caroline S. M. Nakama
@version 2.0 01/02/18
*/

//Funcoes de Interface

/*

*/
Interface::Interface(Eigen::MatrixXd S, bool intern){
	rede = unique_ptr<Network>(new Network(S, intern));
	check_network();
}

/*

*/
Interface::Interface(Eigen::MatrixXd S, bool intern, deque<bool> r){
	rede = unique_ptr<Network>(new Network(S, intern, r));
	check_network();
}

/*

*/
Interface::Interface(string file){
	rede = unique_ptr<DetailedNetwork>(new DetailedNetwork(file));
	check_network();
}

/*
Assigns a new stoichiometry matrix and deletes any existing calculation

@param S stoichiometry matrix
       intern boolean that indicates if the stoichiometry matrix is composed of only internal metabolites
*/
void Interface::addStoichiometryMatrix(Eigen::MatrixXd S, bool intern){
	if (!detailed){
		rede->set_st_matrix(S);
		clearCalculations();
		rede->check_st_matrix();
	}
	else
		throw runtime_error("It is not possible to directly change the stoichiometry matrix. Please update the Network file");
}

/*
Assigns a vector of reversibility of each reation in S and deletes EFM calculations
Throws a runtime error if N is not defined and
       a length error if the number of entries is different that the number of reactions in N

@param r vector of booleans that indicates the reversibility of each reaction
         in the same order it appears in S
*/
void Interface::addReversibilityVector(deque<bool> r){
	if (!detailed){
		try{
			rede->set_rev_vector(r);
			try{
				clearCalculation("EFM");
				clearCalculation("EFMNS");
			}
			catch (runtime_error e){}
		}
		catch (runtime_error e){
			cout << e.what();
		}
		catch (length_error e){
			cout << e.what();
		}
	}
	else
		throw runtime_error("It is not possible to directly change the reversibility vector.Please update the Network file");
}

/*
ARRUMAR
*/
void Interface::addFluxVector(Eigen::VectorXd v){
	/*F.reset(new expdata::ExpData());
	vector<size_t> pos;
	if (N.isZero()){
		throw runtime_error("Please add a Stoichiometry matrix first");
	}
	if (F->setV(N, v)){
		pos = F->getExtReactions();
		cout << "Fluxes successful added for reactions in positions ";
		for (auto it = pos.begin(); it != pos.end(); it++)
			cout << *it << " ";
		cout << "of the stoichiometry matrix in this order";
	}
	else{
		throw length_error("The flux vector provided does not have a entry for all boundary reactions.");
	}*/
}


/*
Deletes the calculation passed as argument
Throws a runtime error if the entry is not found

@param name the name of the stoichiometric calculation to be deleted 
*/
void Interface::clearCalculation(string name){
	transform(name.begin(), name.end(), name.begin(), ::toupper);
	if (stoich.erase(name) == 0)
		throw runtime_error("The name provided is not a valid calculation or calculation has not been performed");
}

/*
Deletes all calculations from stoich_calc
*/
void Interface::clearCalculations(){
	stoich.clear();
}

bool Interface::calcEFV(){
	bool add = false;
	deque<bool> rev;
	rev = rede->get_rev_vector();
	for (int i = 0; i <= 2 * rede->get_num_mes_reaction(); i++)
		rev.push_back(false);
	if (!EV){
		EV = unique_ptr<EFV>(new EFV(rede->get_int_st_matrix(), rede->get_v_exp_rec(), rede->get_num_mes_reaction()));
	}
	try{
		add = (stoich.insert({ "EFV", unique_ptr<StoichAnalysis>(new EFM(EV->getM(), rev)) })).second;
	}
	catch (length_error e){
		cout << e.what() << endl;
	}
	return add;
}

/*
Inserts the Conservation Relation calculation to the list of possible stoichiometric calculations.
If it already exists, it does nothing and returns false.

@return a boolean true if the calculation was successfully added
                  false otherwise
*/
bool Interface::calcPosConsRel(){
	bool add = false;
	try{
		add = (stoich.insert({ "CONSREL", unique_ptr<StoichAnalysis>(new PosConsRel(rede->get_int_st_matrix())) })).second;
	}
	catch (length_error e){
		cout << e.what() << endl;
	}
	return add;
}

/*
Inserts the EFM calculation to the list of possible stoichiometric calculations.
If it already exists, it does nothing and returns false.

@return a boolean true if the calculation was successfully added
false otherwise
*/
bool Interface::calcEFM(){
	bool add = false;
	try{
		add = (stoich.insert({ "EFM", unique_ptr<StoichAnalysis>(new EFM(rede->get_int_st_matrix(), rede->get_rev_vector())) })).second;
	}
	catch (length_error e){
		cout << e.what() << endl;
	}
	catch (runtime_error e){
		cout << e.what() << endl;
	}
	return add;
}

/*
Inserts the EFM calculation with the null space method to the list of possible stoichiometric calculations.
If it already exists, it does nothing and returns false.

@return a boolean true if the calculation was successfully added
false otherwise
*/
bool Interface::calcEFMNS(){
	bool add = false;
	try{
		add = (stoich.insert({ "EFMNS", unique_ptr<StoichAnalysis>(new EFMNS(rede->get_int_st_matrix(), rede->get_rev_vector())) })).second;
	}
	catch (length_error e){
		cout << e.what() << endl;
	}
	catch (runtime_error e){
		cout << e.what() << endl;
	}
	return add;
//	auto add = stoich.insert({ "EFMNS", unique_ptr<StoichAnalysis>(new EFMNS(N, R)) });
//	return add.second;
}

/*

*/
bool Interface::calcMFA(){
	//bool add = false;

	//if (!N.isZero() && !(F != nullptr)){

	//}
	//return add;
	return false;
}


/*
Performs the calculation specified in name using N and R. 

@param name the calculation to be performed {"CONSREL", "EFM", "EFMNS"}
@return a tableau that provides the information calculated
*/
Eigen::MatrixXd Interface::calculate(string name){
	bool added;
	try{
		if (name == "CONSREL")
			added = calcPosConsRel();
		else if (name == "EFM")
			added = calcEFM();
		else if (name == "EFMNS")
			added = calcEFMNS();
		else if (name == "EFV")
			added = calcEFV();
		else
			throw invalid_argument("The name provided for calculation is not valid");
	}
	catch (runtime_error e){
		throw runtime_error(e.what());
	}
	if (added)
		stoich[name]->run();
	return stoich[name]->getT();
}

//Funcoes de InterfaceNetwork

/*
Checks if the network passes the tests for consistency
Throws an exception if the network does not pass consistent and carbon verifications
*/
void Interface::check_network(){
	try{
		rede->check_st_matrix();
	}
	catch (exception e){
		cout << e.what() << endl;
	}
}

/*
updates the metabolic network by creating a new instance of Network and checks it consistency

@param file a string with the complete address of the file that contains the new metabolic network
*/
void Interface::updateNetwork(string file){
	rede.reset(new DetailedNetwork(file));
	check_network();
	clearCalculations();
}

Eigen::MatrixXd Interface::getStoichResult(std::string name){
	if (stoich.find(name) != stoich.end()){
		return stoich[name]->getT();
	}
	else
		return Eigen::MatrixXd::Zero(1, 1);
}

void Interface::printVector(ofstream &file, const vector<size_t> &vector){
	for (auto it = vector.begin(); it != vector.end(); it++){
		file << *it << " ";
	}
}

void Interface::printVectorVector(ofstream &file, const vector<vector<size_t>> &vector){
	int i = 1;
	for (auto it = vector.begin(); it != vector.end(); it++){
		file << "F" << i << ": ";
		printVector(file, *it);
		file << endl;
		i++;
	}
}

/*

*/
void Interface::generateOutputFile(string file, Eigen::MatrixXd EFM){

	// Eigen::MatrixXd T;
	// if (EFM.isZero())
	// 	T = calculate("EFM");
	// else
	// 	T = EFM;

	// if (!EA)
	// 	EA = unique_ptr<EFMAnalysis>(new EFMAnalysis(T, rede->get_mes_st_matrix(), rede->get_ids_upt_reactions()));

	// ofstream output(file + "1-EFM.txt");
	// vector<vector<int>> gr = stoich["EFM"]->get_GR();
	// vector<string> list_reac = rede->get_list_reactions();
	// output << "Stoichiometry matrix" << endl;
	// output << rede->get_st_matrix() << endl << endl;
	// output << "Blocked reactions" << endl << "BR: ";
	// for (int j = 0; j < gr[0].size(); j++){
	// 	output << list_reac[gr[0][j]] << " ";
	// }
	// output << endl << "Coupled reactions";
	// for (int i = 1; i < gr.size(); i++){
	// 	output << endl << "CR" << i << ": ";
	// 	for (int j = 0; j < gr[i].size(); j++){
	// 		output << list_reac[gr[i][j]] << " ";
	// 	}
	// }
	// output << endl << endl << "Reduced stoichiometry matrix" << endl;
	// output << stoich["EFM"]->getNred() << endl;
	// output << "Reversibility of coupled reactions" << endl;
	// deque<bool> crrev = stoich["EFM"]->getRred();
	// for (auto it = crrev.begin(); it != crrev.end(); it++){
	// 	output << *it << " ";
	// }
	// output << endl << endl << "Elementary flux modes" << endl;
	// output << EA->get_ranked_efm(rede->get_ids_prod_reactions()) << endl;
	// output << "This network has " << EA->get_norm_efm().rows() << " EFMs and the first " << stoich["EFM"]->get_rev_efms() << " is(are) reversible." << endl;
	// output << "The order of the reactions follows the order they were declared in the input file. First ENZREV, then ENZIRREV and last ENZMEAS." << endl;
	// output.close();

	// ofstream output2(file + "2-EFMA.txt");
	// output2 << "Thermodynamically infeasible EFMs: ";
	// printVector(output2, EA->get_cycle_efm());
	// output2 << endl << endl << "Families of EFMs with the same overall stoichiometry" << endl;
	// printVectorVector(output2, EA->get_efm_families());
	// output2 << endl << "Overall stoichiometry for each family. Metabolite order are the same as of METEXT in input file." << endl;
	// output2 << EA->get_unique_efm(rede->get_num_mes_reaction()) << endl;
	// output2 << endl << "Maximum contribution of each EFM" << endl;
	// vector<double> mc;
	// mc = EA->get_max_contribution(rede->get_v_exp());
	// for (auto it = mc.begin(); it != mc.end(); it++){
	// 	output2 << "EFM " << it - mc.begin() << ": " << *it << endl;
	// }
	// output2 << endl << "Quantity of EFMs that each reaction is present" << endl;
	// vector<size_t> efm_reac = EA->get_efms_reaction();
	// for (size_t i = 0; i < efm_reac.size(); i++){
	// 	output2 << list_reac[i] << " - " << efm_reac[i] << endl;
	// }
	// output2.close();

	// ofstream output3(file + "3-PCR.txt");
	// output3 << "Positive conservation Relations" << endl;
	// output3 << calculate("CONSREL") << endl;
	// output3 << "The order of the metabolites are the same as that declared under METINT in the input file." << endl;
	// //rede->print_metabolites(output3);
	// output3.close();

	// Eigen::MatrixXd S(3,6);
	// Eigen::VectorXd v(6);
	// S.row(0) = EA->get_unique_efm(rede->get_num_mes_reaction()).row(1).cwiseAbs();
	// S.row(1) = EA->get_unique_efm(rede->get_num_mes_reaction()).row(4).cwiseAbs();
	// S.row(2) = EA->get_unique_efm(rede->get_num_mes_reaction()).row(5).cwiseAbs();
	// //M = unique_ptr<MFA>(new MFA(S.transpose(), rede->get_v_exp()));
	// v << 1.97, 0, 7.14, 4.26, 1.6, 0;
	// M = unique_ptr<MFA>(new MFA(S.transpose(), v));
	// M->run();
	// ofstream output4(file + "4-MFA.txt");
	// output4 << M->getN() << endl << endl;
	// output4 << M->getV() << endl << endl;
	// output4 << M->getX() << endl << endl;
	// output4 << S.transpose()*M->getX() << endl;

	ofstream output5(file + "5-EFV.txt");
	output5 << "Elementary flux vector" << endl;
	output5 << endl << calculate("EFV") << endl;
	//output5 << endl << stoich["EFV"]->augTableauToMatrix() << endl;
	output5 << endl << "Augmented matrix" << endl;
	output5 << EV->getM() << endl;
	
}

/*

*/
void Interface::analyseEFM(){
	if (!EA)
		EA = unique_ptr<EFMAnalysis>(new EFMAnalysis(stoich["EFM"]->getT(), rede->get_mes_st_matrix(), rede->get_ids_upt_reactions()));
	cout << EA->get_unique_efm(rede->get_num_mes_reaction()) << endl;
	cout << EA->get_efm_no_cycle() << endl;
}


