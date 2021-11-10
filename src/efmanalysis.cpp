#include "stdafx.h"
#include "efmanalysis.h"

using namespace efmanalysis;
using namespace std;

// Functions for processing and analysing the EFMs

/*

*/
EFMAnalysis::EFMAnalysis(Eigen::MatrixXd efm, Eigen::MatrixXd sext, vector<int> ids_upt){
	EFM = efm;
	Sext = sext;
	int j = -1;
	cout << EFM.rows();
	for (int i = 0; i < EFM.rows(); i++){
		for (auto it = ids_upt.begin(); it != ids_upt.end(); it++){
			if (EFM(i, *it) != 0){
				j = *it;
				break;
			}
		}
		if (j > 0){
			EFM.row(i) = EFM.row(i) / EFM(i, j);
			j = -1;
		}
	}
	cout << EFM.rows();
}

/*

*/
void EFMAnalysis::gen_cycle_efm(){
	if (G.isZero())
		G = EFM*Sext.transpose();
	int size = G.rows();
	for (size_t i = 0; i < size; i++){
		if (G.row(i).isZero(1e-12)){
			cycle_efm.push_back(i);
		}
	}
}

/*

*/
void EFMAnalysis::gen_efm_families(){
	if (G.isZero())
		G = EFM*Sext.transpose();
	//cout << endl << G << endl;                                        //print para teste
	efm_families = group_same_rows(G);
	//for (int i = 0; i < efm_families.size(); i++){                    //print para teste
	//	for (int j = 0; j < efm_families[i].size(); j++){             //print para teste
	//		cout << efm_families[i][j] << " ";                        //print para teste
	//	}                                                             //print para teste
	//	cout << endl;                                                 //print para teste
	//}                                                                 //print para teste
}

/*

*/
std::vector<std::size_t> EFMAnalysis::get_cycle_efm(){
	if (cycle_efm.empty())
		gen_cycle_efm();
	return cycle_efm;
}

/*

*/
vector<vector<size_t>> EFMAnalysis::get_efm_families(){
	if (efm_families.empty())
		gen_efm_families();
	return efm_families;
}

vector<size_t> EFMAnalysis::get_efm_rank(vector<int> ids){
	if (efm_rank.empty())
		gen_efm_rank(ids);
	return efm_rank;
}

/*

*/
vector<size_t> EFMAnalysis::get_efms_reaction(){
	vector<size_t> efm_reac;
	for (size_t i = 0; i < EFM.cols(); i++){
		efm_reac.push_back((EFM.col(i).array() != 0).count());
	}
	return efm_reac;
}

/*

*/
void EFMAnalysis::group_efm_by_prod(vector<int> ids){           //a funcao vai considerar que a ordem das reacoes e a ordem de importancia
	size_t size = ids.size();
	group.resize(size+1);
	for (size_t j = 0; j < EFM.rows(); j++){
		for (size_t i = 0; i < size; i++){
			if (EFM(j, ids[i]) != 0){
				group[i].push_back(make_pair(j,EFM(j,ids[i])));
				break;
			}
		}
	}
	for (int i = 0; i < group.size(); i++){                                            //print para teste
		size = group[i].size();
		for (int j = 0; j < size; j++){                                    //print para teste
			cout << "(" << group[i][j].first << ", " << group[i][j].second << ") ";   //print para teste
		}                                                                             //print para teste
		cout << endl << endl;                                                                 //print para teste
	}                                                                                 //print para teste
}

/*

*/
void EFMAnalysis::gen_efm_rank(vector<int> ids){
	if (cycle_efm.empty())
		gen_cycle_efm();
	efm_rank.insert(efm_rank.end(), cycle_efm.begin(), cycle_efm.end());
	group_efm_by_prod(ids);
	for (auto it = group.rbegin(); it != group.rend(); it++){
		sort(it->begin(), it->end(), [](pair<size_t, double> &left, pair<size_t, double> &right) {
			return left.second < right.second;
		});
		for (auto itt = it->begin(); itt != it->end(); itt++){
			efm_rank.push_back(itt->first);
		}
	}	
}


/*

*/
Eigen::MatrixXd EFMAnalysis::get_unique_efm(const int mes_reac){     //ARRUMAR
	if (efm_families.empty())
		gen_efm_families();
	size_t size = efm_families.size();
	Eigen::MatrixXd M(size, mes_reac);
	for (size_t i = 0; i < size; i++){
		M.row(i) = G.row(efm_families[i][0]).tail(mes_reac);
	}
	return M;
}

/*

*/
Eigen::MatrixXd EFMAnalysis::get_efm_no_cycle(){
	Eigen::MatrixXd EFMnc;
	EFMnc = EFM;
	if (cycle_efm.empty())
		gen_cycle_efm();
	for (auto it = cycle_efm.rbegin(); it != cycle_efm.rend(); it++){
		remove_row(EFMnc, *it);
	}
	return EFMnc;
}

/*

*/
Eigen::MatrixXd EFMAnalysis::get_ranked_efm(vector<int> ids){
	Eigen::MatrixXd EFMn(EFM.rows(), EFM.cols());
	if (efm_rank.empty())
		gen_efm_rank(ids);
	for (size_t i = 0; i < efm_rank.size(); i++){
		EFMn.row(i) = EFM.row(efm_rank[i]);
	}
	EFM = EFMn;
	G = EFM*Sext.transpose();
	cout << G;
	return EFMn;
}

vector<double> EFMAnalysis::get_max_contribution(Eigen::VectorXd v, bool families){
	if (v.isZero())
		return { 0.0 };
	int size = EFM.rows();
	int reac = EFM.cols();
	int exp = v.size();
	double cont;
	vector<double> max_cont;
	v = v / (v(0) + v(1));     //TODO: A DIVISAO DEVE SER PELA SOMA DE FLUXOS DE TODOS OS SUBSTRATOS QUE SAO FONTE DE CARBONO
	for (int i = 0; i < size; i++){
		cont = 1;
		for (int j = 0; j < exp; j++){
			if (abs(EFM(i, reac - exp + j)) > abs(v(j))){
				if (v(j) / EFM(i, reac - exp + j) < cont){
					cont = v(j) / EFM(i, reac - exp + j);
				}
			}
		}
		max_cont.push_back(cont);
	}
	return max_cont;
}