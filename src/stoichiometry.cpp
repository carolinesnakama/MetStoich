#include "stdafx.h"
#include "stoichiometry.h"

using namespace std;
using namespace stoichiometry;

/**
stoichanalysis.cpp
Purpose: perform stoichometric analysis on a stoichiometry matrix provided. 
Analysis includes positive conservation relations, EFM by Schuster's method and EFM using the nullspace method.

@author Caroline S. M. Nakama
@version 1.2 04/07/18
*/

// Funcoes da classe StoichAnalysis

/*
Returns a vector with (column) indices of N corresponding to blocked reactions

@param remove if columns corresponding to blocked reactions must be removed
@return a vector with indices of blocked reactions
*/
vector<int> StoichAnalysis::get_blocked_reactions(const bool remove){
	vector<int> b;
	Eigen::FullPivLU<Eigen::MatrixXd> LU;
	LU = N.fullPivLu();
	//LU.setThreshold(tol);
	K = LU.kernel();
	//cout << "Kernel" << endl << K << endl << endl;    //print para teste
	int size = K.rows();
	double min;
	for (int i = 0; i < size; i++){
		if (K.row(i).isZero(tol)){
			b.push_back(i);
			CRcoeff.push_back(0);
		}
		else{
			min = minAbsCoeff(K.row(i).transpose());
			K.row(i) = K.row(i) / min;
			CRcoeff.push_back(min);
		}
	}
	if (remove)
		remove_rows(K, b);
	// << "Kernel sem blocked reactions" << endl << K << endl << endl;   //print para teste
	//cout << "Vetor com coupled reactions" << endl;                     //print para teste
	for (int i = 0; i < (int)CRcoeff.size(); i++)
		cout << CRcoeff[i] << endl;
	cout << endl;    //print para teste
	return b;
}

/*
Generates a vector with r entries in the same order as N. 
Same coefficient means the reactions are in the same group
*/
void StoichAnalysis::gen_coupled_reactions(){
	vector<int> b;
	int n = 0;
	b = get_blocked_reactions(true);
	int bsize = b.size();
	int size = K.rows();
	CR.resize(size, 0);
	for (int i = 0; i < size; i++){
		if (CR[i] == 0){
			n++;
			CR[i] = n;
			for (int k = i + 1; k < size; k++){
				if (K.row(i).isApprox(K.row(k), tol)){
					if (CR[k] == 0) CR[k] = n;
				}
			}
		}
	}
	if (!b.empty()){
		for (int i = 0; i < bsize; i++){
			if (b[i] > size)
				CR.push_back(0);
			else
				CR.insert(CR.begin() + b[i], 0);
		}
	}
	cout << "Vetor com coupled reactions" << endl;    //print para teste
	for (int i = 0; i < (int)CR.size(); i++)
		cout << CR[i] << endl;
	cout << endl;    //print para teste
}

/*
Returns the r vector with coupled reactions indicated by the value of the coefficients.

@return vector with coupled reations
*/
vector<int> StoichAnalysis::get_coupled_reactions(){
	if (CR.empty())
		gen_coupled_reactions();
	return CR;
}

/*
Code structure that calculates the chosen stoichiometric analysis. 
*/
void StoichAnalysis::run(){
	Tableau Tt;
	vector<int> ind;
	int size;
	Eigen::VectorXd v;
	cout << "aqui" << endl;
	//genTableau(N);
	// printT(); cout << endl << endl;  //print para teste
	Tt.allocate(T.size() * 1000);
	for (int k = 0; k < cols; k++){
		cout << "Tableau: " << k << endl;
		//int size = T.size();
		//for (int i = 0; i < size; i++){
		for (auto it = T.begin(); it != T.end(); it++){
			if ((*it)(k) > -tol && (*it)(k) < tol)
				addRow(it, Tt);                               
			else
				ind.push_back(it - T.begin());
		}
		size = ind.size();
		for (int i = 0; i < size; i++){
			for (int j = i + 1; j < size; j++){
				//v = combineRows(k, ind[i], ind[j], Tt.begin(), Tt.end(), back_inserter(Tt), inserter(Tt, Tt.begin()));
				//int n = ind[j];
				//int n = T.end() - T.begin();
				v = combineRows(k, T.begin() + ind[i], T.begin() + ind[j]);
				//cout << v.transpose() << endl;
				//cout << v.isZero() << endl;
				//bool b = v.isZero();
				//int n = v.size();
				if (v.size() > 0 && (Tt.begin() == Tt.end() || !check_condition(Tt.begin(), Tt.end(), v)))
					addRowToEnd(v, Tt, T.begin() + ind[i], T.begin() + ind[j]);
			}
		}
		updateT(Tt);
		ind.clear();
	}
	int n = T.size();
	int m = FluxR.size();
	//final_check();
	//changeValues();                                //para melhor vizualizacao
	//cout << endl;                                  //print para teste
	//for (int i = 0; i < (int)T.size(); i++){       //print para teste
	//	cout << T[i].transpose() << endl;          //print para teste 
	//}                                              //print para teste
	//cout << endl;                                  //print para teste
}

/*
Changes tableau coefficients so all entries are integers.
*/
void StoichAnalysis::changeValues(){
	//int size = T.size();
	double min;
	//for (int i = 0; i < size; i++){
	for (auto it = T.begin(); it != T.end(); it++){
		min = minAbsCoeff((*it));
		(*it) = (*it) * (1 / min);
		//		cout << endl << min << endl;
	}
}

/*
Separate coupled reactions into groups

@return a vector with indices of reactions that can be deleted to reduce N
        (keeps the index of the first reaction in each group)
*/
vector<int> StoichAnalysis::group_reactions(){
	vector<int> del;
	int size = CR.size();
	int g = *max_element(CR.begin(), CR.end());
	int a = 0;
	GR.resize(g + 1);
	for (int i = 0; i < size; i++){
		GR[CR[i]].push_back(i);
		if (CR[i] <= a)
			del.push_back(i);
		else
			a++;
	}
	for (int i = 0; i < (int)GR.size(); i++){             //print para teste
		cout << "grupo " << i << ": ";                    //print para teste
		for (int j = 0; j < (int)GR[i].size(); j++)       //print para teste
			cout << GR[i][j] << ", ";                     //print para teste
		cout << endl;                                     //print para teste
	}                                                     //print para teste
	return del;
}

/*
Reduces N by adding the columns of coupled reactions and deleting rows where all metabolites are canceled out
*/
void StoichAnalysis::reduceN(){
	vector<int> del;
	Nred = N;
	gen_coupled_reactions();
	del = group_reactions();
	size_t sizeg = GR.size();
	size_t size;
	for (size_t i = 1; i < sizeg; i++){
		size = GR[i].size();
		for (size_t j = 1; j < size; j++){
			Nred.col(GR[i][0]) += (CRcoeff[GR[i][j]] / CRcoeff[GR[i][0]]) * Nred.col(GR[i][j]);
		}
	}
	remove_cols(Nred, del);
	//cout << "N interno completo" << endl;          //print para teste
	//cout << N << endl;                             //print para teste
	//cout << "N com grupos de reacoes" << endl;     //print para teste
	//cout << Nred << endl << endl;                  //print para teste
	del.clear();
	size = Nred.rows();
	for (size_t i = 0; i < size; i++){
		if (Nred.row(i).isZero(tol))
			del.push_back(i);
	}
	remove_rows(Nred, del);
	//rows = Nred.rows(); cols = Nred.cols();     //TO DO: VERIFICAR SE ESSA LINHA E NECESSARIA
	cout << "N reduzida" << endl;               //print para teste
	cout << Nred << endl << endl;               //print para teste
}

/*
Reduces R by determining whether each group of reactions is reversible or not
If one reaction is irreversible in a group, then the whole group is irreversible
*/
void StoichAnalysis::reduceR(){
	//if (CR.empty())
	//	gen_coupled_reactions();
	//if (GR.empty())
	//	group_reactions();
	//if ((int)R.size() != N.cols()){
	//	cout << "Vector with reversibility is not consistent." << endl;
	//	return;
	//}
	int size;
	int g = GR.size();
	Rred.resize(g - 1);
	for (int i = 1; i < g; i++){
		// cout << GR[i][0] << endl;                  //print para teste
		size = GR[i].size();
		Rred[i - 1] = true;
		for (int j = 0; j < size; j++){
			if (!R[GR[i][j]]){
				Rred[i - 1] = false;
				break;
			}
		}
	}
// 	for (int i = 0; i < (int)R.size(); i++)         //print para teste
// 		cout << R[i] << endl;                         //print para teste
// 	cout << endl;   
// 	cout << "here" << endl;            
// 	for (int i = 0; i < (int)Rred.size(); i++)      //print para teste
// 		cout << Rred[i] << endl;                      //print para teste
// 	cout << "here" << endl;
// }

/*
checks if the new combined row is not equal or equivalent to the rows already in the tableau

@param begin an iterator corresponding to the first entry of the tableau that v will be compared to
       end an iterator corresponding to the location right after the last entry of the tableau that v will be compared to
       v the new combined row
@return false if the row is not equal or equivalent to a row already in the temporary tableau
        true if the row must not be added
*/
bool StoichAnalysis::check_condition(const Eigen::VectorXd *begin, const Eigen::VectorXd *end, const Eigen::VectorXd &v){
	for (auto it = begin; it != end; it++){
		//cout << v.segment(cols, rows).transpose() << endl;
		if (condition((*it).segment(cols, rows), v.segment(cols, rows))){
			return true;
		}
	}
	return false;
}

/*
checks if the new combined row is not equal or equivalent to the rows already in the tableau

@param *a pointer to the fisrt element of the temporary tableau
size size of the temporary tableau
v the new combined row
@return false if the row is not equal or equivalent to a row already in the temporary tableau
true if the row must not be added
*/
//bool StoichAnalysis::check_condition(Eigen::VectorXd *a, const size_t size, Eigen::VectorXd &v){
//	for (int i = 0; i < size; i++){
//		if (condition(a[i].segment(cols, rows), v.segment(cols, rows))){
//			return true;
//		}
//	}
//	return false;
//}

/*
checks if the new combined row has zeros in different positions than the rows already in the
temporary tableau

@param t pointer to a row already in the temporary tableau
       v pointer to the new combined row
@return true if v has zero in the same positions as another existing row
        false if v has at least one zero in a new position

updated: 04/07/2018
*/
bool StoichAnalysis::condition(const Eigen::VectorXd &t, const Eigen::VectorXd &v) const{
	int size = v.size();
	for (int j = 0; j < size; j++){
		if (v(j) > -tol && v(j) < tol){
			if (t(j) < -tol || t(j) > tol){
				return false;
			}
		}
	}
	return true;
}

/*
Updates the original tableau T for conservation relations calculation.

updated: 04/07/2018
*/
void StoichAnalysis::updateT(Tableau &Tt){
	T = Tt;
	//printT(); cout << endl << endl; // print para teste
	Tt.free();
}

/*
Transforms the tableau from a std vector of Eigen vectors in a Eigen Matrix

@return the calculated tableau as a Eigen matrix
*/
Eigen::MatrixXd StoichAnalysis::augTableauToMatrix(){
	int size = T.size();
	int a = rows;
	int b = cols;
	Eigen::MatrixXd Tb(size, a);
	for (int i = 0; i < size; i++){
		//Tb.row(i) = T[i].segment(b,a).transpose();
		Tb.row(i) = (*(T.begin() + i)).segment(b, a).transpose();
	}
	return Tb;
}

/*
Returns the calculated tableau with all reactions (complete N)

@return the right side of T that corresponds to the calculated properties 
*/
Eigen::MatrixXd StoichAnalysis::getT(){        
	Eigen::MatrixXd Tb, M;
	int t = T.size();
	int r = N.cols();
	Tb = augTableauToMatrix();
	M.resize(t, r);
	GR.size();
	int size = GR[0].size();
	for (int j = 0; j < size; j++)
		M.col(GR[0][j]) = Eigen::VectorXd::Zero(t).transpose();
	for (int i = 0; i < rows; i++){
		size = GR[TO[i]+1].size();
		for (int j = 0; j < size; j++){
			M.col(GR[TO[i] + 1][j]) = Tb.col(TO[i]) * (CRcoeff[GR[TO[i] + 1][j]] / CRcoeff[GR[TO[i] + 1][0]]);              
		}
	}
	return M;
}

/*
Generates the tableau by augmenting the stoichiometry matrix with the identity matrix

@param S Eigen matrix to be augmented
*/
void StoichAnalysis::genTableau(const Eigen::MatrixXd &S){
	Eigen::VectorXd v(S.rows() + S.cols());
	for (size_t i = 0; i < S.rows(); i++){
		v << S.row(i).transpose(), Eigen::VectorXd::Zero(S.rows());
		v(S.cols() + i) = 1;
		addRowToTableau(v, i);
		//T.push_back(v);
		//T[i](S.cols() + i) = 1;
	}
}

/*

*/
void StoichAnalysis::addRowToTableau(const Eigen::VectorXd &v, const size_t i){
	T.push_back(v);
}

void StoichAnalysis::clear(){
	delete this;
}

int StoichAnalysis::get_rev_efms() const{
	int n = 0;
	for (auto it = FluxR.begin(); it != FluxR.end(); it++){
		if (*it) n++;
	}
	return n;
}

// Functions for Positive Conservation Relations

/*
Calculates the Positive Conservation Relations of a metabolic network 
*/
void PosConsRel::run(){
	genTableau(N);
	this->StoichAnalysis::run();
}

/*
Checks if two rows of the tableau can be combined for conservation relations calculation 
and do so if possible.

@param pos column of tableau that is being analysed
       row1 first row to be combined
	   row2 second row to be combined
	   begin an iterator that points to the begining of the tableau
	   end an iterator that points to the last row in the tableau to be compared
	   back an insert iterator that appends a new vector at the end of the tableau
*/
Eigen::VectorXd PosConsRel::combineRows(const int pos, const Eigen::VectorXd *row1, const Eigen::VectorXd *row2){
	Eigen::VectorXd v;
	if ((*row1)(pos) * (*row2)(pos) < -tol){
		v = (*row1) + (*row2) * abs((*row1)(pos) / (*row2)(pos));
		//if (Tt.empty() || !check_condition(&Tt[0], Tt.size(), v))
		//if (begin == end || !check_condition(begin, end, v))
			//*inserter = v;
	}
	return v;
}

void PosConsRel::addRowToEnd(const Eigen::VectorXd &row, Tableau &Tt, const Eigen::VectorXd *row1, const Eigen::VectorXd *row2){
	Tt.push_back(row);
}

/*
Add a new row to the temporary tableau for conservation relations calculation

@param row row to be added
Tt temporary tableau that stores the rows to be added to the next tableau
*/
void PosConsRel::addRow(const Eigen::VectorXd *row, Tableau &Tt){
	Tt.push_back((*row));
}

/*
Returns the positive conservation relations of N calculated by Schuster's method

@return an Eigen matrix with the right side of the calculated tableau
*/
Eigen::MatrixXd PosConsRel::getT(){
	return augTableauToMatrix();
}

// Functions for Elementary Flux Modes (EFM)

/*
Creates an EFM object

@param S stoichiometry matrix of internal metabolites 
       r a vector indicating reversibility of reactions in the same order they appear in S
	     true if reversible and false if irreversible
*/
EFM::EFM(Eigen::MatrixXd S, deque<bool> r){
	this->StoichAnalysis::setNR(S, r);
	reduceN();
	reduceR();
}

/*
Generates tableau for calculating EFM using reduced N and R

@param S the stoichiometry matrix that will generate the tableau
*/
void EFM::genTableau(const Eigen::MatrixXd &S){
	rows = S.cols();
	cols = S.rows();
	this->StoichAnalysis::genTableau(S.transpose());
	//int a = rows;
	//rows = cols; cols = a;
	//Eigen::VectorXd v(rows + cols);
	//for (int i = 0; i < rows; i++){
	//	v << Nred.col(i), Eigen::VectorXd::Zero(rows);
	//	v(cols + i) = 1;
	//	addRowToTableau(v, i);
	//	//if (!Rred[i]){
	//	//	T.push_back(v);
	//	//	TO.push_back(i);
	//	//}
	//	//else{
	//	//	T.insert(T.begin(), v);
	//	//	TO.insert(TO.begin(), i);
	//	//}
	//}
	//cout << endl;                                      //print para teste
	//for (int i = 0; i < rows; i++){                    //print para teste
	//	//T[i](cols + i) = 1;                              
	//	cout << TO[i] << endl;                         //print para teste
	//}
	//rev = count(Rred.begin(), Rred.end(), true);
	FluxR = Rred;
	//revt = rev;              
	//revtt = 0;
	//cout << endl << Nred << endl;                      //print para teste
	//cout << endl;                                      //print para teste
	//for (int i = 0; i < (int)T.size(); i++)            //print para teste
	//	cout << T[i].transpose() << endl;              //print para teste
	//cout << endl;                                      //print para teste
}

/*

*/
void EFM::addRowToTableau(const Eigen::VectorXd &v, const size_t i){
	//if (!Rred[i]){
		T.push_back(v);
		TO.push_back(i);
	//}
	//else{
	//	T.insert(T.begin(), v);
	//	TO.insert(TO.begin(), i);
	//}
}

/*
Calculates the EFM of a metabolic network using Schuster's method
*/
void EFM::run(){
	genTableau(Nred);
	this->StoichAnalysis::run();
	final_check();
	changeValues();                                //para melhor vizualizacao
}

/*
Add a new row to the temporary tableau for EFM calculation

@param row row to be added
       Tt temporary tableau that stores the rows to be added to the next tableau
*/
void EFM::addRow(const Eigen::VectorXd *row, Tableau &Tt){
	Tt.push_back(*row);
	FluxRt.push_back(FluxR[row - T.begin()]);
	/*if (row < revt){
		Tt.insert(Tt.begin(), T[row]);
		revtt++;
	}
	else
		Tt.push_back(T[row]);*/
}

/*
Checks if two rows of the tableau can be combined for EFM calculation
and do so if possible.

@param pos column of tableau that is being analysed
       row1 first row to be combined
       row2 second row to be combined
       begin an iterator that points to the begining of the tableau
	   end an iterator that points to the last row in the tableau to be compared
	   back an insert iterator that appends a new vector at the end of the tableau
*/
Eigen::VectorXd EFM::combineRows(const int pos, const Eigen::VectorXd *row1, const Eigen::VectorXd *row2){
	Eigen::VectorXd v;
	short int a = 0;
	short int b = 1;
	if ((*row1)(pos) * (*row2)(pos) < -tol){
		a = 1;
	}
	else if ((*row1)(pos) * (*row2)(pos) > tol){
		if (FluxR[row1 - T.begin()])                                                     //ARRUMAR
			a = -1;
		else if (FluxR[row2 - T.begin()]){
			a = 1; b = -1;
		}
	}
	if (a != 0){
		v = a * (*row1) + b * (*row2) * abs((*row1)(pos) / (*row2)(pos));
		//if (Tt.empty() || !check_condition(&Tt[0], Tt.size(), v)){
		//if (begin == end || !check_condition(begin, end, v)){
		//	if ((FluxR[row1 - T.begin]) && (FluxR[row2 - T.begin])){
		//		*front = v;                                                  //TO DO: ARRUMAR AQUI
		//		revtt++;
		//	}
		//	else
		//		*back = v;
		//}
	}
	return v;
}

void EFM::addRowToEnd(const Eigen::VectorXd &row, Tableau &Tt, const Eigen::VectorXd *row1, const Eigen::VectorXd *row2){
	Tt.push_back(row);
	//int n = row1 - T.begin();
	FluxRt.push_back((FluxR[row1 - T.begin()]) && (FluxR[row2 - T.begin()]));

}

/*
Updates the original tableau T for EFM calculation.
*/
void EFM::updateT(Tableau &Tt){
	FluxR.clear();
	FluxR = FluxRt;
	FluxRt.clear();
	//revt = revtt;
	//revtt = 0;
	this->StoichAnalysis::updateT(Tt);
}

/*
Checks in the last tableau if a row added to it is simpler than a row added priorly
*/
void EFM::final_check(){
	size_t size = T.size();
	vector<size_t> del;
	for (auto it = T.begin(); it != T.end() - 1; it++){
		if (check_condition(it + 1, T.end(), *it))
			del.push_back(it - T.begin());
	}
	//for (size_t i = 0; i < size - 1; i++){
	//	if (check_condition(&T[i + 1], size - i - 1, T[i]))
	//		del.push_back(i);
	//}
	size = del.size();
	T.erase(del);
	for (int i = size - 1; i >= 0; i--){
		//T.erase(T.begin() + del[i]);
		FluxR.erase(FluxR.begin() + del[i]);
		//if (del[i] < revt)
		//	revt--;
	}
}

// Funcoes da classe EFMNS

/*
Creates an EFMNS object

@param S stoichiometry matrix of internal metabolites
       r a vector indicating reversibility of reactions in the same order they appear in S
       true if reversible and false if irreversible
*/
EFMNS::EFMNS(Eigen::MatrixXd S, deque<bool> r){
	this->StoichAnalysis::setNR(S, r);
}

/*
Transform a Eigen matrix into a std vector of Eigen Vectors using the matrix rows

@param M the matrix to be transformed
*/
void EFMNS::matrixToTableau(const Eigen::MatrixXd &M){
	for (size_t i = 0; i < M.rows(); i++){
		FluxR.push_back(checkFluxRev(M.row(i).transpose()));
		T.push_back(M.row(i).transpose());
	//	cout << FluxR[i] << endl;                             //print para teste
	}
	//cout << endl;                                             //print para teste
}

/*
Generates the tableau for calculation EFM using the nullspace approach

@param S the stoichiometry matrix considering only internal metabolites
*/
void EFMNS::genTableau(const Eigen::MatrixXd &S){
	reduceN();
	reduceR();
	rows = Nred.cols();
	cols = Nred.rows(); 
	//cout << Nred << endl;                                                                //print para teste
	Eigen::SparseMatrix<double> M = Nred.fullPivLu().kernel().transpose().sparseView(1, tol);
	// cout << endl << M << endl;                                                           //print para teste
	Eigen::MatrixXd Mt(M.rows(), M.cols());
	int id = 0, in = 0;
	for (size_t i = 0; i < rows ; i++){
		if (M.col(i).nonZeros() == 1){
			//Mt.col(rows - 1 - id) = M.col(i).toDense();
			Mt.col(rows - 1 - id) = M.col(i);
			//T.insert(T.begin(), M.col(i));
			TO.insert(TO.end() - id, i);
			id++;
		}
		else{
			//Mt.col(i - id) = M.col(i).toDense();
			Mt.col(i - id) = M.col(i);
			//T.push_back(M.col(i));
			TO.insert(TO.begin() + in, i);
			in++;
		}
	}
	matrixToTableau(Mt);
	//for (auto it = T.begin(); it != T.end(); it++){    //print para teste
	//	cout << it->transpose() << endl;               //print para teste
	//}                                                  //print para teste
	//cout << endl;                                      //print para teste
}

/*
Checks if the given flux is reversible by verifying the reversibility of each 
reaction envolved

@param v the flux that will be check for reversibility
*/
bool EFMNS::checkFluxRev(const Eigen::VectorXd &v){
	for (size_t i = 0; i < v.size(); i++){
		if ((v(i) > tol || v(i) < -tol) && !Rred[TO[i]]){
			return false;
		}
	}
	return true;
}

/*
Calculates the EFM of a metabolic network using the nullspace approach
*/
void EFMNS::run(){
	genTableau(Nred);
	Eigen::VectorXd v;
	size_t size;
	for (size_t i = 0; i < rows; i++){
		size = T.size();
		for (auto j = T.begin(); j != T.end() - 1; j++){
			for (auto k = j + 1; k != T.end(); k++){
		//for (size_t j = 0; j < size - 1; j++){
			//for (size_t k = j + 1; k < size; k++){
				v = combineRows(i, j, k);
			}
		}
	}
	final_check();
	changeValues();
//	for (auto it = T.begin(); it != T.end(); it++){     //print para teste
//		cout << it->transpose() << endl;                //print para teste
//	}                                                   //print para teste
//	cout << endl;                                       //print para teste
}

/*
Checks if two rows of the tableau can be combined for conservation relations calculation
and do so if possible.

@param pos column of tableau that is being analysed
       row1 first row to be combined
       row2 second row to be combined
       begin an iterator that points to the begining of the tableau
	   end an iterator that points to the last row in the tableau to be compared
	   back an insert iterator that appends a new vector at the end of the tableau
*/
Eigen::VectorXd EFMNS::combineRows(const int pos, const Eigen::VectorXd *row1, const Eigen::VectorXd *row2){
	Eigen::VectorXd v;
	int a = 0, b = 1;
	if ((*row1)(pos) * (*row2)(pos) < -tol){
		a = 1;
	}
	else if ((*row1)(pos) * (*row2)(pos) > tol){
		//if (FluxR[row1]) a = -1;
		//else if (FluxR[row2]) b = -1;
		if (FluxR[row1 - T.begin()]) a = -1;
		else if (FluxR[row2 - T.begin()]) b = -1;
	}
	if (a != 0){
		v = a * (*row1) + b * (*row2) * abs((*row1)(pos) / (*row2)(pos));
		if (!check_condition(T.begin(), T.end(), v)){
			T.push_back(v);
			FluxR.push_back((FluxR[row1 - T.begin()] && FluxR[row2 - T.begin()]));
		}
	}
	return v;
}

/*
checks if the new combined row is not equal or equivalent to the rows already in the tableau

@param begin an iterator corresponding to the first entry of the tableau that v will be compared to
       end an iterator corresponding to the location right after the last entry of the tableau that v will be compared to
       v the new combined row
@return false if the row is not equal or equivalent to a row already in the temporary tableau
        true if the row must not be added
*/
bool EFMNS::check_condition(const Eigen::VectorXd *begin, const Eigen::VectorXd *end, const Eigen::VectorXd &v) const{
	for (auto it = begin; it != end; it++){
		if (condition(*it, v)){
			return true;
		}
	}
	return false;
}

/*
Removes any flux mode with negative irreversible reactions 
and checks in the last tableau if a row added to it is simpler than a row added priorly
*/
void EFMNS::final_check(){
	size_t size = T.size();
	vector<size_t> del;
	for (auto it = T.begin(); it != T.end(); it++){
		if (checkReacRev(*it) || check_condition(it + 1, T.end(), *it))
			del.push_back(it - T.begin());
	}
	size = del.size();
	for (int i = size; i != 0; i--){
		//cout << del[i - 1] << endl;                               //print para teste
		T.erase(del);
	}
	//cout << endl;                                             //print para teste
}

/*
Checks if a flux mode has a reaction with negative flux

@param v a vector representing a flux mode
@return true if there is a irreversible reaction has a negative flux
        false otherwise
*/
bool EFMNS::checkReacRev(const Eigen::VectorXd &v){
	for (size_t i = 0; i < v.size(); i++){
		if (v(i) < -tol && !Rred[TO[i]]){
			return true;
		}
	}
	return false;
}

/*
Transforms the tableau from a std vector of Eigen vectors in a Eigen Matrix

@return the calculated tableau as a Eigen matrix
*/
Eigen::MatrixXd EFMNS::augTableauToMatrix(){
	int size = T.size();
	int sizeb = (*T.begin()).size();
	Eigen::MatrixXd Tb(size, sizeb);
	for (int i = 0; i < size; i++){
		Tb.row(i) = (*(T.begin() + i)).transpose();
	}
	return Tb;
}


