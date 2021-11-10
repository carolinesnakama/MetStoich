#include "stdafx.h"
#include "auxmatrixop.h"

using namespace std;

void remove_row(Eigen::MatrixXd& matrix, unsigned int rowToRemove){
	unsigned int numRows = matrix.rows() - 1;
	unsigned int numCols = matrix.cols();
	if (rowToRemove < numRows)
		matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);
	matrix.conservativeResize(numRows, numCols);
}

void remove_col(Eigen::MatrixXd& matrix, unsigned int colToRemove){
	unsigned int numRows = matrix.rows();
	unsigned int numCols = matrix.cols() - 1;
	if (colToRemove < numCols)
		matrix.block(0, colToRemove, numRows, numCols - colToRemove) = matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);
	matrix.conservativeResize(numRows, numCols);
}

void remove_el(Eigen::VectorXi& vector, unsigned int elToRemove){
	unsigned int numEls = vector.size() - 1;
	if (elToRemove < numEls)
		vector.segment(elToRemove, numEls - elToRemove) = vector.segment(elToRemove + 1, numEls - elToRemove);
	vector.conservativeResize(numEls);
}

double minAbsCoeff(Eigen::VectorXd const & v){   //menor coeficiente diferente de zero valor absoluto TERMINAR
	Eigen::VectorXd vector;
	vector = v.cwiseAbs();
	double a = vector(0);
	int size = vector.size();
	for (int i = 1; i < size; i++){
		if ((vector(i) < a || a < 1e-12) && vector(i) != 0){
			a = vector(i);
		}
	}
	return a;
}

void remove_cols(Eigen::MatrixXd& M, std::vector<int>& ind){
	int size = ind.size();
	for (int i = size - 1; i >= 0; i--)
		remove_col(M, ind[i]);
}

void remove_rows(Eigen::MatrixXd& M, std::vector<int>& ind){
	int size = ind.size();
	for (int i = size - 1; i >= 0; i--)
		remove_row(M, ind[i]);
}

vector<vector<size_t>> group_same_rows(Eigen::MatrixXd &M){
	vector<vector<size_t>> groups;
	int size = M.rows();
	bool added = false;
	groups.push_back({});
	groups[0].push_back(0);
	for (unsigned int i = 1; i < size; i++){
		for (unsigned int j = 0; j < groups.size(); j++){
			if (M.row(i).isApprox(M.row(groups[j][0]), 1e-8)){
				groups[j].push_back(i);
				added = true;
				break;
			}
		}
		if (added)
			added = false;
		else{
			groups.push_back({ i });
		}
	}
	return groups;
}
