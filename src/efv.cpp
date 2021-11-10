#include "stdafx.h"
#include "efv.h"

using namespace std;
using namespace efv;

EFV::EFV(Eigen::MatrixXd S, Eigen::VectorXd v, int mes){
	gen_matrix(S, v, mes);
}

void EFV::gen_matrix(Eigen::MatrixXd S, Eigen::VectorXd v, int mes){
	cout << endl << v << endl;
	if (v.size() == mes){
		//v = v / abs(v(0));
		M = Eigen::MatrixXd::Zero(S.rows() + 2 * mes, S.cols() + 1 + 2 * mes);
		M.block(S.rows(), S.cols() + 2 * mes, mes, 1) = - v;
		M.block(S.rows() + mes, S.cols() + 2 * mes, mes, 1) = v;
		M.block(0, 0, S.rows(), S.cols()) = S;
		M.block(S.rows(), S.cols() - mes, mes, mes) = Eigen::MatrixXd::Identity(mes, mes);
		M.block(S.rows() + mes, S.cols() - mes, mes, mes) = Eigen::MatrixXd::Identity(mes, mes) * -1;
		M.block(S.rows(), S.cols(), 2 * mes, 2 * mes) = Eigen::MatrixXd::Identity(2 * mes, 2 * mes);
	}


}