#include "stdafx.h"
#include "mfa.h"

using namespace std;
using namespace mfa;

/*

*/
MFA::MFA(Eigen::MatrixXd S, Eigen::VectorXd v){
	N = S;
	V = v;
}

/*

*/
void MFA::setN(Eigen::MatrixXd S){
	N = S;
}

/*

*/
void MFA::setV(Eigen::VectorXd v){
	V = v;
}

/*

*/
void MFA::run(){
	Eigen::ColPivHouseholderQR<Eigen::MatrixXd> S;
	S = N.colPivHouseholderQr();
	X = S.solve(V);
	if (S.isInvertible())
		exact_sol = true;
}

/*

*/
void MFA::error(){
	if (exact_sol)
		e = Eigen::VectorXd::Zero(N.cols());
	else
		e = N*X - V;
}