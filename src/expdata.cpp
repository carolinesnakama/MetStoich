#include "stdafx.h"
#include "expdata.h"

using namespace std;
using namespace expdata;

void ExpData::setExtReactions(Eigen::MatrixXd S){
	for (size_t i = 0; i < S.cols(); i++){
		if (!(S.col(i).array() > 0.0).any() || !(S.row(i).array() < 0.0).any())
			p.push_back(i);
	}
}

bool ExpData::setV(Eigen::MatrixXd S, Eigen::VectorXd v){
	setExtReactions(S);
	if (v.size() == p.size()){
		V = v;
		return true;
	}
	return false;
}

void ExpData::set_var(Eigen::MatrixXd W){
	if (!W.isZero())
		Var = W;
	else{
		Var = V.asDiagonal() * 0.05;
	}
}

