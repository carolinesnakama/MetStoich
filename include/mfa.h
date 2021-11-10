#include <Eigen/Dense>

#ifndef MFA_H
#define MFA_H

namespace mfa{

	class MFA{
	protected:
		Eigen::MatrixXd N;
		Eigen::VectorXd V, X, e;
		bool exact_sol = false;
	public:
		MFA(){};
		MFA(Eigen::MatrixXd S, Eigen::VectorXd v);    //TODO: colocar verificacoes
		~MFA(){};
		void setN(Eigen::MatrixXd S);      //TODO: colocar verificacoes
		void setV(Eigen::VectorXd v);      //TODO: colocar verificacoes
		void run();
		void error();
		Eigen::VectorXd getX(){ return X; };
		Eigen::MatrixXd getN(){ return N; };
		Eigen::VectorXd getV(){ return V; };
		Eigen::VectorXd getError(){ return e; };
	};

}

#endif