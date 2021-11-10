#include <Eigen/Dense>
#include <vector>
#include <iostream>

#ifndef EXPDATA_H
#define EXPDATA_H

namespace expdata{

	class ExpData{
	private:
		std::vector<std::size_t> p;
		Eigen::VectorXd V, Vrec;
		Eigen::MatrixXd Var;
		void reconcile(bool);
		void setExtReactions(Eigen::MatrixXd S);
	public:
		ExpData(){};
		~ExpData(){};
		std::vector<std::size_t> getExtReactions(){ return p; };
		Eigen::VectorXd getV(){ return V; };
		bool setV(Eigen::MatrixXd S, Eigen::VectorXd v);
		// void set_var(Eigen::MatrixXd W = Eigen::MatrixXd::Zero(1));
		void set_var(Eigen::MatrixXd W);
	};

}

#endif