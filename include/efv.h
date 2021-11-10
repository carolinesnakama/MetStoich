#include <Eigen/Dense>
#include <vector>
#include <iostream>

#ifndef EFV_H
#define EFV_H

namespace efv{

	class EFV{
	private:
		Eigen::MatrixXd M;
		void gen_matrix(Eigen::MatrixXd S, Eigen::VectorXd v, int mes);
	public:
		EFV(){};
		EFV(Eigen::MatrixXd S, Eigen::VectorXd v, int mes);
		Eigen::MatrixXd getM(){ return M; };
	};

}

#endif