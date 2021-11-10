#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include <numeric>
#include <deque>
#include <utility> 
#include "auxmatrixop.h"

#ifndef EFMANALYSIS_H
#define EFMANALYSIS_H

namespace efmanalysis{

	class EFMAnalysis{
	private:
		Eigen::MatrixXd G, Sext, EFM;
		std::vector<std::size_t> efm_rank;
		std::vector<std::size_t> cycle_efm;
		std::vector<std::vector<std::size_t>> efm_families;
		std::vector<std::vector<std::pair<std::size_t, double>>> group;
		void gen_cycle_efm();
		void gen_efm_families();
		void gen_efm_rank(std::vector<int> ids);
		void group_efm_by_prod(std::vector<int> ids);
	public:
		EFMAnalysis() = delete;
		~EFMAnalysis(){};
		EFMAnalysis(Eigen::MatrixXd efm, Eigen::MatrixXd sext, std::vector<int> ids_upt);
		std::vector<std::size_t> get_cycle_efm();
		std::vector<std::vector<std::size_t>> get_efm_families();
		std::vector<std::size_t> get_efms_reaction();
		std::vector<std::size_t> get_efm_rank(std::vector<int> ids);
		Eigen::MatrixXd get_unique_efm(const int mes_reac);
		Eigen::MatrixXd get_efm_no_cycle();
		Eigen::MatrixXd get_ranked_efm(std::vector<int> ids);
		Eigen::MatrixXd get_norm_efm(){ return EFM; };
		std::vector<double> get_max_contribution(Eigen::VectorXd v, bool families = false);
	};

}

#endif