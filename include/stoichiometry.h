#include <iostream>
#include <vector>
#include <deque>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <time.h>
#include "auxmatrixop.h"
#include "tableau.h"
//#include "efmanalysis.h"

#ifndef STOICHIOMETRY_H
#define STOICHIOMETRY_H

namespace stoichiometry{
	class StoichAnalysis{
	private:
		virtual void addRowToTableau(const Eigen::VectorXd &v, const size_t i);
		virtual void addRow(const Eigen::VectorXd *row, Tableau &Tt){};
		virtual void addRowToEnd(const Eigen::VectorXd &row, Tableau &Tt, const Eigen::VectorXd *row1 = nullptr, const Eigen::VectorXd *row2 = nullptr){};
		virtual Eigen::VectorXd combineRows(const int pos, const Eigen::VectorXd *row1, const Eigen::VectorXd *row2){ return *row1; };
		virtual void final_check(){};
	protected:
		const double tol = 1e-10;
		size_t rows, cols;     //, rev, revt, revtt;
		Eigen::MatrixXd N, K, Nred;
		//std::vector<Eigen::VectorXd> T;
		Tableau T;
		std::vector<int> CR, TO;
		std::vector<std::vector<int>> GR;
		std::vector<double> CRcoeff;
		std::deque<bool> FluxR;
		std::deque<bool> R, Rred;
		virtual bool check_condition(const Eigen::VectorXd *begin, const Eigen::VectorXd *end, const Eigen::VectorXd &v);
		bool condition(const Eigen::VectorXd &t, const Eigen::VectorXd &v) const;
		void gen_coupled_reactions();
		void changeValues();
		void init(Eigen::MatrixXd S){ N = S; rows = N.rows(); cols = N.cols(); };
		virtual void reduceN();
		virtual void reduceR();
		virtual void updateT(Tableau &Tt);
		virtual void genTableau(const Eigen::MatrixXd &S);
		std::vector<int> group_reactions();
		std::vector<int> get_blocked_reactions(const bool remove_rows = false);     //sendo N a matrix estequimetrica de metabolitos internos
		std::vector<int> get_coupled_reactions();
	public:
		virtual void run();
		void clear();
		void setNR(Eigen::MatrixXd S, std::deque<bool> r){ this->StoichAnalysis::init(S); R = r; };
		Eigen::MatrixXd getN() const { return N; };
		Eigen::MatrixXd getK() const { return K; };
		Eigen::MatrixXd getNred() const { return Nred; };
		virtual Eigen::MatrixXd getT();
		virtual Eigen::MatrixXd augTableauToMatrix();
		std::deque<bool> getR() const { return R; };
		std::deque<bool> getRred() const { return Rred; };
		std::vector<int> get_CR() const { return CR; };
		std::vector<std::vector<int>> get_GR() const { return GR; };
		int get_rev_efms() const;
	};

	class PosConsRel : public StoichAnalysis{
	private:
		virtual void addRow(const Eigen::VectorXd * row, Tableau &Tt);
		virtual Eigen::VectorXd combineRows(const int pos, const Eigen::VectorXd *row1, const Eigen::VectorXd *row2);
		void addRowToEnd(const Eigen::VectorXd &row, Tableau &Tt, const Eigen::VectorXd *row1, const Eigen::VectorXd *row2);
		virtual void reduceN(){};
		virtual void reduceR(){};
	public:
		PosConsRel(){};
		~PosConsRel(){};
		PosConsRel(Eigen::MatrixXd S){ this->StoichAnalysis::init(S); };
		virtual void run();
		virtual Eigen::MatrixXd getT(); 
	};

	class EFM : public StoichAnalysis{
	private:
		std::deque<bool> FluxRt;
		virtual void addRowToTableau(const Eigen::VectorXd &v, const size_t i);
		virtual void addRow(const Eigen::VectorXd * row, Tableau &Tt);
		virtual void addRowToEnd(const Eigen::VectorXd &row, Tableau &Tt, const Eigen::VectorXd *row1 = nullptr, const Eigen::VectorXd *row2 = nullptr);
		virtual Eigen::VectorXd combineRows(const int pos, const Eigen::VectorXd *row1, const Eigen::VectorXd *row2);
		virtual void updateT(Tableau &Tt);
		virtual void final_check();
		virtual void genTableau(const Eigen::MatrixXd &S);
	public:
		EFM(){};
		~EFM(){};
		EFM(Eigen::MatrixXd S, std::deque<bool> r);
		virtual void run();
	};

	class EFMNS : public StoichAnalysis{
	private:
		//std::deque<bool> FluxR;
		bool checkFluxRev(const Eigen::VectorXd &v);
		bool checkReacRev(const Eigen::VectorXd &v);
		virtual Eigen::MatrixXd augTableauToMatrix();
		virtual void genTableau(const Eigen::MatrixXd &S);
		virtual Eigen::VectorXd combineRows(const int pos, const Eigen::VectorXd *row1, const Eigen::VectorXd *row2);
		virtual bool check_condition(const Eigen::VectorXd *begin, const Eigen::VectorXd *end, const Eigen::VectorXd &v) const;
		virtual void final_check();
		void matrixToTableau(const Eigen::MatrixXd &M);
		Eigen::MatrixXd gauss();  //std::vector<Eigen::VectorXd> gauss();
	public:
		EFMNS(){};
		~EFMNS(){};
		EFMNS(Eigen::MatrixXd S, std::deque<bool> r);
		//void teste(){ genTableau(N); };                       //TESTE
		virtual void run();
	};

}

#endif