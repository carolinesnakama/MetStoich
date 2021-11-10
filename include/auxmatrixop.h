#include <Eigen/Dense>
#include <vector>

#ifndef AUXMATRIXOP_H
#define AUXMATRIXOP_H

void remove_row(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
void remove_col(Eigen::MatrixXd& matrix, unsigned int colToRemove);
void remove_el(Eigen::VectorXi& vector, unsigned int elToRemove);
void remove_cols(Eigen::MatrixXd& M, std::vector<int>& ind);
void remove_rows(Eigen::MatrixXd& M, std::vector<int>& ind);
double minAbsCoeff(Eigen::VectorXd const & v);
std::vector<std::vector<std::size_t>> group_same_rows(Eigen::MatrixXd &M);

#endif