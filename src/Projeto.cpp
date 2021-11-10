// Projeto.cpp : Defines the entry point for the console application.
//

#include "network.h"
#include "metabolite.h"
#include "reaction.h"
#include "stoichiometry.h"
#include "interface.h"
// #include "OsiClpSolverInterface.hpp"
// #include "CoinPackedMatrix.hpp"
// #include "CoinPackedVector.hpp"
#include <fstream>
#include <iterator>
#include <Eigen/Dense>
#include <Eigen/src/Core/IO.h>
//#include <iostream>
//#include <string>
//#include <fstream>

using namespace std;
using namespace stoichiometry;
using namespace network;


int main()
{
	/****************Teste com o exemplo 8.3 da livro do Stephanopoulos*************/

	//Eigen::MatrixXd N(18, 18);
	//deque<bool> r(18);

	//r = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	//N << 1, -0.5, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 1, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0,
	//	0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1,
	//	0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 1, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 1, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 1, -1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0,
	//	-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
	//
	//Interface candida(N, false, r);
	//
	//cout << candida.calculate("EFM") << endl;
	//cout << endl;
	//cout << candida.calculate("CONSREL") << endl;
	//cout << endl;
	//cout << candida.calculate("EFMNS") << endl;
	//cout << endl;

	//candida.analyseEFM();

	//Eigen::VectorXd V(7);
	//V << 1, 0.3, 0, 0.2, 0.1, 0.2, 0.2;

	/****************REDE DA B SACCHARI COM FLUXOS DA TESE DA TATIANE*************/

	Interface input("CASESTUDY2.dat");
	// Interface input("ECOLI2.dat");

	// Eigen::MatrixXd EFM(305, 46);

	input.generateOutputFile("CASESTUDY2_out");
	// input.generateOutputFile("ECOLI2out");
	/*input.generateOutputFile("ESTUDOout", EFM);

	respbiophbout << respbiophb.calculate("EFMNS") << endl;
	cout << "This network has " << respbiophb.getStoichResult("EFMNS").rows() << "EFMs" << endl;
	respbiophbout << respbiophb.getStoichResult("EFMNS") << endl;*/

	/******************REDE DA E COLI OBTIDA DO ARTIGO DE LEIGHTY 2013**********************/

	//Interface ecoli("TESTE.dat");

	//const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ",", "\n");

	//ofstream output("TESTE-N.csv");
	//output << ecoli.getStoichMatrix().format(CSVFormat);

	//ofstream output2("TESTE-Nelst.csv");
	//output2 << ecoli.getElstStoichMatrix().format(CSVFormat);
	//output2.close();
	//cout << endl << ecoli.getElstStoichMatrix().rows() << endl << ecoli.getElstStoichMatrix().cols() << endl;

	//ofstream output3("TESTE-Met.txt");
	//ecoli.generateListMetabolites(output3);
	//output3.close();

	//ofstream output4("TESTE-Reac.txt");
	//ecoli.generateListReactions(output4);
	//output4.close();

	//ofstream output5("TESTE-Nsubs.csv");
	//output5 << ecoli.getElstSubsStoichMatrix().cwiseAbs().format(CSVFormat);
	//output5.close();

	//ofstream output6("TESTE-Nrev.csv");
	//output6 << ecoli.getRevStoichMatrix().format(CSVFormat);
	

	/****************** alpha spectrum ***********************/

	//// Create a problem pointer. We use the base class here.
	//OsiSolverInterface *si;
	//// When we instantiate the object, we need a specific derived class.
	//si = new OsiClpSolverInterface;
	//// Read in an mps file. This one�s from the MIPLIB library.
	//si->readMps("EcoliRec");
	//// Solve the (relaxation of the) problem
	//si->initialSolve();
	//// Check the solution
	//if (si->isProvenOptimal())
	//{
	//	std::cout << "Found optimal solution!" << std::endl;
	//	std::cout << "Objective value is " << si->getObjValue() << std::endl;
	//	int n = si->getNumCols();
	//	const double *solution;
	//	solution = si->getColSolution();
	//	// We could then print the solution or examine it.
	//}
	//else
	//{
	//	std::cout << "Did not find optimal solution." << std::endl;
	//	// Could then check other status functions.
	//}


	return 0;
}

