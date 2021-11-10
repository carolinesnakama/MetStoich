//#include <fstream>
//#include <iostream>
//#include <vector>
//#include <string>
//#include <Eigen/Dense>
//#include "network.h"
//#include "stoichanalysis.h"
//
//#ifndef OUTPUT_H
//#define OUTPUT_H
//
//template <typename C>
//class output
//{
//private:
//	std::ofstream file;
//	network::Network net;
//	C stoich;
//public:
//	output(){};
//	output(std::string name){ create_file(name); };
//	void set_network(network::Network rede){ net = rede; };
//	void set_calculation(C op){ calc = op; };
//	void create_file(std::string name);
//	void write_text(std::string text);
//	void write_matrix(const Eigen::MatrixXd& M);
//	template<typename T> void write_vector(std::vector<T> const& v);
//	void close_file(){ file.close(); };
//};
//
//#endif 