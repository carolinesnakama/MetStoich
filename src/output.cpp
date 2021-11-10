// #include "stdafx.h"
//#include "output.h"
//
//using namespace std;
//
//void output::create_file(string name){
//	file.open((name + ".dat"), std::ios::binary);
//	file << name << endl << endl;
//}
//
//void output::write_text(string text){
//	file << endl << text << endl;;
//}
//
//void output::write_matrix(const Eigen::MatrixXd& M){
//	file << M << endl;
//}
//
//template<typename T>
//void output::write_vector(const T& v){
//	copy(v.begin(), v.end(), ostream_iterator<typename T::value_type>(file, "\n"));
//}
//
