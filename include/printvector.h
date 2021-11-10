#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <fstream>

#ifndef PRINTVECTOR_H
#define PRINTVECTOR_H

template<typename T>
void print_vector(std::ofstream &file, const T& t) {
	std::copy(t.cbegin(), t.cend(), std::ostream_iterator<typename T::value_type>(file, ", "));
}


#endif