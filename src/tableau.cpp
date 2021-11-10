#include "stdafx.h"
#include "tableau.h"

using namespace std;

Tableau::~Tableau(){
	free();
}

Tableau::Tableau(const Tableau &s){
	auto newdata = alloc_n_copy(s.begin(), s.end());
	elements = newdata.first;
	first_free = cap = newdata.second;
}

Tableau &Tableau::operator=(const Tableau &rhs){
	free();
	auto data = alloc_n_copy(rhs.begin(), rhs.end());
	elements = data.first;
	first_free = cap = data.second;
	return *this;
}

void Tableau::free(){
	if (elements){
		for (auto p = first_free; p != elements;){
			alloc.destroy(--p);                            //primeiro volta um e depois destroi o objeto que o ponteiro aponta
		}
		alloc.deallocate(elements, cap - elements);
	}
	elements = first_free = cap = nullptr;
}

pair<Eigen::VectorXd*, Eigen::VectorXd*> Tableau::alloc_n_copy(const Eigen::VectorXd *v1, const Eigen::VectorXd *v2){
	auto data = alloc.allocate(v2 - v1);
	return{ data, uninitialized_copy(v1, v2, data) };   
}

void Tableau::reallocate(){
	auto newcapacity = size() ? 2 * size() : 1;
	auto newdata = alloc.allocate(newcapacity);
	auto dest = newdata;
	auto elem = elements;
	for (size_t i = 0; i != size(); ++i)
		alloc.construct(dest++, std::move(*elem++));
	free();
	elements = newdata;
	first_free = dest;
	cap = elements + newcapacity;
}

void Tableau::allocate(size_t n){
	if (!elements){
		auto newdata = alloc.allocate(n);
		elements = first_free = newdata;
		cap = newdata + n;
	}
}

void Tableau::push_back(const Eigen::VectorXd v){
	check_n_alloc();
	alloc.construct(first_free++, v);
}

void Tableau::erase(const vector<size_t> &del){
	auto data = alloc.allocate(first_free - elements - del.size());
	auto dest = data;
	auto elem = elements;
	auto j = 0;
	for (size_t i = 0; i != size(); i++){
		if (j >= del.size() || i != del[j])
			alloc.construct(dest++, std::move(*elem++));
		else{
			elem++;
			j++;
		}
	}
	free();
	elements = data;
	first_free = cap = dest;
}