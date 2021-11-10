#define _SCL_SECURE_NO_WARNINGS

#include <vector>
#include <memory>
#include <Eigen/Dense>

#ifndef TABLEAU_H
#define TABLEAU_H

class Tableau{
private:
	std::allocator<Eigen::VectorXd> alloc;
	Eigen::VectorXd *elements;                          //first element fo the allocated memory
	Eigen::VectorXd *first_free;                        //past the last constructed element
	Eigen::VectorXd *cap;                               //past the last element of the allocated memory
	void check_n_alloc(){ if (size() == capacity()) reallocate(); };
	void reallocate();
	std::pair<Eigen::VectorXd*, Eigen::VectorXd*> alloc_n_copy(const Eigen::VectorXd*, const Eigen::VectorXd*);
public:
	Tableau() : elements(nullptr), first_free(nullptr), cap(nullptr){};
	Tableau(const Tableau &);
	Tableau &operator=(const Tableau&);
	~Tableau();
	Eigen::VectorXd * begin() const { return elements; };
	Eigen::VectorXd * end() const { return first_free; };
	std::size_t size() const { return first_free - elements; };
	std::size_t capacity() const { return cap - elements; };
	void erase(const std::vector<std::size_t> &);
	void allocate(std::size_t);
	void free();
	void push_back(const Eigen::VectorXd);
};

// std::allocator<Eigen::VectorXd> Tableau::alloc;     COLOCAR NO ARQUIVO .CPP

#endif