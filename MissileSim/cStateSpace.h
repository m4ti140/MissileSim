#pragma once
#include <valarray>
#include <vector>
#include <stdexcept>
#include <memory>
#include <functional>
#include "utilities.h"


template<class Child>
class cStateSpace
{
public:
	void StepI() 
	{
		static_cast<Child*>(this)->Step();
	};
};

//using dynamic polymorphism to allow for different SSM engine to be selected with command line arguments, might do away with it if one approach is terribly slow


class cLTI : public cStateSpace<cLTI>
{
public:
	cLTI(
		std::shared_ptr<mtx> A,
		std::shared_ptr<mtx> B,
		std::shared_ptr<mtx> C,
		std::shared_ptr<mtx> D,
		std::shared_ptr<vec> u,
		std::shared_ptr<vec> y,
		vec x_0,
		std::shared_ptr<double> ts);
	~cLTI();

	void Step();
private:
	std::shared_ptr<mtx> A, B, C, D;
	std::shared_ptr<vec> u, y;
	vec x, dx;
	std::shared_ptr<double> ts;
	size_t dimArow, dimAcol, dimBrow, dimBcol, dimCrow, dimCcol, dimDrow, dimDcol, dimx, dimy, dimu;
};

class cGeneralSSM : public cStateSpace<cGeneralSSM>
{
public:
	cGeneralSSM(
		std::function<void(vec&, const vec&, const vec&)> f,
		std::function<void(vec&, const vec&, const vec&)> g,
		std::shared_ptr<vec> u,
		std::shared_ptr<vec> y,
		vec x_0,
		std::shared_ptr<double> ts);
	~cGeneralSSM();

	void Step();
private:
	std::function<void(vec&, const vec&, const vec&)> f;
	std::function<void(vec&, const vec&, const vec&)> g;
	std::shared_ptr<vec> u, y;
	vec x, dx, k1, k2, k3, k4;
	std::shared_ptr<double> ts;
};

inline cLTI::cLTI(
	std::shared_ptr<mtx> A,
	std::shared_ptr<mtx> B,
	std::shared_ptr<mtx> C,
	std::shared_ptr<mtx> D,
	std::shared_ptr<vec> u,
	std::shared_ptr<vec> out,
	vec x_0,
	std::shared_ptr<double> ts):
	A(A), B(B), C(C), D(D), u(u), y(out), ts(ts), x(x_0)
{
	//check if the state space is well formulated
	dimx = x_0.size();
	dimu = u->size();
	//A must by dimx X dimx
	dimArow = A->size();
	dimArow == 0 ? dimAcol = 0 : dimAcol = (*A)[0].size();
	if (dimArow != dimx or dimAcol != dimx) throw std::domain_error("State matrix dimensions do not match the state vector");
	for (auto& i : (*A)) if (i.size() != dimAcol) throw std::domain_error("A is not a matrix");
	
	//B must be dimu X dimx
	dimBrow = B->size();
	dimBrow == 0 ? dimBcol = 0 : dimBcol = (*B)[0].size();
	if (dimBrow != dimx or dimBcol != dimu) throw std::domain_error("Input matrix dimensions do not match the state or input vector");
	for (auto& i : (*B)) if (i.size() != dimBcol) throw std::domain_error("B is not a matrix");

	//C must be dimx X dimy
	dimCrow = C->size();
	dimCrow == 0 ? dimCcol = 0 : dimCcol = (*C)[0].size();
	if (dimCcol != dimx) throw std::domain_error("Observation matrix does not match the state vector");
	for (auto& i : (*C)) if (i.size() != dimCcol) throw std::domain_error("C is not a matrix");
	dimy = dimCrow;

	//D must be dimu X dimy
	dimDrow = D->size();
	dimDrow == 0 ? dimDcol = 0 : dimDcol = (*D)[0].size();
	if (dimDrow != dimy || dimDcol != dimu) throw std::domain_error("State matrix dimensions don't match the state vector");
	for (auto& i : (*A)) if (i.size() != dimAcol) throw std::domain_error("A is not a matrix");

	dx.resize(dimx);
	if (y == nullptr) y = std::make_shared<vec>();
	if (y->size() != dimy)
		y->resize(dimy);
	

}

inline cLTI::~cLTI()
{
}

inline void cLTI::Step()
{
	//euler
	for (size_t i = 0; i < dimx; i++)
	{
		dx[i] = ((*A)[i] * x).sum() + ((*B)[i] * (*u)).sum();
	}
	x += dx*(*ts);
	for (size_t i = 0; i < dimy; i++)
	{
		(*y)[i] = ((*C)[i] * x).sum() + ((*D)[i] * (*u)).sum();
	}
}




inline cGeneralSSM::cGeneralSSM(
	std::function<void(vec&, const vec&, const vec&)> f,
	std::function<void(vec&, const vec&, const vec&)> g,
	std::shared_ptr<vec> u,
	std::shared_ptr<vec> y,
	vec x_0,
	std::shared_ptr<double> ts):
	f(f), g(g) ,u(u), y(y)
{
	x = x_0;
	dx = x_0;

	k1.resize(x.size());
	k2.resize(x.size());
	k3.resize(x.size());
	k4.resize(x.size());
}

inline cGeneralSSM::~cGeneralSSM()
{
}

inline void cGeneralSSM::Step()
{
	//RK4
	f(k1, x, (*u));
	f(k2, x + k1 * (*ts) / 2., (*u));
	f(k3, x + k2 * (*ts) / 2., (*u));
	f(k4, x + k3 * (*ts), (*u));

	x += (*ts) / 6 * (k1 + 2 * k2 + 2 * k3 + k4);

	g((*y), x, (*u));
}
