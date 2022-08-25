#pragma once
#include <valarray>
#include <vector>
#include <stdexcept>
#include "utilities.h"

template<typename T>
class cStateSpace
{
public:
	virtual ~cStateSpace() = 0;
	virtual std::valarray<T>* Step() = 0;
};

template<typename T>
class cLTI : public cStateSpace<T>
{
public:
	cLTI(
		std::valarray<std::valarray<T>>* A,			//using pointers to avoid function calls and copy operations
		std::valarray<std::valarray<T>>* B,
		std::valarray<std::valarray<T>>* C,
		std::valarray<std::valarray<T>>* D,
		std::valarray<T>* u,
		std::valarray<T>* y,
		std::valarray<T> x_0,
		T* ts);
	~cLTI();

	std::valarray<T>* Step() override;
private:
	std::valarray<std::valarray<T>>* A, * B, * C, * D;
	std::valarray<T>* u, * y, x, dx;
	T* ts;
	size_t dimArow, dimAcol, dimBrow, dimBcol, dimCrow, dimCcol, dimDrow, dimDcol, dimx, dimy, dimu;
};

template<typename T>
class cGeneralSSM : public cStateSpace<T>
{
public:
	cGeneralSSM(
		void (*f)(std::valarray<T>&, const std::valarray<T>&, const std::valarray<T>&),
		void (*g)(std::valarray<T>*, const std::valarray<T>&, const std::valarray<T>&),
		std::valarray<T>* u,
		std::valarray<T>* y,
		std::valarray<T> x_0, 
		T* ts);
	~cGeneralSSM();

	std::valarray<T>* Step() override;
private:
	void (*f)(std::valarray<T>&, const std::valarray<T>&, const std::valarray<T>&);
	void (*g)(std::valarray<T>*, const std::valarray<T>&, const std::valarray<T>&);
	std::valarray<T>* u, * y, x, dx, k1, k2, k3, k4;
	T* ts;
};

template<typename T>
inline cLTI<T>::cLTI(
	std::valarray<std::valarray<T>>* A, 
	std::valarray<std::valarray<T>>* B, 
	std::valarray<std::valarray<T>>* C, 
	std::valarray<std::valarray<T>>* D, 
	std::valarray<T>* u,
	std::valarray<T>* y, 
	std::valarray<T> x_0, 
	T* ts) :
	A(A), B(B), C(C), D(D), u(u), y(y), ts(ts)
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
	dimy = dimCcol;

	//D must be dimu X dimy
	dimDrow = D->size();
	dimDrow == 0 ? dimDcol = 0 : dimDcol = (*D)[0].size();
	if (dimDrow != dimy or dimDcol != dimu) throw std::domain_error("State matrix dimensions don't match the state vector");
	for (auto& i : (*A)) if (i.size() != dimAcol) throw std::domain_error("A is not a matrix");

	//x = new std::valarray<T>(x_0);
	//dx = new std::valarray<T>(x_0);
	//temp = nullptr;

	x = x_0;
	dx.resize(dimx);
	if (this->y == nullptr);
	{
		this->y = new std::valarray<T>;
		this->y->resize(dimy);
	}

}

template<typename T>
inline cLTI<T>::~cLTI()
{
}

template<typename T>
inline std::valarray<T>* cLTI<T>::Step()
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
	return y;
}



template<typename T>
inline cGeneralSSM<T>::cGeneralSSM(
	void(*f)(std::valarray<T>&, const std::valarray<T>&, const std::valarray<T>&),
	void(*g)(std::valarray<T>*, const std::valarray<T>&, const std::valarray<T>&),
	std::valarray<T>* u,
	std::valarray<T>* y,
	std::valarray<T> x_0,
	T* ts):
	f(f), g(g) ,u(u), y(y)
{
	x = x_0;
	dx = x_0;
	if (this->y == nullptr)
		this->y = new std::valarray<T>;

	k1.resize(x.size());
	k2.resize(x.size());
	k3.resize(x.size());
	k4.resize(x.size());
}

template<typename T>
inline cGeneralSSM<T>::~cGeneralSSM()
{
}

template<typename T>
inline std::valarray<T>* cGeneralSSM<T>::Step()
{
	//RK4
	f(k1, x, (*u));
	f(k2, x + k1 * (*ts) / 2., (*u));
	f(k3, x + k2 * (*ts) / 2., (*u));
	f(k4, x + k3 * (*ts), (*u));

	x += (*ts) / 6 * (k1 + 2 * k2 + 2 * k3 + k4);

	g(y, x, (*u));
	return y;	
}

template<typename T>
inline cStateSpace<T>::~cStateSpace()
{
}
