#pragma once
#include <memory>
#include <valarray>
#include <stdexcept>
#include <utility>
#include "utilities.h"


//calculates aerodynamic coefficients

class cCoeff
{
public:
	cCoeff(
		std::shared_ptr<std::valarray<std::valarray<double>>> A,
		std::shared_ptr<std::valarray<std::valarray<double>>> B,
		std::shared_ptr<std::valarray<std::valarray<double>>> C,
		std::shared_ptr<std::valarray<std::valarray<double>>> D,
		std::shared_ptr<std::valarray<double>> y,
		std::shared_ptr<std::valarray<double>> u,
		std::shared_ptr<sSharedDataMissile> dat
		) :A(A), B(B), C(C), D(D), y(y), u(u), data(dat)
		
	{
		if (y->size() != C->size()) throw std::domain_error("invalid output vector");
		if (u->size() != 0 && u->size() != (*C)[0].size()) throw std::domain_error("invalid input vector");

		data->lasti = 1;
		data->Cq = 0.5 * rho * data->S;
		data->V_e = data->I_sp * g_0;
		if (data->C_D_Ma[0].size() != data->C_D_Ma[1].size()) throw std::domain_error("drag table invalid");
		(*u)[1] = data->m_dot_0;
	}
	~cCoeff() {}
	void Calculate(const double& t)
	{
		//y = {u, r, m, ay}
		//x = {u, m}
		//A = {	{X_u, 0},
		//		{0,	0}}
		//B = {	{X_CY, X_mdot},
		//		{0		-1	}}
		//C = {	{1,0},
		//		{0,0},
		//		{0,1},
		//D = {	{0,0},
		//		{N_CY,0},
		//		{0,0},

		(*B)[0][1] = data->V_e / (*y)[2]; //X_mdot
		data->Ma = (*y)[0] / a;
		data->qS = data->Cq * (*y)[0] * (*y)[0];
		data->qSmu = data->qS / (*y)[2] / (*y)[0];
		data->C_L_req = g_0 * (*y)[2] / data->qS;

		(*A)[0][0] = -data->qSmu * (interp::Linear(data->C_D_Ma[1], data->C_D_Ma[0], data->Ma, data->lasti) + data->K * (data->C_L_req * data->C_L_req)); //X_u

		(*D)[1][0] = data->qSmu; //N_CY
		(*B)[0][0] = -data->qS/(*y)[2] * (*u)[0] * data->K; //X_CY
		//(*D)[3][0] = (*D)[1][0] * (*y)[0];
		if (t >= 6)
		{
			(*u)[1] = 0; //m_dot;
		}
	}
private:
	std::shared_ptr<mtx> A, B, C, D;
	std::shared_ptr<vec> y, u;
	//double t;
	std::shared_ptr<sSharedDataMissile> data;
};

