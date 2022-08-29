#pragma once
#include "utilities.h"
#include "cStateSpace.h"
#include "cKinematics.h"
#include "cAero.h"
#include "cControl.h"

class cWorldObject
{
public:
	virtual ~cWorldObject() {};
	virtual void Simulate() = 0;
	virtual const vec& GetState() const = 0;
	virtual const sPosition& GetPosition() const = 0;
	virtual void SetTarget(std::shared_ptr<const cWorldObject> other) = 0;
	
};

class cMissileLTI : public cWorldObject
{
public:
	cMissileLTI(vec initial,
		std::shared_ptr<double> ts,
		sSharedDataMissile initialData,
		sPosition initialPosit);
	~cMissileLTI();

	void Simulate();
	const vec& GetState() const;
	const sPosition& GetPosition() const;
	void SetTarget(std::shared_ptr<const cWorldObject> other);

private:
	std::shared_ptr<mtx> A, B, C, D;
	std::shared_ptr<vec> input;
	std::shared_ptr<vec> SSoutput;
	std::shared_ptr<sPosition> position;
	std::shared_ptr<sSharedDataMissile> data;
	std::shared_ptr<double> ts;
	std::shared_ptr<const cWorldObject> target;
	cCoeff aero;
	cLTI state;
	cKinematics kinematics;
	cPN controller;
	double t;
};

class cMissileNonLinear : public cWorldObject
{
public:
	cMissileNonLinear(vec initial,
		std::shared_ptr<double> ts,
		sSharedDataMissile initialData,
		sPosition initialPosit);
	~cMissileNonLinear();
	void Simulate();
	const vec& GetState() const;
	const sPosition& GetPosition() const;
	void SetTarget(std::shared_ptr<const cWorldObject> other);

private:
	static inline void f(vec& k, const vec& x, const vec& u, const std::shared_ptr<sSharedDataMissile>& data, const double& t)
	{
		if (t < data->t_burn) k[1] = -data->m_dot_0;
		else k[1] = 0;

		double m = std::max(data->m_0 + k[1] * t,170.);

		data->qS = data->Cq * x[0] * x[0];
		data->Ma = x[0] / a;
		data->C_L_req = m * g_0/data->qS;
		double C_D = interp::Linear(data->C_D_Ma[1], data->C_D_Ma[0], data->Ma, data->lasti) + data->K * data->C_L_req * data->C_L_req + data->K * u[0] * u[0];
		k[0] = (data->V_e * -k[1] - data->qS * C_D) / m;

 		//if (std::isnan(k[0]) || std::isnan(k[1])) throw std::domain_error("output is NAN");
	}
	static inline void g(vec& y, const vec& x, const vec& u, const std::shared_ptr<sSharedDataMissile>& data, const double& t)
	{
		y[0] = x[0];
		y[1] = data->Cq * x[0] * u[0]/x[1];
		y[2] = x[1];
		//if (std::isnan(y[0]) || std::isnan(y[1]) || std::isnan(y[2])) throw std::domain_error("output is NAN");
	}
	std::shared_ptr<sSharedDataMissile> data;
	std::shared_ptr<vec> input, SSoutput;
	std::shared_ptr<sPosition> position;
	std::shared_ptr<const cWorldObject> target;
	std::shared_ptr<double> ts, t;
	cGeneralSSM<sSharedDataMissile> state;
	cKinematics kinematics;
	cPN controller;
};

class cTarget : public cWorldObject
{
public:
	cTarget(vec initialState, sPosition initialPosition, sSharedDataTarget initialData, std::shared_ptr<double> ts);
	~cTarget();
	void Simulate();
	const vec& GetState() const;
	const sPosition& GetPosition() const;
	void SetTarget(std::shared_ptr<const cWorldObject> other){} //useless for targets
private:
	std::shared_ptr<vec> state;
	std::shared_ptr<sPosition> position;
	std::shared_ptr<double> ts, t;
	std::shared_ptr<sSharedDataTarget> data;
	cKinematics kinematics;
	cRandomPath controller;
};

namespace range
{
	inline double GetRange(const sPosition& pos1, const sPosition& pos2)
	{
		return sqrt(pow(pos1.pos[0] - pos2.pos[0], 2) + pow(pos1.pos[1] - pos2.pos[1], 2));
	}
};
