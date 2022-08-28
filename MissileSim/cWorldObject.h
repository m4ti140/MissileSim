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

//template <class Child>
//class cMissile : public cWorldObject
//{
//public:
//	cMissile(vec initial,
//		sSharedDataMissile initialData,
//		sPosition initialPosit)
//		:SSoutput(std::make_shared<vec>(initial)),
//		position(std::make_shared<sPosition>(initialPosit)),
//		data(std::make_shared<sSharedDataMissile>(initialData)),
//		target(nullptr)		//no target unless set
//	{
//		this->InitI();
//	}
//	void InitI()
//	{
//		static_cast<Child*>(this)->Init(SSoutput, data, position);
//	}
//	void Simulate()
//	{
//		static_cast<Child*>(this)->SimulateImpl(target);
//	}
//	std::shared_ptr<const cWorldObject> GetTarget()
//	{
//		return target;
//	}
//	const vec& GetState() const;
//	const sPosition& GetPosition() const;
//	void SetTarget(std::shared_ptr<const cWorldObject> other);
//
//private:
//	std::shared_ptr<vec> SSoutput;
//	std::shared_ptr<sPosition> position;
//	std::shared_ptr<sSharedDataMissile> data;
//	std::shared_ptr<const cWorldObject> target;
//};

class cMissileLTI : public cWorldObject
{
public:
	cMissileLTI(vec initial,
		std::shared_ptr<double> ts,
		sSharedDataMissile initialData,
		sPosition initialPosit);
	~cMissileLTI();
	//void Init(std::shared_ptr<vec> output,
	//	std::shared_ptr<sSharedDataMissile> data,
	//	std::shared_ptr<sPosition> posit);
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
	std::shared_ptr<const cWorldObject> target;
	std::unique_ptr<cCoeff> aero;
	std::unique_ptr<cLTI> state;
	std::unique_ptr<cKinematics> kinematics;
	std::unique_ptr<cPN> controller;
	std::shared_ptr<double> ts;
	double t;
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
