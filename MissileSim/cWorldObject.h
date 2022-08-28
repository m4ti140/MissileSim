#pragma once
#include "utilities.h"
#include "cStateSpace.h"
#include "cKinematics.h"
#include "cAero.h"

class cWorldObject
{
public:
	virtual ~cWorldObject() {};
	virtual void Simulate() = 0;
	virtual const vec& GetState() const = 0;
	virtual const sPosition& GetPosition() const = 0;
	
};

template <class Child>
class cMissile : public cWorldObject
{
	void Simulate()
	{
		static_cast<Child*>(this)->SimulateImpl();
	}
};

class cMissileLTI : public cMissile<cMissileLTI>
{
public:
	cMissileLTI(
		vec initial,
		std::shared_ptr<double> ts,
		sSharedDataMissile initial_data,
		sPosition initial_posit);
	~cMissileLTI();
	void SimulateImpl();
	const vec& GetState() const;
	const sPosition& GetPosition() const;
	void SetTarget(std::shared_ptr<const cWorldObject>);

private:
	std::shared_ptr<mtx> A, B, C, D;
	std::shared_ptr<vec> SSoutput,input;
	std::shared_ptr<sPosition> position;
	std::shared_ptr<const cWorldObject> target;
	std::shared_ptr<double> ts;
	std::unique_ptr<cCoeff> aero;
	std::unique_ptr<cLTI> state;
	std::unique_ptr<cKinematics> kinematics;
	std::shared_ptr<sSharedDataMissile> data;
	double t;
};

class cTarget : public cWorldObject
{
public:
	cTarget(vec initial, std::shared_ptr<double> ts);
	~cTarget();
	void Simulate();
	const vec& GetState() const;
	const sPosition& GetPosition() const;
private:
	std::shared_ptr<vec> input, position, state;
	std::shared_ptr<double> ts;
	std::shared_ptr<azimuth_vec> direction;
	std::unique_ptr<cKinematics> kinematics;
};