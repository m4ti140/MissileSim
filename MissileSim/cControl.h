#pragma once
#include "utilities.h"

template<class Child>
class cControl
{
public:
	void StepI() 
	{
		static_cast<Child*>(this)->Step();
	}
};

class cPN : public cControl<cPN>
{
public:
	cPN(std::shared_ptr<vec> state,
		std::shared_ptr<sPosition> position,
		std::shared_ptr<vec> input,
		std::shared_ptr<sSharedDataMissile> data,
		std::shared_ptr<double> ts,
		double gain):
		state(state), position(position),input(input), data(data), lambda(0), prevLambda(0), error(0), V_c(0), ts(ts),gain(gain)
	{}
	~cPN(){}
	void Setpoint(const vec& targetState, const sPosition& targetPosit)
	{
		this->targetPosit = targetPosit;
		this->targetState = targetState;
	}
	void Step()
	{
		//input = {CY,m_dot}
		//state = {u,r,m,ay}
		lambda = atan2(targetPosit.pos[0]- position->pos[0], targetPosit.pos[1] - position->pos[1]);
		V_c = (*state)[0] * cos(position->Psi - lambda) - targetState[0] * cos(targetPosit.Psi - lambda);
		error = (lambda - prevLambda) / (*ts);
		(*input)[0] = (*state)[2] / data->qS * gain * error * V_c;
	}

private:
	std::shared_ptr<sSharedDataMissile> data;
	std::shared_ptr<sPosition> position;
	std::shared_ptr<vec> state, input;
	vec targetState;
	sPosition targetPosit;
	cAzimuth<double> lambda, prevLambda;
	double V_c,error,gain;
	std::shared_ptr<double> ts;
};

class cRandom : public cControl<cRandom>
{

};