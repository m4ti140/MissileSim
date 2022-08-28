#pragma once
#include "utilities.h"
#include <deque>
#include <random>
#include <algorithm>

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
		V_c = (*state)[0] * cos(lambda.Delta(position->Psi)) - targetState[0] * cos(lambda.Delta(targetPosit.Psi));
		error = prevLambda.Delta(lambda) / (*ts);
		(*input)[0] = std::clamp((*state)[2] / data->qS * gain * error * V_c,-data->C_Y_max,data->C_Y_max);
		prevLambda = lambda;
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

class cRandomPath : public cControl<cRandomPath>
{
public:
	cRandomPath(std::shared_ptr<vec> state,
		std::shared_ptr<sSharedDataTarget> data,
		std::shared_ptr<double> ts,
		std::shared_ptr<double> t):
		state(state),data(data),ts(ts),t(t)
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<double> g_turn(-(data->g_t_max), data->g_t_max);
		std::uniform_real_distribution<double> t_turn(0, 10);
		std::uniform_int_distribution<int> turns(3, 8);

		for (int i = 0; i < turns(gen); i++)
		{
			schedule[0].push_back(t_turn(gen));
			schedule[0].push_back(t_turn(gen));
			schedule[1].push_back(g_turn(gen));
			schedule[1].push_back(0);
		}
		std::sort(schedule[0].begin(), schedule[0].end());
	};
	void Step()
	{
		if (!schedule[0].empty() && *t > schedule[0].front())
		{
			if (schedule[1].front() == 0.) (*state)[1] = 0;
			else (*state)[1] = schedule[1].front() * g_0 / (*state)[0];

			schedule[0].pop_front();
			schedule[1].pop_front();
		}
	}
private:
	std::deque<double> schedule[2];
	std::shared_ptr<vec> state;
	std::shared_ptr<sSharedDataTarget> data;
	std::shared_ptr<double> ts, t;
};