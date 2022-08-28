#pragma once
#include "utilities.h"
#include <valarray>
#include <stdexcept>
#include <memory>

class cKinematics
{
public:
	cKinematics(
		std::shared_ptr<vec> y, 
		std::shared_ptr<sPosition> pos,
		std::shared_ptr<double> ts):
		y(y), pos(pos), ts(ts), l(0), ur(0)
	{
		if (y->size() < 2) 
			throw std::domain_error("Insufficient state information");
		if (pos == nullptr)	
			pos = std::make_shared<sPosition>();
	}

	void Step()
	{
		// y = {u, r}
		// pos = {x, y}
		// Psi
		// dPsi
		if ((*y)[1] != 0.)			//arc
		{
			ur = (*y)[0] / (*y)[1];
			pos->dPsi = (*y)[1] * (*ts);
			l = 2 * ur * sin(pos->dPsi / 2.);
			pos->Psi += pos->dPsi;
			pos->pos[0] += l * sin(pos->Psi / 2.);
			pos->pos[1] += l * cos(pos->Psi / 2.);
		}
		else						//straight
		{
			l = (*y)[0] * (*ts);
			pos->dPsi = 0.;
			pos->pos[0] += l * sin(pos->Psi);
			pos->pos[1] += l * cos(pos->Psi);
		}
	}

	~cKinematics() {
	}

private:
	std::shared_ptr<vec> y;
	std::shared_ptr<sPosition> pos;
	std::shared_ptr<double> ts;
	double ur, l;

};

