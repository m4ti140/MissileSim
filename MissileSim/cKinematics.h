#pragma once
#include "utilities.h"
#include <valarray>
#include <stdexcept>

template <typename T>
class cKinematics
{
public:
	cKinematics(std::valarray<T>* y, std::valarray<T>* pos, std::valarray<cAzimuth<T>>* dir, T* ts)
		y(y), pos(pos), dir(dir), ts(ts)
	{
		if (y->size() < 2) 
			throw std::domain_error("Insufficient state information");
		if (pos == nullptr)	
			pos = new std::valarray<T>{ 0,0 };
		if (pos->size() < 2)
			pos->resize(2);
		if (dir == nullptr)
			dir = new std::valarray<cAzimuth<T>>{ 0,0 };
		if (dir->size() < 2)
			dir->resize(2);
			
	}

	std::valarray<T>* Step()
	{
		// y = {u, r}
		// pos = {x, y}
		// dir = {Psi, deltaPsi}
		if ((*y)[1] != 0.)
		{
			ur = (*y)[0] / (*y)[1];
			(*dir)[1] = (*y)[1] * (*ts);
			l = 2 * ur * sin((*dir)[1] / 2.);
			(*dir)[0] += (*dir)[1];
			(*pos)[0] += l * sin((*dir)[0] / 2.);
			(*pos)[1] += l * cos((*dir)[0] / 2.);
		}
		else
		{
			l = (*y)[0] * ts;
			(*dir)[1] = 0.;
			(*pos)[0] += l * sin((*dir)[0]);
			(*pos)[1] += l * cos((*dir)[0]);
		}
	}

	~cKinematics() {
		//we own nothing
	}

private:
	std::valarray<T>* y, * pos;
	std::valarray<T>* dir;
	T* ts, ur, l;

};

