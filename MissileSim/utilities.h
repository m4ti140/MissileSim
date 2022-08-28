#pragma once
#include<math.h>
#include<valarray>
#include<memory>
#include<array>
#include<vector>

//constants

constexpr double PI = 3.14159265359;
constexpr double _2PI = 2 * PI;
constexpr double g_0 = 9.80665;
constexpr double r_min = 15;
constexpr double rho = 1.00649; //2km at isa
constexpr double a = 332.53;

using vec = std::valarray<double>;
using mtx = std::valarray<std::valarray<double>>;


namespace interp {

	template<typename T>
	inline T Linear(const std::vector<T>& y, const std::vector<T>& x, const T& xt, size_t& lasti)
	{
		size_t i = lasti, j = lasti;
		bool jcond = x[j] > xt;
		while (x[i] < xt && jcond)
		{
			i++; j--;
			if (i == y.size())
			{
				while (x[j] > xt) j--;
				i = ++j;
				break;
			}
			else if (j == 0)
			{
				while (x[i] > xt) i++;
				break;
			}
			jcond = x[j] > xt;
		}
		if (!jcond) i = ++j;
		lasti = i;
		return y[i - 1] + (xt - x[i - 1]) * (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
	}

}

template <typename T>
class cAzimuth;
namespace az 
{
	template <typename T>
	T Delta(const cAzimuth<T>& lhs, const cAzimuth<T>& rhs)
	{
		return fmod(rhs.val - lhs.val, _2PI);
	}
};
template <typename T>
class cAzimuth
{
public:
	cAzimuth():val(0){};
	cAzimuth(const T& val):val(val) {
		_Norm();
	}
	operator T() const { return val; }
	cAzimuth<T>& operator=(const T& other)
	{
		this->val = other;
		_Norm();
		return *this;
	}
	cAzimuth<T>& operator+=(const cAzimuth<T>& rhs)
	{
		this->val += rhs.val;
		_Norm();
		return *this;
	}
	cAzimuth<T>& operator+=(const T& rhs)
	{
		this->val += rhs;
		_Norm();
		return *this;
	}
	friend cAzimuth<T> operator+(cAzimuth<T> lhs, const cAzimuth<T>& rhs)
	{
		lhs += rhs;
		return lhs;
	}
	friend cAzimuth<T> operator+(cAzimuth<T> lhs, const T& rhs)
	{
		lhs += rhs;
		return lhs;
	}
	cAzimuth<T>& operator-=(const cAzimuth<T>& rhs)
	{
		this->val -= rhs.val;
		_Norm();

		return *this;
	}
	cAzimuth<T>& operator-=(const T& rhs)
	{
		this->val -= rhs;
		_Norm();
		return *this;
	}
	friend cAzimuth<T> operator-(cAzimuth<T> lhs, const cAzimuth<T>& rhs)
	{
		lhs -= rhs;
		return lhs;
	}
	friend cAzimuth<T> operator-(cAzimuth<T> lhs, const T& rhs)
	{
		lhs -= rhs;
		return lhs;
	}
	friend T az::Delta(const cAzimuth<T>& lhs, const cAzimuth<T>& rhs);
	T Delta(const cAzimuth<T>& rhs)
	{
		return az::Delta(*this, rhs);
	}

	cAzimuth<T>& operator*=(const T& rhs)
	{
		this->val *= rhs;
		_Norm();
		return *this;
	}
	friend cAzimuth<T> operator*(cAzimuth<T> lhs, const T& rhs)
	{
		lhs *= rhs;
		return lhs;
	}	
	cAzimuth<T>& operator/=(const T& rhs)
	{
		this->val /= rhs;
		_Norm();
		return *this;
	}
	friend cAzimuth<T> operator/(cAzimuth<T> lhs, const T& rhs)
	{
		lhs /= rhs;
		return lhs;
	}
	T Degrees()
	{
		return val / PI * 180;
	}

private:
	void _Norm() {
		while (val < 0) val += _2PI;
		while (val >= _2PI) val -= _2PI;
	}

	T val;
};

struct sSharedDataMissile
{
	double Ma, qS, qSmu, C_L_req, V_min, I_sp, m_0, t_burn, K, g_max, S, m_dot_0, Psi_0, x_0, y_0, C_Y_max, Cq, V_e;
	std::vector<double> C_D_Ma[2];
};

struct sSharedDataTarget
{
	double V, g_t_max;
};

struct sPosition
{
	vec pos;
	cAzimuth<double> Psi;
	double dPsi;
};

using azimuth_vec = std::valarray < cAzimuth<double>>;