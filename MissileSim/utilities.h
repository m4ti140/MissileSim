#pragma once
#include<math.h>
#include<valarray>
#include<memory>

constexpr double PI = 3.14159265359;
constexpr double _2PI = 2 * PI;

//constants

constexpr double g_0 = 9.80665;
//constexpr double V_m_min = 300;
constexpr double r_min = 15;
//constexpr double I_sp = 235;
//constexpr double m_0 = 230;
//constexpr double t_burn = 6;
//constexpr double C_D_Ma[2][7] = { {0.5, 0.9, 1.2, 1.5, 2., 3., 4.},{0.012,0.015,0.046,0.044,0.038,0.030,0.026} };
//constexpr double K = 1.5;
//constexpr double g_m_max = 30;
constexpr double g_t_max = 6;
//constexpr double S = 0.9;
//constexpr double m_dot_0 = 10;
constexpr double rho = 1.00649; //2km at isa
constexpr double a = 332.53;
//constexpr double psi_m_0 = 0;
constexpr double psi_t_0 = 180;
//constexpr double Cq = 0.5 * rho * S;
//constexpr double V_e = I_sp * g_0;
//constexpr double C_Y_max = 0.1689029567786894;

using vec = std::valarray<double>;
using mtx = std::valarray<std::valarray<double>>;

struct sSharedDataMissile
{
	double Ma, qS, qSmu, C_L_req, V_min, I_sp, m_0, t_burn, K, g_max, S, m_dot_0, Psi_0, x_0, y_0, C_Y_max, Cq, V_e;
	double C_D_Ma[2][7] = { {0.5, 0.9, 1.2, 1.5, 2., 3., 4.},{0.012,0.015,0.046,0.044,0.038,0.030,0.026} };
};

struct sPosition
{
	vec pos = { 0,0 };
	cAzimuth<double> Psi = 0;
	cAzimuth<double> dPsi = 0;
};

namespace interp {

	template<typename T, size_t S>
	inline T Linear(const T(&y)[S], const T(&x)[S], const T& xt, size_t& lasti)
	{
		size_t i = lasti, j = lasti;
		bool jcond = x[j] > xt;
		while (x[i] < xt && jcond)
		{
			i++; j--;
			if (i == S)
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

namespace integrators
{
	//template <typename T, typename U>
	//inline void RK4(T(*) f(T*x,T*u),T* x, U* t, U* ts)
	//{
	//	T k = f(x, u, t);
	//}
};

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
	cAzimuth() {};
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

using azimuth_vec = std::valarray < cAzimuth<double>>;