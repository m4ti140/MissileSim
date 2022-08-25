#pragma once
#include<math.h>

constexpr double PI = 3.14159265359;
constexpr double _2PI = 2 * PI;

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
