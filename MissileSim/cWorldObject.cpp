#include "cWorldObject.h"

cMissileLTI::cMissileLTI(vec initial, std::shared_ptr<double> ts, sSharedDataMissile initialData, sPosition initialPosit)
	: ts(ts), t(0),
	A(std::make_shared<mtx>(mtx{
			{ 0, 0 },
			{ 0, 0 } })),
	B(std::make_shared<mtx>(mtx{
			{0, 0},
			{0,-1} })),
	C(std::make_shared<mtx>(mtx{
					{1, 0},
					{0, 0},
					{0, 1},
					{0, 0} })),
	D(std::make_shared<mtx>(mtx{
						{0, 0},
						{0, 0},
						{0, 0},
						{0, 0} })),
	SSoutput(std::make_shared<vec>(initial)),
	position(std::make_shared<sPosition>(initialPosit)),
	data(std::make_shared<sSharedDataMissile>(initialData)),
	target(nullptr),	//no target unless set
	input(std::make_shared<vec>(vec{ 0,0 })),
	aero(std::make_unique<cCoeff>(
		A, B, C, D, SSoutput, input, data)),
	state(std::make_unique<cLTI>(
	A, B, C, D,
	input,
	SSoutput,
	vec{ (*SSoutput)[0],(*SSoutput)[2] },
	ts)),
	kinematics(std::make_unique<cKinematics>(SSoutput, position, ts)),
	controller(std::make_unique<cPN>(SSoutput, position, input, data, ts, 5))
{}

cMissileLTI::~cMissileLTI()
{}

//void cMissileLTI::Init(std::shared_ptr<vec> output,
//	std::shared_ptr<sSharedDataMissile> data,
//	std::shared_ptr<sPosition> posit)
//{
//	//A = std::make_shared<mtx>(mtx{
//	//		{ 0, 0 },
//	//		{ 0, 0 } });
//	//B = std::make_shared<mtx>(mtx{
//	//			{0, 0},
//	//			{0,-1} });
//	//C = std::make_shared<mtx>(mtx{
//	//				{1, 0},
//	//				{0, 0},
//	//				{0, 1},
//	//				{0, 0} });
//	//D = std::make_shared<mtx>(mtx{
//	//					{0, 0},
//	//					{0, 0},
//	//					{0, 0},
//	//					{0, 0} });
//	aero = std::make_unique<cCoeff>(
//		A, B, C, D, output, input, data);
//	state = std::make_unique<cLTI>(
//		A, B, C, D,
//		input,
//		output,
//		vec{ (*output)[0],(*output)[2] },
//		ts);
//	kinematics = std::make_unique<cKinematics>(output, posit, ts);
//	controller = std::make_unique<cPN>(output, posit, input, data, ts, 5);
//}

void cMissileLTI::Simulate()
{
	aero->StepI(t);
	state->StepI();
	kinematics->Step();
	if (target != nullptr)
	{
		controller->Setpoint(target->GetState(), target->GetPosition());
	}
	controller->StepI();
	t += *ts;
}

//template<class Child>
inline const vec& cMissileLTI::GetState() const
{
	return *SSoutput;
}

//template<class Child>
inline const sPosition& cMissileLTI::GetPosition() const
{
	return *position;
}

//template<class Child>
inline void cMissileLTI::SetTarget(std::shared_ptr<const cWorldObject> other)
{
	target = other;
}

cTarget::cTarget(vec initialState, sPosition initialPosition, sSharedDataTarget initialData, std::shared_ptr<double> ts) :
	state(std::make_shared<vec>(initialState)),
	position(std::make_shared<sPosition>(initialPosition)),
	data(std::make_shared<sSharedDataTarget>(initialData)),
	ts(ts), t(std::make_shared<double>(0)),
	kinematics(state,position,ts),
	controller(state,data,ts,t)
{}

cTarget::~cTarget()
{}

void cTarget::Simulate()
{
	kinematics.Step();
	controller.StepI();
	*t += *ts;
}

const vec& cTarget::GetState() const
{
	return *state;
}

const sPosition& cTarget::GetPosition() const
{
	return *position;
}
