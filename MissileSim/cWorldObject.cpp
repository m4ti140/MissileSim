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
					{0, 1}})),
	D(std::make_shared<mtx>(mtx{
						{0, 0},
						{0, 0},
						{0, 0} })),
	SSoutput(std::make_shared<vec>(initial)),
	position(std::make_shared<sPosition>(initialPosit)),
	data(std::make_shared<sSharedDataMissile>(initialData)),
	target(nullptr),	//no target unless set
	input(std::make_shared<vec>(vec{ 0,0 })),
	aero(
		A, B, C, D, SSoutput, input, data),
	state(
	A, B, C, D,
	input,
	SSoutput,
	vec{ (*SSoutput)[0],(*SSoutput)[2] },
	ts),
	kinematics(SSoutput, position, ts),
	controller(SSoutput, position, input, data, ts, 5)
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
	aero.Calculate(t);
	if (target != nullptr)
	{
		controller.Setpoint(target->GetState(), target->GetPosition());
	}
	controller.StepI();
	state.Step();
	kinematics.Step();
	t += *ts;
}

inline const vec& cMissileLTI::GetState() const
{
	return *SSoutput;
}

inline const sPosition& cMissileLTI::GetPosition() const
{
	return *position;
}

inline void cMissileLTI::SetTarget(std::shared_ptr<const cWorldObject> other)
{
	target = other;
}


cMissileNonLinear::cMissileNonLinear(
	vec initial, 
	std::shared_ptr<double> ts, 
	sSharedDataMissile initialData, 
	sPosition initialPosit):
	input(std::make_shared<vec>(vec{0.,data->m_dot_0})),
	SSoutput(std::make_shared<vec>(initial)),
	position(std::make_shared<sPosition>(initialPosit)),
	data(std::make_shared<sSharedDataMissile>(initialData)),
	target(nullptr),
	ts(ts),t(std::make_shared<double>(0.)),
	state(cMissileNonLinear::f, cMissileNonLinear::g, input, SSoutput, vec{(*SSoutput)[0],(*SSoutput)[2]}, data, ts, t),
	kinematics(SSoutput,position,ts),
	controller(SSoutput,position,input,data,ts,5)
{
	data->lasti = 1;
	data->Cq = 0.5 * rho * data->S;
	data->V_e = data->I_sp * g_0;
}

cMissileNonLinear::~cMissileNonLinear()
{
}

void cMissileNonLinear::Simulate()
{
	if (*t >= data->t_burn) (*input)[1] = 0;
	if (target != nullptr)
	{
		controller.Setpoint(target->GetState(), target->GetPosition());
	}
	controller.StepI();
	//(*input)[0] = data->C_Y_max;
	state.Step();
	kinematics.Step();
	*t += *ts;
}

const vec& cMissileNonLinear::GetState() const
{
	return *SSoutput;
}

const sPosition& cMissileNonLinear::GetPosition() const
{
	return *position;
}

void cMissileNonLinear::SetTarget(std::shared_ptr<const cWorldObject> other)
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