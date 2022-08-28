#include "cWorldObject.h"

cMissileLTI::cMissileLTI(
	vec initial,
	std::shared_ptr<double> ts,
	sSharedDataMissile initial_data,
	sPosition initial_posit)
	:data(std::make_shared<sSharedDataMissile>(initial_data)),
	position(std::make_shared<sPosition>(initial_posit)),
	SSoutput(std::make_shared<vec>(initial)),
	ts(ts), 
	t(0),
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
	input(std::make_shared<vec>(vec{0,0})),
	target(nullptr),		//no target unless set
	aero(std::make_unique<cCoeff>(
		A,B,C,D,SSoutput,input)),
	state(std::make_unique<cLTI>(
		A, B, C, D, 
		input, 
		SSoutput, 
		vec{(*SSoutput)[0],(*SSoutput)[2]},
		ts)),
	kinematics(std::make_unique<cKinematics>(SSoutput,position,ts))
{}

cMissileLTI::~cMissileLTI()
{}

void cMissileLTI::SimulateImpl()
{
	aero->StepI();
	state->StepI();
	kinematics->Step();
}

const vec& cMissileLTI::GetState() const
{
	return *SSoutput;
}

const sPosition& cMissileLTI::GetPosition() const
{
	return *position;
}

void cMissileLTI::SetTarget(std::shared_ptr<const cWorldObject> other)
{
	target = other;
}

cTarget::cTarget(vec initial, std::shared_ptr<double> ts):
	state(std::make_shared<vec>(vec{ initial })),	//
	input(std::make_shared<vec>(vec{ initial }))
{}

cTarget::~cTarget()
{
}

void cTarget::Simulate()
{
}

const vec& cTarget::GetState() const
{
	// TODO: insert return statement here
}

const sPosition& cTarget::GetPosition() const
{
	// TODO: insert return statement here
}
