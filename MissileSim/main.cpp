#include "cWorldObject.h"
#include <chrono>

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::time_point;
using std::chrono::milliseconds;

int main()
{
	sSharedDataMissile missileData;
	missileData.V_min = 300;
	missileData.I_sp = 235;
	missileData.m_0 = 230;
	missileData.t_burn = 6;
	missileData.C_D_Ma[0] = { 0.5, 0.9, 1.2, 1.5, 2., 3., 4. };
	missileData.C_D_Ma[1] = {0.012,0.015,0.046,0.044,0.038,0.030,0.026};
	missileData.K = 1.5;
	missileData.g_max = 30;
	missileData.S = 0.9;
	missileData.m_dot_0 = 10;
	missileData.C_Y_max = 0.1689029567786894; //obtained by running the sim with control disabled to obtain max speed, then calculating C_Y=g_max*g_0 * m / qS

	sPosition missilePosit;
	missilePosit.dPsi = 0;
	missilePosit.Psi = 0;
	missilePosit.pos = { 0,0 };

	vec missileState = { missileData.V_min,0,missileData.m_0,0 };

	sSharedDataTarget targetData;
	targetData.g_t_max = 6;
	targetData.V = 300;

	sPosition targetPosit;
	targetPosit.dPsi = 0;
	targetPosit.pos = { 0,10000 };
	targetPosit.Psi = PI;

	vec targetState = { targetData.V,0 };

	std::shared_ptr<double> ts = std::make_shared<double>(1. / 120.);

	std::vector<std::shared_ptr<cWorldObject>> objects;
	objects.push_back(std::make_shared<cMissileLTI>(missileState, ts, missileData, missilePosit));
	objects.push_back(std::make_shared<cTarget>(targetState, targetPosit, targetData, ts));

	objects[0]->SetTarget(objects[1]);

	double r = range::GetRange(missilePosit, targetPosit);
	double V_m = objects[0]->GetState()[0];

	FILE* f;
	errno_t err = fopen_s(&f,"out.csv", "w");
	if (f == NULL)
	{
		printf("failed to open file");
		return 0;
	}

	double worldTime = 0;
	double v_max = 0;
#ifndef _DEBUG
	fprintf(f,"time,xm,ym,xt,yt\n");
	std::vector < std::array< double,5> > results;
	results.reserve(4000);
#else
	fprintf(f,"time,xm,ym,xt,yt,Psim,dPsim,m_u,m_r,Psit,dPsit,t_u,t_r\n");
	std::vector < std::array< double, 13> > results;
	results.reserve(4000);
#endif

	time_point startTime = high_resolution_clock::now();

	while (r > r_min && V_m >= missileData.V_min)
	{
		for (auto& object : objects)
			object->Simulate();
		missilePosit = objects[0]->GetPosition();
		targetPosit = objects[1]->GetPosition();
#ifndef _DEBUG
		results.push_back({ worldTime, missilePosit.pos[0], missilePosit.pos[1], targetPosit.pos[0], targetPosit.pos[1] });
#else
		missileState = objects[0]->GetState();
		targetState = objects[1]->GetState();
		results.push_back({ worldTime, missilePosit.pos[0], missilePosit.pos[1], targetPosit.pos[0], targetPosit.pos[1], missilePosit.Psi, missilePosit.dPsi,
			missileState[0],missileState[1],targetPosit.Psi,targetPosit.dPsi,targetState[0],targetState[1]});
#endif
		worldTime += *ts;
		r = range::GetRange(missilePosit, targetPosit);
		V_m = objects[0]->GetState()[0];
		v_max = std::max(v_max, V_m);
	}

	time_point endTime = high_resolution_clock::now();

	printf("Maximum missile speed: %f m/s, Sim time: %d ms, writing file",v_max, duration_cast<milliseconds>(endTime - startTime));

	for (auto& i : results)
	{
#ifndef _DEBUG
		fprintf(f, "%f,%f,%f,%f,%f\n", i[0], i[1], i[2], i[3], i[4]);
#else
		fprintf(f, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10], i[11], i[12]);
#endif
	}


	fclose(f);
	
	return 0;

}