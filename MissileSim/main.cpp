#include "cWorldObject.h"
#include <chrono>

using namespace std::chrono;

int main(int argc, char* argv[])
{
	

	//timestep

	std::shared_ptr<double> ts = std::make_shared<double>(1. / 120.);

	//output file
	std::string path="out.txt";

	//use nonlinear model?
	bool nonlinear=false;

	//process arguments

	double r_0, V_t;

	//help string
	const char* helpstring = R"(
MissileSim target_velocity target_range [options]
or
MissileSim -help
options:
	-help - display this help
	-nonlinear - use nonlinear model for the missile
		Default model is linear time invariant state space representation with euler solver and coefficients recalculated every run of main loop
		with -nonlinear full equations of motions and RK4 solver are used
	-ts - set simulation timestep. Must be > 0
	-path - full path to output file. If not provided, out.txt is created in executable's directory
)";


	if (argc == 2 && !strcmp("-help",argv[1]))
	{
		printf(helpstring);
		return 0;
	}
	else if (argc < 3)
	{
		printf("Not enough arguments. See MissileSim -help for arguments\n");
		return 0;
	}
	else
	{
		V_t = strtod(argv[1],nullptr);
		if (errno == ERANGE || V_t<0)
		{
			printf("Target velocity out of acceptable range");
			return 0;
		}
		r_0 = strtod(argv[2], nullptr);
		if (errno == ERANGE || r_0<0)
		{
			printf("Target range out of acceptable range");
			return 0;
		}
		if (argc > 3)
		{
			for (int i = 3; i < argc; i++)
			{
				if (!strcmp(argv[i],"-nonlinear"))
				{
					nonlinear = true;
				}
				else if (!strcmp(argv[i],"-path"))
				{
					if (argc == i + 1)
					{
						printf("-path used but path not specified");
						return 0;
					}
					else
						path = argv[++i];
				}
				else if (!strcmp(argv[i],"-ts"))
				{
					if (argc == i + 1)
					{
						printf("-ts used but timestep not specified");
						return 0;
					}
					else
					{
						*ts = strtod(argv[++i], nullptr);
						if (*ts <= 0)
						{
							printf("Timestep cannot be 0 or negative");
							return 0;
						}
					}
				}
				else if (!strcmp(argv[i], "-help"))
				{
					printf(helpstring);
					return 0;
				}
			}
		}
	}


	//missile definitions

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

	vec missileState = { missileData.V_min,0,missileData.m_0 };

	//target definitions

	sSharedDataTarget targetData;
	targetData.g_t_max = 6;
	targetData.V = 300;

	sPosition targetPosit;
	targetPosit.dPsi = 0;
	targetPosit.pos = { 0,10000 };
	targetPosit.Psi = PI;

	vec targetState = { targetData.V,0 };


	//instantiating world objects

	std::vector<std::shared_ptr<cWorldObject>> objects;
	if (nonlinear)
		objects.push_back(std::make_shared<cMissileNonLinear>(missileState, ts, missileData, missilePosit));
	else
		objects.push_back(std::make_shared<cMissileLTI>(missileState, ts, missileData, missilePosit));
	objects.push_back(std::make_shared<cTarget>(targetState, targetPosit, targetData, ts));

	//setting target for missile object

	objects[0]->SetTarget(objects[1]);

	//initial calculation of stop conditions

	double r = range::GetRange(missilePosit, targetPosit);

	//preparing output file

	FILE* f;
	errno_t err = fopen_s(&f,path.c_str(), "w");
	if (f == nullptr)
	{
		printf("failed to open file");
		return 0;
	}

#ifndef _DEBUG
	fprintf(f,"time; missie velocity; missile X position; missile Y position; target X position; target Y position\n");
	std::vector < std::array< double,6> > results;
	results.reserve(4000);
#else
	fprintf(f,"time,xm,ym,xt,yt,Psim,dPsim,m_u,m_r,Psit,dPsit,t_u,t_r\n");
	std::vector < std::array< double, 13> > results;
	results.reserve(4000);
#endif

	//initializing model time and counter variable for max missile velocity

	double v_max = 0;
	double worldTime = 0;

	//preparing execution time measurments

	auto startTime = high_resolution_clock::now();
	double frameN = 0;
	double totalTime = 0.;

	//main loop

	while (r > r_min && missileState[0]>= missileData.V_min)
	{
		auto frameStart = high_resolution_clock::now();

		//loop through all objects
		for (auto& object : objects)
			object->Simulate();

		//get execution time of simulation step
		auto frameTime = high_resolution_clock::now() - frameStart;
		totalTime += duration_cast<microseconds>(frameTime).count();
		frameN+=1;

		//dump telemetry to buffer
		missilePosit = objects[0]->GetPosition();
		targetPosit = objects[1]->GetPosition();
		missileState = objects[0]->GetState();
#ifndef _DEBUG
		results.push_back({ worldTime, missileState[0], missilePosit.pos[0], missilePosit.pos[1], targetPosit.pos[0], targetPosit.pos[1] });
#else
		targetState = objects[1]->GetState();
		results.push_back({ worldTime, missilePosit.pos[0], missilePosit.pos[1], targetPosit.pos[0], targetPosit.pos[1], missilePosit.Psi, missilePosit.dPsi,
			missileState[0],missileState[1],targetPosit.Psi,targetPosit.dPsi,targetState[0],targetState[1]});
#endif
		
		worldTime += *ts;
		r = range::GetRange(missilePosit, targetPosit);
		v_max = std::max(v_max, missileState[0]);
	}

	//display timing and dump results to file
	double simTime = duration_cast<microseconds>(high_resolution_clock::now()-startTime).count();
	simTime /= 1000;
	double averageFrame = totalTime / frameN;
	printf("Maximum missile speed: %f m/s, Sim time: %f ms, Average step time: %f us, writing file",v_max, simTime, averageFrame);

	for (auto& i : results)
	{
#ifndef _DEBUG
		fprintf(f, "%f;%f;%f;%f;%f;%f\n", i[0], i[1], i[2], i[3], i[4], i[5]);
#else
		fprintf(f, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10], i[11], i[12]);
#endif
	}

	fclose(f);
	
	return 0;

}