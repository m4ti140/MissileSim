### Overview
This program simulates target interception by a missile in 2D plane. The missile is represented with a state space model (by default a linear time-invariant model with parameter lookup at the start of every loop, optionally a full non-linear model is provided) and guided with proportional navigation. The target follows a randomized path. The simulation concludes if either the target is intercepted or missile runs out of energy (missile velocity drops below initial velocity of 300 m/s) The results of each simulation is saved by default in an `out.txt' file in program directory in semi-colon-separated format

### Usage
Called with
`MissileSim target_velocity target_range [options]`
Where `target_velocity` and `target_range` are required parameters initializing target velocity (constant) and initial range respectively.
Following additional options are available
```
	-help - display this help
	-nonlinear - use nonlinear model for the missile
		Default model is linear time invariant state space representation with euler solver and coefficients recalculated every run of main loop
		with -nonlinear full equations of motions and RK4 solver are used
	-ts - set simulation timestep. Must be > 0
	-path - full path to output file. If not provided, out.txt is created in executable's directory
```
`-help` can be called without required parameters

Upon successful execution the program displays maximum velocity achieved by the missile, total execution time (minus file output) and average loop time (minus telemetry output)
Unless otherwised specified using `-path`, the results are dumped into an 'out.txt' file in executable directory. The file uses the following format:
```
time; missie velocity; missile X position; missile Y position; target X position; target Y position
```

### Design assumptions
The following information is known about the missile
  - $C_{D_0} (Ma)$ table
  - Symmetrical drag curve with coefficient of $1.5$
  - Reference area of $0.9 m^2$
  - constant altitude of $2km$ (corresponding to ISA density of $1.00649 kg/m^3$
  - initial velocity of $300m/s$
  - maximium $n_y=30g$ - this corresponds to coefficient of side force of $C_{Y_{max}} \approx 0.168903$ at a velocity of $~809m/s$ which was determined experimentally by running the simulation without guidance. The limit was assumed as an aerodynamic limit rather than a structural one (although the simulation could be further simplified otherwise) as no information was given about $C_{Y_{max}}$
  - Initial mass of $230kg$, propellant mass of $60kg$ and engine $I_sp=235$

Additionally, the following assumptions were taken based on available data
  - Since no information about mass distribution, lift polar or aerodyanmic configuration of the missile is given, the conservation of angular momentum simulation was skipped - it is assumed that the missile always points along velocity vector
  - To facilitate the above the missile configuration is assumed to be wing controlled, with 4 big control surfaces near the center of mass and 4 static fins for stability.
  - Because of these assumptions the local horizon, body and stability frames of referance are alligned and thrust is aligned to velocity vector, simplifying equations of motion in body axes to:
  
  $$m{du \over dt} = I_{sp} {g_0} {dm \over dt} - {1\over2}\rho{u^2}S{C_D}$$
  
  $${dm \over dt} = \dot{m}\\ \\ r = {{1\over2}\rho{u}S{C_Y}\over{m}}$$
  
  - And kinematic equations in world reference frame for piecewise-circular motion can be given by:
  
  $$d\Psi = rdt$$
  
  $$dx = 2{u \over r}sin({{d\Psi} \over 2})sin(\Psi)$$
  
  $$dy = 2{u \over r}sin({{d\Psi} \over 2})cos(\Psi)$$

  - Where:
    * $u = {dx \over dt}$
    * $r = {d{\psi} \over dt}$
    * $C_D = {C_{D_0} (Ma)} + K {C_L}^2 + K {C_Y}^2$ - total drag coefficient
    * $C_L = {2mg_0 \over {\rho}u^2S}$ - required lift coefficient
    * $\dot{m} = m_{propelant}/t_{burn}$ at $t\leq6$
    * $\dot{m} = 0$ at $t>6$
  - Since there is no information about actuators or critical angle of deflection of control surfaces, the $C_Y$ is used directly as a control variable

The state equations can be linearized as

![CodeCogsEqn](https://user-images.githubusercontent.com/19916764/187320403-fe12d9dc-7eb2-48cd-b86c-ff1a89dbfe8d.png)

![CodeCogsEqn2](https://user-images.githubusercontent.com/19916764/187320557-b61ef448-a84d-48c5-85bf-2aa48beefe02.png)

Where:
  - $X_u = -{qSC_D \over mu}$
  - $X_{C_Y} = -{qS \over m}KC_Y$
  - $X_{\dot{m}} = {V_e \over m}$
  - $N_{C_Y} = {qS \over mu}$
  - $q = {1\over2}{\rho}u^2$ - dynamic pressure
  
Finally, the guidance law is given by
$$ C_Y = {m \over qS}k_p\dot{\lambda}V_c
Where
  - $\lambda = atan2({x_{target}-x_{missile} \over y_{target}-y_{missile}})$
  - $V_c = u_{missile}cos(\Psi_{missile}-\lambda)-u_{target}cos(\Psi_{target}-\lambda)$

### Structure
The signal flow for the missile model is represented by the following diagram

![Signal_missile_sim drawio](https://user-images.githubusercontent.com/19916764/187311999-601ecc81-6bea-4310-8224-3a1db65ff076.png)

The signal flow for target model is simpler and represented by:

![signal_target drawio](https://user-images.githubusercontent.com/19916764/187313128-988310ce-38f1-43dd-8709-d71752696a1c.png)

The following diagram shows the class structure. `cMissileNonLinear` is not shown as it is very similar to `cMissileLTI`

![MissileSimstructure](https://user-images.githubusercontent.com/19916764/187315562-73e5ad39-1c42-4e12-9ee4-c8a821000ac9.png)

Notable feature is shared ownership of data structures (state space matrices, input/output/observation vectors etc.) between parent `cWorldObject` and all worker objects it is composed of, to avoid passing references on every call of individual worker methods. Dynamic polymorphism is exclusively used between `cWorldObjects` to avoid unnecessary virtual function calls while still allowing for an arbitrary amount of missiles and targets to be simulated if need be.

Additionally, `cAzimuth` class (see utilities.h) was created to handle operations on angles.

`cLTI` uses Euler solver to integrate state equations, while `cGeneralSSM` uses Runge-Kutta 4th order solver and is supplied with arbitrary functions $f(x,u,t,...)$ to allow for full equations of motion without linearization.

In debug environment there is a 50% performance loss when using non-linear solver for minimal gain in quality. With more complex systems howeer, pre-generated lookup tables would likely be used for LTI system, which would further increase performance over full non-linear model.
