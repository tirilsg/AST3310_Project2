# AST3310_Project2

This repository contains python code constructing a model for energy transportation within a sun `EnergyTransportation` in `project2.py`, making use of the model for energy production defined in `stellar_core.py`, the class `StellarCore`. The code can be ran by just calling `project2.py`, as it imports all important functions in the documents `atellar_core.py`, `cross_section.py` and `cross_section_sub.py` as modules. 


-------------------------------------

### `stellar_core.py`: 
Contains a class `StellarCore` which takes a density constant $\rho$ and temperature T, and estimates the energy production at an arbitrary, fixed time. This energy production is a result of fusion of hydrogen into helium, through four different cycles of fusion reactions PP1, PP2, PP3 and CNO, derived and explained in the report. The class contains the following methods:

* `einsteins_principle(m_0,m_1):` which estimates energy from the einstein's mass-energy equivalency principle

* `number_density(mass_frac, Z):` estimates the number density of an arbitrary element Z

* `calculate_densities():` calls on the function `number_density` to calculate the number density of the elements relevant to the fusion reaction taking place within the core. The densities and isotopes are coupled and stored within a directory `n` so that the information can be extracted whenever

* `proportionality():` defines the relevant proportionality functions, for each of the relevant reactions of fusion. The proportionality is saved within a directory
  
* `reaction_rates(alter):` calculates the rates of reactions for of the relevant reactions of fusion from the calculated number densities and proportionality functions calculated and stored by `number_density()` and `proportionality()`. Takes 'alter' as an argument, which if alter=True is grounds for an alteration to the rate-calculation which takes into account realistic element production and consumption.

* `estimated_energy():` estimates the energy output of each reaction of fusion taking place within the stellar core from `einsteins_principle()`

* `cycle_relation():` calculates the total energy outputs from the reaction rates and energies calculated by `reaction_rates()` and `estimated_energy()`, and relates them to each of the cycles of fusion

* `sanity_check(t_c):` performs a test of the correctness of the model, for the temperatures $T=10^8$ or $T=1.57\times10^7$, and prints to the terminal if the sanity check is failed to be met 

---------------------------------------
### `project2.py`: 
Contains a redefinition of the class defined in `stellar_core.py`, which adds a calculation of the total energy-production to the function relating the production of energy to each branch of fusion occurring. 

Furthermore, a class `EnergyTransportation(L0, R0, M0, rho0, T0, print_sanity)` is implemented, which takes the surface defining variables L0, R0, M0, rho0, T0, and the boolean print_sanity as arguments, and defines the contents; 

* `initialized():` takes the surface-defining variables and defines L, r, m, rho and T, as well as significant constants nabla_ad, mean_molecular, the heat capacity cp, pressure P, as well as all variables calculated by the function `calculate_variables()`. If the boolean print_sanity is True, the function also prints the temperature gradients, as well as the mean molecular weight to the terminal in order to perform sanity checks.

  
* `read_opacity():` opens and reads the file `opacity.txt`, and creates a function for opacity by use of interpolation from the data within the txt-file.
  
* `opacity(rho,T):` calculates R from rho, and uses the function created by `read_opacity()` to return a value for opacity

* `calculate_variables():` calls on all functions, to calculate the variables pressure_rad, pressure_gas, rho, energies, K, Hp, lm, g, U, nabla_stable, nabla_star, nabla_p, xi, F_con and F_rad
  
* `dr_dm(m,r):` takes mass m and radius r as arguments, calculates the governing equation, representing radius change with respect to mass. 
  
* `dP_dm(m,P):` takes mass m and pressure P as arguments, calculates the governing equation, representing pressure change with respect to mass.
  
* `dL_dm(m,L):` takes mass m and luminosity L as arguments, calculates the governing equation, representing luminosity change with respect to mass.
  
* `dT_dm(m,T):` takes mass m and temperature T as arguments, calculates the governing equation, representing temperature change with respect to mass. 
  
* `calculate_Hp(P,m,r):` takes pressure P, mass m and radius r as arguments, and calculates the pressure scale height $H_P$
  
* `density(P,T):` calculates the density from the arguments pressure P and temperature T
  
* `calculate_g(M,R):` calculates the gravitational acceleration constant g for the star, from the arguments mass M and radius R
  
* `calculate_U(T, K, rho, Hp,g):` calculates the argument U from the arguments temperature T, opacity K, density $\rho$, pressure scale height $H_P$ and the gravitational acceleration constant g
  
* `nbl_star(nbl_stbl,lm,xi,U):` calculates the temperature gradient $\nabla^*$ 
  
* `nbl_stable(L,rho,T,K,r,Hp):` calculates the temperature gradient $\nabla_{stable}$ 
  
* `nbl_p(nbl_str,xi):` calculates the temperature gradient $\nabla_{p}$ 
  
* `nbl_xi(lm,U,nbl_ad,nbl_stbl):` calculates the temperature gradient $\xi$ 
  
* `F_rad_cal(T,K,rho,Hp,nbl_str):` calculates the flux of energy due to energy transportation by radiation 
  
* `F_con_cal(T,g,rho,lm,Hp,xi,nbl_ad,P,nbl_stbl):` calculates the flux of energy due to energy transportation by convection
  
* `variable_step(p):` calculates the step length by the variables r, P, L and T, divided by their respective governing equations, scaled by an arbitrary p which can be adjusted from need. Returns the smallest absolute value as the steplength dm.
  
* `step(dm):` evolves the system by a step of length dm, and calculates all nessecary variables for further simulation by calling `calculate_variables()`.
  
* `run_simulation(init,p=0.01):` runs the simulation and stores data in lists later turned into arrays, while the mass and radius stays positive.
  
* `opacity_sanity(tol):` performs a sanity check for the opacity function 
  
* `print_expected():` prints a number of values in order to perform a further sanity check 

Finally, a bunch of functions are implemented, used to create temperature gradient-radius plots by `gradient_plot(r,nabla_star, nabla_stable,nabla_ad,p,name="name")`, energy-radius plots by `energies_plot(r,Q_PP1_values, Q_PP2_values,Q_PP3_values,Q_CNO_values,Q_total,p,name="name")`, flux-radius plots `flux_plot(r,F_con, F_rad,p,name="name")`, variable-radius plots `plot_all_variables(m,r,rho,P,L,T,p,name="name")`, and scaled variable-radius plots `plot_all_variables(m,r,rho,P,L,T,p,name="name")`. A function that plots the structural implications of changing the surface defining variables is also `implemented simulation_parameters(p):`, and all figures are saved by the code in the directory "report/Figures/". At the end, all these functions are called by testing the parameters simulation_parameters(p) for p=0.1, testing the sanity of the opacity function for a class defined by the variables `simulation=EnergyTransportation(1,1,1,1.42*10**(-7),5770, print_sanity=True)`, as well as cross section plots and temperature plots calling their respective functions. The sanity of the model itself are further tested by `simulation=RedefOpacity(L0=1,R0=0.84,M0=0.99, rho0=55.9/(1.408*10**3 ),T0= 0.9*10**6, print_sanity=True)`, and lastly all plots are created by calling their respective functions for a class defined by `L0=1.7,R0=1,M0=1.1, rho0=70*1.42*10**(-7),T0= 2*5770`. 


---------------------------------------------
### `cross_section.py`: 
Defines a function `cross_section(R, L, F_C,cr,cv, sanity=False, savefig=False, name="name")` that creates cross-section plots in its own figure 


---------------------------------------------
### `cross_section_sub.py`: 
Defines a function `cross_section_mod(R, L, F_C, ax=None, sanity=False, savefig=False, name="name")` that creates cross-section plots which can be added to an arbitrary space in a grid of subplots

