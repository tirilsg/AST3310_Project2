#program that implements a class simulating energy transportation through a star, starting from the surface moving towards the core
#the class is used to run simulations with varying surface defining variables and performing analysis on data and the structural implications
#creates a final model, creating a star with preferable transportation-shape based on the testing, limiting the size of a radiative energy transportation zone at the surface of the star
#creates temperature gradient, energy-production, transportation zone, energy flux and varable behaviour plots for the model using these surface defining parameters

#to run this python program, opacity.txt, the files stellar_core.py, cross_section_sub.py and cross_section.py must reside within the same map

from stellar_core import StellarCore #imports class simulating energy production 
from cross_section_sub import cross_section_mod #imports a function creating cross_section plots that can be implemented in subplots
from cross_section import cross_section #cross_section function given with the project desctiption, containing small changes
import scipy.interpolate as intp #interpolation function
import matplotlib.pyplot as plt #plotting module
import numpy as np 

# a couple changes needed to be done to the energy production calculator in the StellarCore class in order to make accurate calculations
# this is done by making the changes directly in this python doc
class Energy(StellarCore):
    def cycle_relation(self): #relation between each fusion-cycle and its energy output calculated by Equation 6 in the report
        self.E_PP0 = (self.E_l1 + self.E_l1d)
        self.E_CNO = self.E_lp12 + self.E_l13 + self.E_lp13 + self.E_lp14 + self.E_lp15 + self.E_l15
        self.Q_PP0 = self.r_lpp * (self.E_PP0)
        self.Q_PP1 = self.r_l33 * (self.E_l33 + 2*self.E_PP0) 
        self.Q_PP2 = self.r_l34 * (self.E_l34 + self.E_PP0) + self.r_le7*self.E_le7 + self.r_l17l*self.E_l17l
        self.Q_PP3 = self.r_l17 * (self.E_l17 + self.E_dc) + self.r_l34 * (self.E_l34 + self.E_PP0)
        self.Q_CNO = self.r_lp14 * self.E_CNO
        self.tot   = self.r_lpp * self.E_PP0 + self.r_l33 * self.E_l33 + self.r_l34 * self.E_l34+  self.r_le7 * self.E_le7 \
        + self.r_l17l * self.E_l17l + self.r_l17 * (self.E_l17 + self.E_dc) + self.r_lp14 * (self.E_CNO)   


class EnergyTransportation:
    def __init__(self,L0, R0, M0, rho0, T0, print_sanity):
        self.L0=L0
        self.R0=R0
        self.M0=M0
        self.rho0=rho0
        self.T0=T0
        
        self.print_sanity=print_sanity
        
        #variables defining the surface of the sun
        self.lum_sun=3.846*10**26       # units of W, solar luminosity
        self.rad_sun=6.96*10**8         # units meters, solar radius
        self.mass_sun=1.989*10**30      # units of kg, solar mass
        self.avrg_dens_sun=1.408*10**3  # units of kg/m^3, average density of the sun
        
        #mass fractions needed for calculations of molecular mass
        self.X=0.7 ; self.Y3He=10**(-10) ; self.Y=0.29
        self.Z7Li=10**(-7) ; self.Z7Be=10**(-7) ; self.Z14N=10**(-11)
        
        #definitions of constants and units
        self.u=1.660539*10**(-27)    # unit kg, atomic mass unit u convertion
        self.kb=1.3806*10**(-23)     # unit m^2kg/(s^2K), Boltzmann's constant 
        self.stfbltz=5.6704*10**(-8) # unit of W/(m^2K^4), Stefan Boltzmann's constant
        self.G=6.6742*10**(-11)      # unit of Nm^2/kg^2, gravitational constant for the sun
        self.c=2.9979*10**8          # unit of m/s, speed of light
        
        #defining the function needed calculating the opacity from density and temperature
        self.read_opacity() 
        
        #calculating each variable needed for calculation in this initial condition for the system snapshot
        self.initialize()
    
    #function initializing the system
    def initialize(self): 
        #calculating the initial values of L, r, m, rho and T
        self.L = self.L0*self.lum_sun
        self.r = self.R0*self.rad_sun
        self.m = self.M0*self.mass_sun
        self.rho = self.rho0*self.avrg_dens_sun
        self.T = self.T0
        
        #defining the adiabatic temperature gradient
        self.nabla_ad=2/5
        
        #calculating the mean molecular weight, equation 1 in the report
        self.mean_molecular = 1/(2*self.X+self.Y3He+3*self.Y/4+4*self.Z7Li/7+5*self.Z7Be/7+8*self.Z14N/14)
        
        #initial pressure calculations
        self.pressure_gas = self.kb*self.T*self.rho/(self.u*self.mean_molecular) #equation 3 in the report
        self.pressure_rad = 4*self.stfbltz*self.T**4/(3*self.c) #equation 4 in the report
        self.P = self.pressure_gas + self.pressure_rad #equation 5 in the report
        
        #defining the heat capacity used for other calculations
        self.cp = self.kb/(self.mean_molecular*self.u)/self.nabla_ad #equation 20 in the report
        
        #defining every variable needed for further calculations in this initial condition 
        self.calculate_variables()
        
        #prints temperature gradients, and meal molecular weight
        if self.print_sanity==True:
            print(f"Parameter test for L0={self.L0}, R0={self.R0}, M0={self.M0}, T0={self.T0}, rho0={self.rho0:.3e}:")
            print(f"Mean molecular weight: {self.mean_molecular:.3f}")
            print(f"nabla_ad             : {self.nabla_ad}")
            print(f"nabla_p              : {self.nabla_p}")
            print(f"nabla^*              : {self.nabla_star}")
            print(f"nabla_stable         : {self.nabla_stable}")
            print(f"")
        
        
    #function opening the file containing R, T and K data, creating an interpolation of the data and creating function
    def read_opacity(self): 
        file=np.genfromtxt("opacity.txt")
        #saving the data
        self.logT=np.array(file[1:,0])
        self.logR=np.array(file[0,1:])
        self.logK=np.array(file[1:,1:])
        #creating function from interpolation
        self.opacityfunc = intp.interp2d(self.logR,self.logT,self.logK)
        self.warning = False
        return self.opacityfunc #returning the function 

    #function that takes an arbitry value of density and temperature and returns an estimated opacity
    def opacity(self,rho,T): 
        logT = np.log10(T) #log-value for T
        logR = np.log10(rho/10**3/(T*1e-6)**3) #takes the density an transforms it into the log-value for R
        #uses the interpolation function to calculate opacity and returns to non-logarithmic function
        K = 10**(self.opacityfunc(logR,logT)[0])/10 
        return K
    
    #function that calculates every relevant variable from L,R,M,rho, T and P
    def calculate_variables(self): 
        self.pressure_rad = 4*self.stfbltz/(3*self.c)*self.T**4                        #pressure from radiation
        self.pressure_gas = self.P - self.pressure_rad                                 #pressure from gas
        self.rho          = self.density(self.P,self.T)                                #calculates density from pressure and temperature
        self.energies     = Energy(self.rho,self.T)                                    #creates an instance of the energy-production class
        self.K            = self.opacity(self.rho,self.T)                              #calculates the opacity from density and temperature
        self.Hp           = self.calculate_Hp(self.P,self.m,self.r)                    #calculates pressure scale height
        self.lm           = self.Hp                                                    #defines mixing lenght
        self.g            = self.calculate_g(self.m,self.r)                            #calculates g
        self.U            = self.calculate_U(self.T, self.K, self.rho, self.Hp,self.g) #calculates U

        #temperature gradient calculations 
        #self.nabla_stable = 3*self.L*self.K*self.rho*self.Hp/(64*np.pi*self.r**2*self.stfbltz*self.T**4)
        self.nabla_stable= self.nbl_stable(self.L,self.rho,self.T,self.K,self.r,self.Hp)
        self.xi           = self.nbl_xi(self.lm,self.U,self.nabla_ad,self.nabla_stable)
        self.nabla_star   = self.nbl_star(self.nabla_stable,self.lm,self.xi,self.U)
        self.nabla_p      = self.nbl_p(self.nabla_star,self.xi)
        if self.nabla_stable > self.nabla_ad:
            self.nabla_star += 0
        else:
            self.nabla_star = self.nabla_stable
        #energy flux calculations
        self.F_con        = self.F_con_cal(self.T,self.g,self.rho,self.lm,self.Hp,self.xi,self.nabla_ad,self.P,self.nabla_stable)
        self.F_rad        = self.F_rad_cal(self.T,self.K,self.rho,self.Hp,self.nabla_star)
    
    
    #################################################################################
    ##################### definition of the governing equations #####################
    #################################################################################
    
    #equation 6 in the report
    def dr_dm(self,m,r): 
        return 1/(4*np.pi*r**2*self.rho)

    #equation 7 in the report   
    def dP_dm(self,m,P): 
        return -self.G*m/(4*np.pi*self.r**4)

    #equation 8 in the report
    def dL_dm(self,m,L):
        return self.energies.tot

    def dT_dm(self,m,T):
        if self.nabla_stable > self.nabla_ad: #checks what definition of dT_dm we can use
            #equation 9 in the report
            dTdm = self.nabla_star*(T/self.P)*self.dP_dm(m,self.P)
        else:
            #equation 6 in the report
            dTdm = -3*self.K*self.L/(256*np.pi**2*self.stfbltz*self.r**4*T**3)
        return dTdm
    
    
    ################################################################################
    #################### definition of other relevant equations ####################
    ################################################################################

    #equation 15 in the report
    def calculate_Hp(self,P,m,r):
        return -P*(self.dr_dm(m,r)/self.dP_dm(m,P)) 
    
    def density(self,P,T):
        return self.u*self.mean_molecular/(self.kb*T)*(P-4*self.stfbltz*T**4/(3*self.c))
    
    #gravitational acceleration for the star
    def calculate_g(self,M,R):
        return self.G*M/R**2 
    
    #introduced in the report right before equation 26
    def calculate_U(self, T, K, rho, Hp,g):
        return 64*self.stfbltz*T**3/(3*K*rho**2*self.cp)*np.sqrt(Hp/g)
    

    ################################################################################
    ##################### calculation of temperature gradients #####################
    ################################################################################
    
    #equation 30 in the report
    def nbl_star(self,nbl_stbl,lm,xi,U):
        return -lm**2*xi**3/U+nbl_stbl
    
    #equation 29 in the report
    def nbl_stable(self,L,rho,T,K,r,Hp):
        return L*3*K*rho*Hp/(64*np.pi*r**2*self.stfbltz*T**4)
    
    #equation 31 in the report
    def nbl_p(self,nbl_str,xi):
        return nbl_str-xi**2
    
    #solves the third degree polynomial that is xi defined by equation 28 in the report
    def nbl_xi(self,lm,U,nbl_ad,nbl_stbl):
        a=lm**2/U
        b=1
        c=U*(4/lm**2)
        d=-(nbl_stbl-nbl_ad)
        roots = np.roots([a,b,c,d])
        xi=roots[np.argwhere(roots.imag==0)].real[0][0] #extracts the real root
        return xi #returns the value for xi 
    
    
    ###############################################################################
    ######################### calculation of energy flux ##########################
    ###############################################################################
    
    #equation 14 in the report
    def F_rad_cal(self,T,K,rho,Hp,nbl_str):
        return 16*self.stfbltz*T**4/(3*K*rho*Hp)*nbl_str
    
    #equation 23 in the report
    def F_con_cal(self,T,g,rho,lm,Hp,xi,nbl_ad,P,nbl_stbl):
        if nbl_stbl > nbl_ad:
            val = rho*self.cp*T*np.sqrt(g/Hp**3)*(lm/2)**2*xi**3
        else:
            val = 0
        return val
    
    #implements the variable step length, defined by p
    def variable_step(self,p):
        dm_r = p*self.r/self.dr_dm(self.m,self.r)
        dm_P = p*self.P/self.dP_dm(self.m,self.P)
        dm_L = p*self.L/self.dL_dm(self.m,self.L)
        dm_T = p*self.T/self.dT_dm(self.m,self.T)
        dm = np.min(abs(np.array([dm_r,dm_P,dm_L,dm_T])))
        return dm

    #evolves the system a step dm, and calculates all the new variables, and defines them as variables in the class
    def step(self,dm):
        #uses the governing equations 
        #subtracts the change in radius from the previous r
        self.r = self.r - self.dr_dm(self.m, self.r)*dm
        #subtracts the change in pressure from the previous pressure calculation
        self.P = self.P - self.dP_dm(self.m, self.P)*dm 
        #subtracts the change in luminosity from the previous calculation
        self.L = self.L - self.dL_dm(self.m, self.L)*dm
        #subtracts the change in temperature from the previous temperature definition
        self.T = self.T - self.dT_dm(self.m, self.T)*dm
        #estimates the new mass 
        self.m = self.m - dm
        self.calculate_variables() #calculates all new variables
    
    #runs the simulation for steps of variable length, until the mass hits 0, and stores all relevant data
    def run_simulation(self,init,p=0.01):
        #defines empty lists
        Q_PP1_values=[]       ;  Q_PP2_values=[]         ;  Q_PP3_values=[]          
        Q_CNO_values=[]       ; Q_total=[]
        nabla_star=[]  ;  nabla_stable=[]  ;  nabla_ad=[] 
        m_values=[]    ;  r_values=[]      ;  rho_values=[]    ;  P_values=[] 
        T_values=[]    ;  L_values=[]      ;  F_con_values=[]  ;  F_rad_values=[] 
        core= False
        while init.m > 0 and init.r > 0:
            #appends all relevant values to the relevant lists
            m_values.append(init.m)
            r_values.append(init.r)
            rho_values.append(init.rho)
            P_values.append(init.P)
            L_values.append(init.L)
            T_values.append(init.T)
            F_con_values.append(init.F_con)
            F_rad_values.append(init.F_rad)
            Q_total.append(init.energies.tot)
            Q_PP1_values.append(init.energies.Q_PP1)
            Q_PP2_values.append(init.energies.Q_PP2)
            Q_PP3_values.append(init.energies.Q_PP3)
            Q_CNO_values.append(init.energies.Q_CNO)
            nabla_ad.append(init.nabla_ad)
            nabla_star.append(init.nabla_star)
            nabla_stable.append(init.nabla_stable)
            #calculates the new step length
            dm = init.variable_step(p)
            #evolves the system a step by using governing equations 
            init.step(dm)
            if L_values[-1]>0.995*L_values[0]:
                self.core_size=init.r/r_values[0]
            if F_con_values[-1]!=0 and init.F_con==0 and core==False:
                self.conv_zone=((r_values[0]-init.r)/r_values[0])
                core = True
            if init.r==r_values[-1]:
                break
        Q_PP1_values=np.array(Q_PP1_values)     ;  Q_PP2_values=np.array(Q_PP2_values)     ;  Q_PP3_values=np.array(Q_PP3_values)          
        Q_CNO_values=np.array(Q_CNO_values)       ; Q_total=np.array(Q_total) 
        nabla_star=np.array(nabla_star)   ;  nabla_stable=np.array(nabla_stable)  ;  nabla_ad=np.array(nabla_ad) 
        m_values=np.array(m_values)    ;  r_values=np.array(r_values)       ;  rho_values=np.array(rho_values)     ;  P_values=np.array(P_values)  
        T_values=np.array(T_values)     ;  L_values=np.array(L_values)       ;  F_con_values=np.array(F_con_values)   ;  F_rad_values=np.array(F_rad_values) 
        return Q_PP1_values, Q_PP2_values,Q_PP3_values,Q_CNO_values,Q_total,nabla_star, nabla_stable, nabla_ad,m_values,r_values,rho_values,P_values,T_values,L_values,F_con_values,F_rad_values

    
    #performs a sanity check for the opacity estimation 
    def opacity_sanity(self, tol):
        #definitions of expectation values
        log_10_T=np.array([3.750, 3.755, 3.755, 3.755, 3.755, 3.770, 3.780, 3.795, 3.770, 3.775, 3.780, 3.795, 3.800])
        log_10_R=np.array([-6.00, -5.95, -5.80, -5.70, -5.55, -5.95, -5.95 , -5.95, -5.80, -5.75, -5.70, -5.55, -5.50])
        log_10_k=np.array([2.84*10**(-3), 3.11*10**(-3), 2.68*10**(-3), 2.46*10**(-3), 2.12*10**(-3), 4.70*10**(-3), 6.25*10**(-3), 
                           9.45*10**(-3), 4.05*10**(-3), 4.43*10**(-3), 4.94*10**(-3), 6.89*10**(-3), 7.69*10**(-3)])
        T = 10**log_10_T
        rho = (T/10**(6))**3*10**log_10_R*10**3
        print("Sanity check for Opacity estimation")
        print("Expected K [SI]:   Estimated K [SI]:    Relative Error:")
        for i in range(len(log_10_T)):
            K = self.opacity(rho[i],T[i])
            print(f"    {log_10_k[i]:.3e}          {K:.3e}           {abs(K-log_10_k[i]):.3e}")
            if abs(K-log_10_k[i])>tol:
                print("The sanity check for opacity fails")
    
    #function that prints values for a check of correctness if needed
    def print_expected(self):
        print(f"")
        print(f"The values definig behaviour L0={self.L0}, R0={self.R0}, M0={self.M0}, T0={self.T0}, rho0={self.rho0:.3e}:")
        print(f"Hp                       :  {self.Hp:.3f}")
        print(f"U                        :  {self.U:.3f}")
        print(f"Opacity                  :  {self.K:.3f}")
        print(f"xi                       :  {self.xi:.3e}")
        print(f"F_con/(F_con+F_rad)      :  {self.F_con/(self.F_con+self.F_rad):.3f}")
        print(f"F_rad/(F_con+F_rad)      :  {self.F_rad/(self.F_con+self.F_rad):.3f}")
        print(f"{self.nabla_ad:10f} < {self.nabla_p:.10f} < {self.nabla_star:.10f} < {self.nabla_stable:.10f}")
        print(f"")


#############################################################################################
########################## definition of functions creating plots ###########################
#############################################################################################
def gradient_plot(r,nabla_star, nabla_stable,nabla_ad,p,name="name"):
    plt.plot(r/r[0],nabla_stable,label = r"$\nabla_{stable}$")
    plt.plot(r/r[0],nabla_star,label = r"$\nabla^*$")
    plt.plot(r/r[0],nabla_ad, label = r"$\nabla_{ad}$")
    plt.yscale("log")
    plt.xlabel(r"$R/R_\odot$")
    plt.ylabel(r"$\nabla$")
    plt.title(f"Temperature gradient plot for p={p}")
    plt.legend()
    plt.savefig(f"report/Figures/TemperatureGradients{name}.pdf")
    plt.show()

def energies_plot(r,Q_PP1_values, Q_PP2_values,Q_PP3_values,Q_CNO_values,Q_total,p,name="name"):
    scaling=Q_PP1_values+Q_PP2_values+Q_PP3_values+Q_CNO_values
    plt.plot(r/r[0], Q_PP1_values/scaling ,label = r"$\epsilon_{PP1}$")
    plt.plot(r/r[0],Q_PP2_values/scaling  ,label = r"$\epsilon_{PP2}$")
    plt.plot(r/r[0],Q_PP3_values/scaling   , label = r"$\epsilon_{PP3}$") 
    plt.plot(r/r[0],Q_CNO_values/scaling  , label = r"$\epsilon_{CNO}$") 
    plt.plot(r/r[0], Q_total/np.max(Q_total) , label = r"$\epsilon/\epsilon_{total}$")
    plt.xlabel(r"$R/R_\odot$")
    plt.ylabel("Scaled Energy")
    plt.title(f"Energy scaled by total production for p={p}")
    plt.legend()
    plt.savefig(f"report/Figures/EnergyProduction{name}.pdf")
    plt.show()

def flux_plot(r,F_con, F_rad,p,name="name"):
    flux = F_con + F_rad
    plt.plot(r/r[0], F_con/flux , label = r"Convection Flux, $F_{con}$") 
    plt.plot(r/r[0], F_rad/flux , label = r"Radiation Flux, $F_{rad}$") 
    plt.xlabel(r"$R/R_\odot$")
    plt.ylabel("Flux")
    plt.title(f"Flux of energy scaled by total flux, for p={p}")
    plt.legend()
    plt.savefig(f"report/Figures/Flux{name}.pdf")
    plt.show()

def plot_all_variables(m,r,rho,P,L,T,p,name="name"):
    plt.plot(r/r[0], L/np.max(L),label = r"Luminosity L") 
    plt.plot(r/r[0], P/np.max(P),label = r"Pressure P")
    plt.plot(r/r[0], T/np.max(T),label = r"Temperature T")
    plt.plot(r/r[0], rho/np.max(rho),label = r"Density $\rho$")
    plt.plot(r/r[0], m/np.max(m),label = r"Mass m")
    plt.xlabel(r"$R/R_\odot$")
    plt.ylabel("Fraction of Max Value")
    plt.title(f"All varriables scaled, for p={p}")
    plt.legend()
    plt.savefig(f"report/Figures/Variables{name}.pdf")
    plt.show()

def variables_special_scaling(m,r,rho,P,L,T,p,name="name"):
    fig,axs = plt.subplots(2,2, figsize=(12, 6))
    axs[0,0].plot(r/r[0], L/L[0],label = r"$L/L_\odot$")
    axs[0,0].plot(r/r[0], m/m[0],label = r"$m/m_\odot$")  
    axs[0,0].set_ylabel("Fraction of Max Value")
    axs[0,0].legend()
    axs[1,0].plot(r/r[0], P ,label = r"Pressure P")
    axs[1,0].set_xlabel(r"$R/R_\odot$")
    axs[1,0].set_ylabel("Fraction of Max Value")
    axs[1,0].set_yscale("log")
    axs[1,0].set_title(f"Pressure P")
    axs[0,1].plot(r/r[0], T,label = r"Temperature T")
    axs[0,1].set_title(f"Temperature T")
    axs[1,1].plot(r/r[0], rho/rho[0],label = r"$\rho /\rho_\odot$") 
    axs[1,1].set_xlabel(r"$R/R_\odot$")
    axs[1,1].set_yscale("log")
    axs[1,1].set_title(r"$\rho /\rho_\odot$")
    plt.savefig(f"report/Figures/VariablesSpecial{name}.pdf")
    plt.show()

#function that simulations different variables and plots their impact on the energy transportation zone structure
def simulation_parameters(p):
    L0=1 ; R0=1 ; T0=5770
    M0=1 ; rho0= 1.42*10**(-7)
    ############################## DENSITY #############################
    fig,axs = plt.subplots(4,2, figsize=(10, 20), layout = "tight")
    axs[0,0].set_title(r"$\rho_0/10$")
    simulation=EnergyTransportation(L0, R0  , M0, rho0/10, T0,print_sanity=False)
    Q_PP1_values, Q_PP2_values,Q_PP3_values,Q_CNO_values,Q_total,nabla_star, nabla_stable, nabla_ad,m_values,r_values1,rho_values,P_values,T_values,L_values,F_con_values,F_rad_values=simulation.run_simulation(simulation,p)
    star_zone1=cross_section_mod(r_values1,L_values,F_con_values, ax=axs[0,0])
    axs[0,0].text(0,0.9,f"Core size: {simulation.core_size:.3f}")
    axs[0,0].text(0,0.5,f"Outer Conv zone: {simulation.conv_zone:.3f}")
    simulation=EnergyTransportation(L0, R0 , M0, rho0*10, T0,print_sanity=False) 
    axs[0,1].set_title(r"$\rho_0\times10$")
    Q_PP1_values, Q_PP2_values,Q_PP3_values,Q_CNO_values,Q_total,nabla_star, nabla_stable, nabla_ad,m_values,r_values2,rho_values,P_values,T_values,L_values,F_con_values,F_rad_values=simulation.run_simulation(simulation,p)
    star_zone2=cross_section_mod(r_values2,L_values,F_con_values, ax=axs[0,1])
    axs[0,1].text(0,0.9,f"Core size: {simulation.core_size:.3f}")
    axs[0,1].text(0,0.5,f"Outer Conv zone: {simulation.conv_zone:.3f}")
    ############################## TEMPERATURE #############################
    axs[1,0].set_title(r"$T_0/2$")
    simulation=EnergyTransportation(L0  , R0  , M0, rho0, T0/2,print_sanity=False)
    Q_PP1_values, Q_PP2_values,Q_PP3_values,Q_CNO_values,Q_total,nabla_star, nabla_stable, nabla_ad,m_values,r_values3,rho_values,P_values,T_values,L_values,F_con_values,F_rad_values=simulation.run_simulation(simulation,p)
    star_zone3=cross_section_mod(r_values3,L_values,F_con_values, ax=axs[1,0])
    axs[1,0].set_title(r"$T_0/2$")
    axs[1,0].text(0,0.9,f"Core size: {simulation.core_size:.3f}")
    axs[1,0].text(0,0.5,f"Outer Conv zone: {simulation.conv_zone:.3f}")
    simulation=EnergyTransportation(L0  , R0  , M0, rho0, T0*2,print_sanity=False) 
    Q_PP1_values, Q_PP2_values,Q_PP3_values,Q_CNO_values,Q_total,nabla_star, nabla_stable, nabla_ad,m_values,r_values4,rho_values,P_values,T_values,L_values,F_con_values,F_rad_values=simulation.run_simulation(simulation,p)
    star_zone4=cross_section_mod(r_values4,L_values,F_con_values, ax=axs[1,1])
    axs[1,1].text(0,0.9,f"Core size: {simulation.core_size:.3f}")
    axs[1,1].text(0,0.5,f"Outer Conv zone: {simulation.conv_zone:.3f}")
    axs[1,1].set_title(r"$T_0\times2$")
    
    ############################# LUMINOCITY #############################    
    axs[2,0].set_title(r"$L_0/2$")
    simulation=EnergyTransportation(L0/2 , R0  , M0, rho0, T0,print_sanity=False)
    Q_PP1_values, Q_PP2_values,Q_PP3_values,Q_CNO_values,Q_total,nabla_star, nabla_stable, nabla_ad,m_values,r_values5,rho_values,P_values,T_values,L_values,F_con_values,F_rad_values=simulation.run_simulation(simulation,p)
    star_zone5=cross_section_mod(r_values5,L_values,F_con_values, ax=axs[2,0])
    axs[2,0].text(0,0.9,f"Core size: {simulation.core_size:.3f}")
    axs[2,0].text(0,0.5,f"Outer Conv zone: {simulation.conv_zone:.3f}")
    axs[2,1].set_title(r"$L_0$")
    simulation=EnergyTransportation(L0, R0  , M0, rho0, T0,print_sanity=False) 
    Q_PP1_values, Q_PP2_values,Q_PP3_values,Q_CNO_values,Q_total,nabla_star, nabla_stable, nabla_ad,m_values,r_values6,rho_values,P_values,T_values,L_values,F_con_values,F_rad_values=simulation.run_simulation(simulation,p)
    star_zone6=cross_section_mod(r_values6,L_values,F_con_values, ax=axs[2,1])
    axs[2,1].text(0,0.9,f"Core size: {simulation.core_size:.3f}")
    axs[2,1].text(0,0.5,f"Outer Conv zone: {simulation.conv_zone:.3f}")
   
    ############################## RADIUS #############################
    axs[3,0].set_title(r"$R_0$")
    simulation=EnergyTransportation(L0, R0, M0, rho0, T0,print_sanity=False)
    Q_PP1_values, Q_PP2_values,Q_PP3_values,Q_CNO_values,Q_total,nabla_star, nabla_stable, nabla_ad,m_values,r_values7,rho_values,P_values,T_values,L_values,F_con_values,F_rad_values=simulation.run_simulation(simulation,p)
    star_zone7=cross_section_mod(r_values7,L_values,F_con_values, ax=axs[3,0])
    axs[3,0].text(0,0.9,f"Core size: {simulation.core_size:.3f}")
    axs[3,0].text(0,0.5,f"Outer Conv zone: {simulation.conv_zone:.3f}")
    axs[3,1].set_title(r"$R_0\times2$")
    simulation=EnergyTransportation(L0 , R0*2, M0, rho0, T0,print_sanity=False) 
    Q_PP1_values, Q_PP2_values,Q_PP3_values,Q_CNO_values,Q_total,nabla_star, nabla_stable, nabla_ad,m_values,r_values8,rho_values,P_values,T_values,L_values,F_con_values,F_rad_values=simulation.run_simulation(simulation,p)
    star_zone8=cross_section_mod(r_values8,L_values,F_con_values, ax=axs[3,1])
    axs[3,1].text(0,0.9,f"Core size: {simulation.core_size:.3f}")
    axs[3,1].text(0,0.5,f"Outer Conv zone: {simulation.conv_zone:.3f}")
    plt.savefig("report/Figures/simulation_part1.pdf")
    
    ############################## PRESSURE #############################
    fig,axs = plt.subplots(1,2, figsize=(10, 5), layout = "tight")
    axs[0].set_title(r"$P_0/10$")
    simulation=EnergyTransportation(L0, R0, M0, rho0, T0,print_sanity=False)
    simulation.P=simulation.P/10
    Q_PP1_values, Q_PP2_values,Q_PP3_values,Q_CNO_values,Q_total,nabla_star, nabla_stable, nabla_ad,m_values,r_values9,rho_values,P_values,T_values,L_values,F_con_values,F_rad_values=simulation.run_simulation(simulation,p)
    star_zone9=cross_section_mod(r_values9,L_values,F_con_values, ax=axs[0])
    axs[0].text(0,0.9,f"Core size: {simulation.core_size:.3f}")
    axs[0].text(0,0.5,f"Outer Conv zone: {simulation.conv_zone:.3f}")
    axs[1].set_title(r"$P\times10$")
    simulation=EnergyTransportation(L0 , R0, M0, rho0, T0,print_sanity=False) 
    simulation.P=simulation.P*10
    Q_PP1_values, Q_PP2_values,Q_PP3_values,Q_CNO_values,Q_total,nabla_star, nabla_stable, nabla_ad,m_values,r_values10,rho_values,P_values,T_values,L_values,F_con_values,F_rad_values=simulation.run_simulation(simulation,p)
    star_zone10=cross_section_mod(r_values10,L_values,F_con_values, ax=axs[1])
    axs[1].text(0,0.9,f"Core size: {simulation.core_size:.3f}")
    axs[1].text(0,0.5,f"Outer Conv zone: {simulation.conv_zone:.3f}")
    plt.savefig("report/Figures/simulation_part2.pdf")
    plt.show()
    
    #aditional plot showing off the active zones as functions of the radius for all the different simulations ran
    height=np.linspace(0,3,101)
    core=[] ; surface=[]
    for i in range(len(height)):
        core.append(0.1)
        surface.append(1)
    core=np.array(core) ; surface=np.array(surface)
    plt.plot(r_values1/r_values1[0],star_zone1,linewidth=0.8, label=rf"$\rho_0/10$")
    plt.plot(r_values2/r_values2[0],star_zone2,linewidth=0.8, label=rf"$\rho_0\times10$")
    plt.plot(r_values3/r_values3[0],star_zone3,linewidth=0.8, label=rf"$T_0/2$")
    plt.plot(r_values4/r_values4[0],star_zone4,linewidth=0.8, label=rf"$T_0\times2$")
    plt.plot(r_values5/r_values5[0],star_zone5,linewidth=0.8, label=rf"$L_0/2$")
    plt.plot(r_values6/r_values6[0],star_zone6,linewidth=0.8, label=rf"$L_0$")
    plt.plot(r_values7/r_values7[0],star_zone7,linewidth=0.8, label=rf"$R_0$")
    plt.plot(r_values8/r_values8[0],star_zone8,linewidth=0.8, label=rf"$R_0\times2$")
    plt.plot(r_values9/r_values9[0],star_zone9,linewidth=0.8, label=rf"$P_0/10$")
    plt.plot(r_values10/r_values10[0],star_zone10,linewidth=0.8, label=rf"$P_0\times10$")
    plt.plot(core,height,color="k", linestyle= "--", linewidth = 2, alpha=0.3, label="Core")
    plt.plot(surface,height,color="k", linestyle= "--", linewidth = 2, alpha=0.3, label="Surface")
    plt.legend(loc="lower center",fontsize=8)
    plt.xlabel(r"$R/R_\odot$")
    plt.ylabel(r"Zones")
    plt.savefig("report/Figures/construction.pdf")
    plt.show()


p=0.1    
#function simulationing the behaviour of the different parameters 
simulation_parameters(p)     

#Sanity checks for the surface condition described by the parameters in the task
simulation=EnergyTransportation(1,1,1,1.42*10**(-7),5770, print_sanity=True)  
print("Opacity Sanity Check for Task Parameters")
simulation.opacity_sanity(tol=0.005)  
Q_PP1_values, Q_PP2_values,Q_PP3_values,Q_CNO_values,Q_total,nabla_star, nabla_stable, nabla_ad,m_values,r_values,rho_values,P_values,T_values,L_values,F_con_values,F_rad_values=simulation.run_simulation(simulation,p)
core_sz=simulation.core_size
conv_sz=simulation.conv_zone
ss=cross_section(r_values,L_values,F_con_values,cr=core_sz,cv=conv_sz, sanity=False, savefig=True, name="_sanity_check")
gradient_plot(r_values,nabla_star, nabla_stable,nabla_ad,p,name="_sanity_check")

#redefinition of a class returning the spesific opacity value for the testing we wish to perform
class RedefOpacity(EnergyTransportation):
    def opacity(self, rho, T):
        return 3.98 
    
#test of the variables from example 5.1 in the notes from AST3310, presented in the table 3 in report
simulation=RedefOpacity(L0=1,R0=0.84,M0=0.99, rho0=55.9/(1.408*10**3 ),T0= 0.9*10**6, print_sanity=True)  
simulation.print_expected()

#data analysis for fitted model, making use of the best parameters
simulation=EnergyTransportation(L0=1.7,R0=1,M0=1.1, rho0=70*1.42*10**(-7),T0= 2*5770, print_sanity=False)       
Q_PP1_values, Q_PP2_values,Q_PP3_values,Q_CNO_values,Q_total,nabla_star, nabla_stable, nabla_ad,m_values,r_values,rho_values,P_values,T_values,L_values,F_con_values,F_rad_values=simulation.run_simulation(simulation,p)
core_sz=simulation.core_size ; conv_sz=simulation.conv_zone
star_zone=cross_section(r_values,L_values,F_con_values,cr=core_sz, cv=conv_sz,sanity=False, savefig=True, name="final")
gradient_plot(r_values,nabla_star, nabla_stable,nabla_ad,p,name="final")
energies_plot(r_values,Q_PP1_values, Q_PP2_values,Q_PP3_values,Q_CNO_values,Q_total,p,name="final")
flux_plot(r_values,F_con_values, F_rad_values,p,name="final")
variables_special_scaling(m_values,r_values,rho_values,P_values,L_values,T_values,p,name="final")

