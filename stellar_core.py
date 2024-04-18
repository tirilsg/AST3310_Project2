#program defining a class that models energy production within a stellar core, grounded in theory presented in the attatched report
#implements a number of relevant constants, as well as sizes like the mass of relevant isotopes created within the core
#defines a number of functions for number density of elements, reaction rates and energy output as a result of fusion
#a function testing the accuracy of the model is also implemented, sanity_check(), which can be called if needed

from math import*
import numpy as np
import matplotlib.pyplot as plt
import sys

class StellarCore:
    def __init__(self,rho,T):
        #definitions of relevant constants
        self.u_VJ=931.4943
        self.Na=6.022141*10**(23)
        self.g_si=1/(10**6*self.Na)
        self.MeV_J=10**6*1.602176*10**(-19)
        self.u=1.660539*10**(-27)
        self.rho=rho 
        self.T=T
        self.T9=T*10**(-9)

        #definition of relevant mass fractions
        self.X=0.7 ; self.Y3He=10**(-10) ; self.Y=0.29
        self.Z7Li=10**(-7) ; self.Z7Be=10**(-7) ; self.Z14N=10**(-11)

        #masses of isotopes relevant for the cycle reactions, units u
        self.mH1=1.007825 ; self.mD2=2.014102 ; self.mHe3=3.016029 ; self.mHe4=4.002603 
        self.mBe7=7.016929 ; self.mBe8=8.005305 ; self.mLi7=7.016004 ; self.mB8=8.024607 
        self.mC12=12.000000 ; self.mC13=13.003355 ; self.mN13=13.005739 ; self.mN14=14.003074
        self.mN15=15.000109 ; self.mO15=15.003066

        #energy lost to neutrinos per cycle, unit MeV
        self.E_PP0=0.265 ; self.E_PP2=0.815 ; self.E_PP3=6.711
        self.E_CNO1=0.707 ; self.E_CNO2=0.997

        #calculations of the energy output for each reaction cycle
        self.calculate_densities()
        self.proportionality()
        self.reaction_rates(True)
        self.estimated_energy()
        self.cycle_relation()

    def einsteins_principle(self,m_0,m_1): #Einstein's mass-energy equivalency principle, Equation 1 in the report
        return (m_0-m_1)*self.u_VJ

    def number_density(self, mass_frac, Z): #number density of an arbitrary element Z, Equation 2 in the report
        return mass_frac*self.rho/(Z*self.u)
    
    def calculate_densities(self): #calculation of the number density for all relevant isotopes
        #calls number_density for each of the relevant isotopes existing within the stellar core
        self.n1H=self.number_density(self.X,1)
        self.n3He=self.number_density(self.Y3He,3)
        self.n4He=self.number_density(self.Y,4)
        self.n7Li=self.number_density(self.Z7Li,7)
        self.n7Be=self.number_density(self.Z7Be,7)
        self.n14N=self.number_density(self.Z14N,14)
        self.ne=(self.n1H+2*self.n3He+2*self.n4He)+(3*self.n7Li+3*self.n7Be+7*self.n14N)
    
    def proportionality(self): #definition of all relevant proportionality functions, presented in the Table III in the report
        #definition of constants needed in the estimation of proportionality functions
        T9=self.T9 
        T9S=T9/(1+4.95*10**(-2)*T9)
        T9SS=T9/(1+0.759*T9)
        self.NA_lambd_lpp=self.g_si*4.01*10**(-15)*T9**(-2/3)*np.exp(-3.380*T9**(-1/3))*(1+0.123*T9**(1/3)+1.09*T9**(2/3)+0.938*T9)
        self.NA_lambd_l33=self.g_si*6.04*10**10*T9**(-2/3)*np.exp(-12.276*T9**(-1/3))*(1+0.034*T9**(1/3)-0.522*T9**(2/3)-0.124*T9+0.353*T9**(4/3)+0.213*T9**(5/3))
        self.NA_lambd_l34=self.g_si*5.61*10**6*T9S**(5/6)*T9**(-3/2)*np.exp(-12.826*T9S**(-1/3))
        lambda_e7=(1.34*10**(-10)*T9**(-1/2)*(1-0.537*T9**(1/3)+3.86*T9**(2/3)+0.0027*T9**(-1)*np.exp(2.515*10**(-3)*T9**(-1))))*self.g_si
        if self.T<10**6: #test that insures the existence of the upper limit for Li_3^7 is considered
            if lambda_e7 > (1.57*10**(-7))*10**-6/(self.ne*self.Na):
                lambda_e7 = (1.57*10**(-7))*10**-6/(self.ne*self.Na)
        self.NA_lambd_le7=lambda_e7
        self.NA_lambd_l17l=self.g_si*(1.096*10**9*T9**(-2/3)*np.exp(-8.472*T9**(-1/3))-4.830*10**8*T9SS**(5/6)*T9**(-3/2)*np.exp(-8.472*T9SS**(-1/3))+1.06*10**10*T9**(-2/3)*np.exp(-30.442*T9**(-1)))
        self.NA_lambd_l17=self.g_si*(3.11*10**5*T9**(-2/3)*np.exp(-10.262*T9**(-1/3))+2.53*10**3*T9**(-3/2)*np.exp(-7.306*T9**(-1)))
        self.NA_lambd_lp14=self.g_si*((4.90*10**7*T9**(-2/3)*np.exp(-15.228*T9**(-1/3)-0.092*T9**2)*(1+0.027*T9**(1/3)-0.778*T9**(2/3)-0.149*T9+0.261*T9**(4/3)+0.127*T9**(5/3)))+2.37*10**3*T9**(-3/2)*np.exp(-3.011*T9**(-1))+2.19*10**4*np.exp(-12.53*T9**(-1)))

    def reaction_rates(self, alter): #definition of rate of reactions, presented and explained as Equation 3 in the report
        self.r_lpp= self.n1H*self.n1H/(self.rho*(1+1))*self.NA_lambd_lpp
        self.r_l33 = self.n3He*self.n3He/(self.rho*(1+1))*self.NA_lambd_l33
        self.r_l34 = self.n3He*self.n4He/(self.rho)*self.NA_lambd_l34
        self.r_le7 = self.ne*self.n7Be/(self.rho)*self.NA_lambd_le7
        self.r_l17l = self.n7Li*self.n1H/(self.rho)*self.NA_lambd_l17l
        self.r_l17 = self.n7Be*self.n1H/(self.rho)*self.NA_lambd_l17
        self.r_lp14 = self.n14N*self.n1H/(self.rho)*self.NA_lambd_lp14
        #if the argument alter is True, reaction rates are altered to make sure energy consumption in the stellar core
        #does not exceed the rate of production of elements, in accordance with expressions presented in Table IV in the report
        if alter==True: 
            const=self.r_lpp/(2*self.r_l33+self.r_l34)
            #if the rate of production of 2*l33 (lpp needs to happen 2 times in order for l33 to happen) and l34 is higher than lpp, rates needs to be adjusted
            if (2*self.r_l33+self.r_l34) > self.r_lpp:
                self.r_l33 = const*self.r_l33
                self.r_l34 = const*self.r_l34
            const=self.r_l34/(self.r_le7+self.r_l17)
            #the fallout of l34 can be either reaction l17 and le7, and if the sum of these consumption rates exceeds the rate of production l34;
            if self.r_le7+self.r_l17 > self.r_l34:
                self.r_le7 = const*self.r_le7
                self.r_l17 = const*self.r_l17
            #if more of l17l happens more frequently than le7, the rates need to be adjusted to be realistic
            if self.r_l17l>self.r_le7 :
                self.r_l17l = self.r_le7

    def estimated_energy(self): #definition of total energy output of each reaction of fusion
        #calculates energy output from einssteins principle for all reactions of fusion occurring in the stellar core
        self.E_l1 = self.MeV_J*(self.einsteins_principle(2*self.mH1 , self.mD2) - self.E_PP0) #neutrino results in energy loss self.E_PP0
        self.E_l1d = self.MeV_J*self.einsteins_principle(self.mD2 + self.mH1 , self.mHe3)
        self.E_l33 = self.MeV_J*self.einsteins_principle(2*self.mHe3 , self.mHe4 + 2*self.mH1)
        self.E_l34 = self.MeV_J*self.einsteins_principle(self.mHe3 + self.mHe4 , self.mBe7)
        self.E_le7 = self.MeV_J*(self.einsteins_principle(self.mBe7 , self.mLi7) - self.E_PP2)
        self.E_l17l = self.MeV_J*self.einsteins_principle(self.mLi7 + self.mH1 , 2*self.mHe4)
        self.E_l17 = self.MeV_J*self.einsteins_principle(self.mBe7 + self.mH1 , self.mB8)
        self.E_dc = self.MeV_J*(self.einsteins_principle(self.mB8 , 2*self.mHe4) - self.E_PP3)
        self.E_lp12 = self.MeV_J*self.einsteins_principle(self.mC12 + self.mH1 , self.mN13)
        self.E_l13 = self.MeV_J*(self.einsteins_principle(self.mN13 , self.mC13) - self.E_CNO1)
        self.E_lp13 = self.MeV_J*self.einsteins_principle(self.mC13 + self.mH1 , self.mN14)
        self.E_lp14 = self.MeV_J*self.einsteins_principle(self.mN14 + self.mH1 , self.mO15)
        self.E_l15 = self.MeV_J*(self.einsteins_principle(self.mO15 , self.mN15) - self.E_CNO2)
        self.E_lp15 = self.MeV_J*self.einsteins_principle(self.mN15 + self.mH1 , self.mC12 + self.mHe4)
    
    def cycle_relation(self): #relation between each fusion-cycle and its energy output calculated by Equation 6 in the report
        self.Q_PP1 = self.r_l33 * (self.E_l33 + 2*(self.E_l1 + self.E_l1d)) 
        self.Q_PP2 = self.r_l34 * (self.E_l34 + (self.E_l1 +self.E_l1d)) + self.r_le7*(self.E_le7 + self.E_l17l)
        self.Q_PP3 = self.r_l17 * (self.E_l17+self.E_dc)+self.r_l34*(self.E_l34+(self.E_l1 + self.E_l1d))
        self.Q_CNO = self.r_lp14 * (self.E_lp12+self.E_l13+self.E_lp13+self.E_lp14+self.E_lp15+self.E_l15)
    
    def sanity_check(self,t_v,err): #a function that will check the accuracy of the developed model
        print(f"Sanity check for temperature T={self.T} and density rho={self.rho}")
        print("Expected value:   Estimated value:   Relative Error:")
        e_v=np.array([
            self.r_lpp*(self.E_l1 + self.E_l1d)*self.rho,
            self.r_l33*self.E_l33*self.rho,
            self.r_l34*self.E_l34*self.rho,
            self.r_le7*self.E_le7*self.rho,
            self.r_l17l*self.E_l17l*self.rho,
            self.r_l17*(self.E_l17+self.E_dc)*self.rho,
            self.r_lp14*(self.E_lp12+self.E_l13+self.E_lp13+self.E_lp14+self.E_lp15+self.E_l15)*self.rho])
        #prints the expected and estimated value for the energy production per unit mass, as well as the relative error is printed to the terminal
        list=[]
        for i in range(len(t_v)):
            rel=abs((e_v[i]-t_v[i])/t_v[i])
            if rel>err:
                list.append(i)
            print(f"{t_v[i]:.3}          {e_v[i]:.3}              {rel:.3}")
        if len(list)!=0:
            for i in range(len(list)):
                print(f"The sanity check is failed for Sanity Check Equation {list[i]+1}")
            #sys.exit() #exits entire program, meaning plots and other code trying to be ran after calling a sanity check that fails will not be able to run
