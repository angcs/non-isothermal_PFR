import pandas as pd
import numpy as np

## read kinetics data
data = pd.read_excel("ibbett.xlsx", sheet_name="summary")
data.set_index("reaction", inplace=True)
data = data.astype("float64")
data

## calculate kinetic constants
def kinetic_constant(dataframe, pH, temperature):
    R = 8.314e-3  # [kJ/mol K]

    for rxn in dataframe.index:
        A = dataframe.loc[rxn,"A"]
#         AH = dataframe.loc[rxn,"AH+"]
        AH = 0
        Ea = dataframe.loc[rxn,"Ea"]
        
        dataframe.loc[rxn,"k"] = (A + AH*10**(-pH))*np.exp(-Ea/(R*T))

pH = 1
T = 433
kinetic_constant(data,pH,T)

data.loc[:,["A","Ea","k"]]  ## working w/o acid

from scipy.integrate import odeint
import matplotlib.pyplot as plt

## Create kinetics constant value for different conditions
def model(initial_conditions, volume, volumetric_flowrate):

    ## Create ODE solver and plotting function
    ks = data.loc["ks","k"]
    k2 = data.loc["k2","k"]
    k3 = data.loc["k3","k"]
    k4 = data.loc["k4","k"]
    kfa = data.loc["kfa","k"]

    def odes(var, V):
        # assign each ODE to a vector element
        # A = hc
        # B = xlog
        # C = xls
        # D = fur
        
        A = var[0]
        B = var[1]
        C = var[2]
        D = var[3]

        # define each ODE
        dAdV = -ks*A
        dBdV = ks*A + (-k2-kfa)*B
        dCdV = k2*B + (-k3)*C
        dDdV = k3*C - k4*D

        return [dAdV, dBdV, dCdV, dDdV]

    # initial conditions [kg/m3]
    if type(initial_conditions) == int or type(initial_conditions) == float:
        var0 = [initial_conditions, 0, 0, 0]
        
    elif type(initial_conditions) == list:
        var0 = initial_conditions
        
    # declare a volume vector (reactor volume window to integrate)
    V = np.linspace(0,volume,volume*500)  ## 500 volume step
    tau = V*60/volumetric_flowrate
    
    ## defining "var" as ODE solver results
    var = odeint(odes,var0,V)
    
    ## retreiving each species profile from ODE solver
    A = var[:,0]
    B = var[:,1]
    C = var[:,2]
    D = var[:,3]

    # plot the results
    plt.style.use("seaborn")
    plt.plot(V,A,"r-", label = "hc")
    plt.plot(V,B,"b-", label = "xlog")
    plt.plot(V,C,"y:", label = "xls")
    plt.plot(V,D,"g:", label = "fur")

    # plt.axis([0,10,0,16])
    plt.legend(loc=9,ncol=3)
    plt.ylabel('varentration [kg/m3]')
    plt.xlabel('volume [m3]')
    plt.title(str(T) +" Celsius at pH " + str(pH))
    plt.show()
    
    # plot the results
    plt.style.use("seaborn")
    plt.plot(tau,A,"r-", label = "hc")
    plt.plot(tau,B,"b-", label = "xlog")
    plt.plot(tau,C,"y:", label = "xls")
    plt.plot(tau,D,"g:", label = "fur")

    # plt.axis([0,10,0,16])
    plt.legend(loc=9,ncol=3)
    plt.ylabel('varentration [kg/m3]')
    plt.xlabel('residence time [min]')
    plt.title(str(T) +" Celsius at pH " + str(pH))
    plt.show()

V = model(86.4, 100, 98.9)
