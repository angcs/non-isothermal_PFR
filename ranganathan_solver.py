import pandas as pd
import numpy as np

## read kinetics data
# A [1/s]
# m [-]
# Ea [kj/mol]
data = pd.read_excel("ranganathan.xlsx", sheet_name="unit")
data.set_index("reaction", inplace=True)
data = data.astype("float64")
data

from scipy.integrate import odeint
import matplotlib.pyplot as plt

R = 8.314e-3  # [kJ/mol K]

def model(hc0, ch, T0, Ta, df, V, vf, U, a, hstrm):
    T0 = T0 + 273  # reactor temp
    Ta = Ta + 273  # jacket temp
    
    def odes(var, V):
#         print("should be last B", b)
        ## assign each ODE to a vector element
        # A = hc
        # B = xlos
        # C = fur
        # T = temp of reactor
        A = var[0]
        B = var[1]
        C = var[2]
        T = var[3]
        
        ## define rate constant for ODE depending on ch, T
        for r in df.index:
            pexpo = df.loc[r,"A"]
            m = df.loc[r,"m"]
            Ea = df.loc[r,"Ea"]
            
            df.loc[r,"k"] = (pexpo*ch**m)*np.exp(-Ea/(R*T))
            
        k1 = df.loc["k1","k"]
        k2 = df.loc["k2","k"]

        # define series of ODE
        dAdV = (-k1*A)
        dBdV = (k1*A + (-k2)*B)
        dCdV = (k2*B)
        dTdV = ((dBdV*28.2+dCdV*(-149))+(U*a*(Ta-T)*3.6))/hstrm
        
        return [dAdV, dBdV, dCdV, dTdV]
      
    ## declare a volume vector (reactor volume window to integrate)
    vol = np.linspace(0,V,V*500)  # 500 volume step
    tau = vol*60/vf  # equivalent residence time step
    
    ## initial conditions
    var0 = [hc0,0,0,T0]
    
    ## defining "var" as ODE solver results
    var = odeint(odes,var0,vol)
    
    ## retreiving each species profile from ODE solver
    A = var[:,0]
    B = var[:,1]
    C = var[:,2]
    T = var[:,3]-273  #temp in celsius
    
    ## plot the results
    # initiate plot
    plt.style.use("default")
    fig, ax1 = plt.subplots()
    plt.title("Inlet of " + str(T0-273) +" °C at " + str(ch) + " wt% acid")
    
    # axis 1 (reactants)
    ax1.plot(vol,A,"tab:olive", label = "hc")
    ax1.plot(vol,B,"tab:orange", label = "xls")
    ax1.plot(vol,C,"tab:blue", label = "fur")
    ax1.set_xlabel('Volume [m3]')
    ax1.set_ylabel('Concentration [kg/m3]')
    
    # axis 2 (temp)
    ax2 = ax1.twinx()
    ax2.plot(vol,T,"k", label = "T")
    ax2.set_ylabel('Temperature [°C]')
    ax2.set_ylim(0,200)
    
    # axis 3 (corresponding residence time)
    ax3 = ax1.twiny()
    ax3.plot(tau,np.zeros(var.shape[0]), alpha=0)
    ax3.set_xlabel('Time [min]')
    
    # compile legends
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc="center right")
    
    # show plot
    plt.show()
    
    return var

# hc0 = initial hemicellulose concentration [kg/m3]
# ch = acid concentration [wt%]
# T0 = inlet temperature [°C]
# Ta = steam jacket temperautre [°C]
# U = overall heat exchange coefficient [W/m2 K]
# a = specific heat exchange area [m2/m3]
# hstrm = stream enthalpy [kJ/hr K]
# df = kinetics data [dataframe]
# V = reactor volume [m3]
# vf = volumetric flowrate [m3/hr]
    
solved = model(hc0=19, 
               ch=0.02, 
               T0=87,
               Ta=160,
               U=4000,
               a=10,
               hstrm=478148,
               df=data,
               V=80, 
               vf=317)

import pandas as pd
import numpy as np

## read kinetics data
# A [1/s]
# m [-]
# Ea [kj/mol]
data = pd.read_excel("ranganathan.xlsx", sheet_name="unit")
data.set_index("reaction", inplace=True)
data = data.astype("float64")
data

from scipy.integrate import odeint
import matplotlib.pyplot as plt

R = 8.314e-3  # [kJ/mol K]

def model(hc0, ch, T0, Ta, df, V, vf, U, a, hstrm):
    T0 = T0 + 273  # reactor temp
    Ta = Ta + 273  # jacket temp
    
    def odes(var, V):
#         print("should be last B", b)
        ## assign each ODE to a vector element
        # A = hc
        # B = xlos
        # C = fur
        # T = temp of reactor
        A = var[0]
        B = var[1]
        C = var[2]
        T = var[3]
        
        ## define rate constant for ODE depending on ch, T
        for r in df.index:
            pexpo = df.loc[r,"A"]
            m = df.loc[r,"m"]
            Ea = df.loc[r,"Ea"]
            
            df.loc[r,"k"] = (pexpo*ch**m)*np.exp(-Ea/(R*T))
            
        k1 = df.loc["k1","k"]
        k2 = df.loc["k2","k"]

        # define series of ODE
        dAdV = (-k1*A)
        dBdV = (k1*A + (-k2)*B)
        dCdV = (k2*B)
        dTdV = ((dBdV*28.2+dCdV*(-149))+(U*a*(Ta-T)*3.6))/hstrm
        
        return [dAdV, dBdV, dCdV, dTdV]
      
    ## declare a volume vector (reactor volume window to integrate)
    vol = np.linspace(0,V,V*500)  # 500 volume step
    tau = vol*60/vf  # equivalent residence time step
    
    ## initial conditions
    var0 = [hc0,0,0,T0]
    
    ## defining "var" as ODE solver results
    var = odeint(odes,var0,vol)
    
    ## retreiving each species profile from ODE solver
    A = var[:,0]
    B = var[:,1]
    C = var[:,2]
    T = var[:,3]-273  #temp in celsius
    
    ## plot the results
    # initiate plot
    plt.style.use("default")
    fig, ax1 = plt.subplots()
    plt.title("Inlet of " + str(T0-273) +" °C at " + str(ch) + " wt% acid")
    
    # axis 1 (reactants)
    ax1.plot(vol,A,"tab:olive", label = "hc")
    ax1.plot(vol,B,"tab:orange", label = "xls")
    ax1.plot(vol,C,"tab:blue", label = "fur")
    ax1.set_xlabel('Volume [m3]')
    ax1.set_ylabel('Concentration [kg/m3]')
    
    # axis 2 (temp)
    ax2 = ax1.twinx()
    ax2.plot(vol,T,"k", label = "T")
    ax2.set_ylabel('Temperature [°C]')
    ax2.set_ylim(0,200)
    
    # axis 3 (corresponding residence time)
    ax3 = ax1.twiny()
    ax3.plot(tau,np.zeros(var.shape[0]), alpha=0)
    ax3.set_xlabel('Time [min]')
    
    # compile legends
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc="center right")
    
    # show plot
    plt.show()
    
    return var

# hc0 = initial hemicellulose concentration [kg/m3]
# ch = acid concentration [wt%]
# T0 = inlet temperature [°C]
# Ta = steam jacket temperautre [°C]
# U = overall heat exchange coefficient [W/m2 K]
# a = specific heat exchange area [m2/m3]
# hstrm = stream enthalpy [kJ/hr K]
# df = kinetics data [dataframe]
# V = reactor volume [m3]
# vf = volumetric flowrate [m3/hr]
    
solved = model(hc0=19, 
               ch=0.02, 
               T0=87,
               Ta=160,
               U=4000,
               a=10,
               hstrm=478148,
               df=data,
               V=80, 
               vf=317)
