"""
Spyder Editor
This is a temporary script file.
"""
# -*- coding: utf-8 -*-
"""
Analysis of Pang and Stoylar, conducts experiments using proposed Scheme B
with initial parameters suggested in paper.  
@authors Fasil Cheema, Ajanthan Mathialagan, Marisa Signorile
"""

def schemeB_fluid_with_A(X_0,Y_0,time,t_dependence = False):
    import numpy as np
    
    #Setup initial parameters
    beta    = 1
    epsilon = 0.2
    gamma   = 2
    dt      = 0.0001
    
    #compute number of time steps
    num_steps = int(time/dt)
    
    #initialize scheme variables
    t    = np.arange(0,time,dt)
    X    = np.zeros(num_steps)
    Y    = np.zeros(num_steps)
    
    #Determines Alpha depending on which case we are using (time dependent or constant)
    if t_dependence == False:
        alpha = np.zeros(num_steps) + 1000
    else:    
        alpha   = 1000 + 200*np.sin((2*np.pi*t/120))
    
    X[0] = X_0
    Y[0] = Y_0
    """
    Based on the following equations compute update conditions:
    X' = - (gamma)*Y' - epsilon*Y
    Y' = beta * X - alpha 
    """
    
    for i in range(1,num_steps):
        dYdt = (beta*X[i-1])-alpha[i-1]
        dXdt = (-(gamma)*dYdt)-(epsilon*Y[i-1])
        
    
        X[i] = X[i-1] + (dXdt*dt)
        Y[i] = Y[i-1] + (dYdt*dt)
        
    return t,X,Y

def schemeB_simulation_with_A(X_0,Y_0,time,t_dependence = False):
    import random
    import numpy as np
    
    #Setup initial parameters
    rnd_seed = 343
    beta     = 1
    epsilon  = 0.2
    gamma    = 2
    dt       = 0.0001
    random.seed(rnd_seed)
    
    #compute number of time steps
    num_steps = int(time/dt)
    
    #initialize scheme variables
    t    = np.arange(0,time,dt)
    X    = np.zeros(num_steps)
    Y    = np.zeros(num_steps)
    
    #Determines Alpha depending on which case we are using (time dependent or constant)
    if t_dependence == False:
        alpha = np.zeros(num_steps) + 1000
    else:    
        alpha   = 1000 + 200*np.sin((2*np.pi*t/120))
    
    X[0] = X_0
    Y[0] = Y_0
    
    
    """
    Models a CTMC with probabilities updating and at every instance an events probability is calculated
    """
    
    
    for i in range(1,num_steps):
        #Compute and update probabilities to be used
        p_arrival  = alpha[i-1]*dt
        p_accept   = beta*X[i-1]*dt
        p_addition = epsilon*abs(Y[i-1])*dt
        
        #Create a random number to represent if event occurs
        curr_arriv = np.random.rand()
        curr_acc   = np.random.rand()
        curr_add   = np.random.rand()
        
        dX = 0
        dY = 0
        
        
        '''
        Following code block checks to see if an event occurs and if it does approriately updates 
        dX and dY for the current step
        '''
        if curr_arriv <= p_arrival:
            dY = dY - 1
            dX = dX + gamma
            
        if curr_acc <= p_accept:
            dY = dY + 1
            if (X[i-1]-gamma) > 0:
                dX = dX - gamma
        
        if curr_add <= p_addition:
            if (X[i-1]) >= 1:
                dX = dX - np.sign(Y[i-1])
            elif (X[i-1]) == 0:
                if ((Y[i-1]) < 0):
                    dX = dX + 1
        
        #Updates X, and Y
        X[i] = X[i-1] + dX
        Y[i] = Y[i-1] + dY

    return t,X,Y

def plot_graphs(X0,Y0,t,X_sim,Y_sim,X_fluid,Y_fluid, t_dependence = False):
    #Plots the figures, and formatting figures
    import numpy as np
    import matplotlib.pyplot as plt
    
    #set up x axis ticks
    if t_dependence == False:
        tau = 5
    else:
        tau = 50
    major_ticks = np.arange(0,t[len(t)-1]+1,tau)
    
    #Format titles,ticks, colors, and legends of plot then plot
    fig_title = "Y(t) v Time with initial parameters X(0) = {val1}, Y(0) = {val2}".format(val1 = X0,val2 = Y0 )
    plt.plot(t,Y_sim,'-r', label="sim")
    plt.plot(t,Y_fluid,'-b', label="fluid")
    plt.legend(loc="lower right")
    plt.xlabel("Time (s)")
    plt.xticks(major_ticks)
    plt.ylabel("Y(t)")
    plt.title(fig_title)
    plt.show()
    
    if t_dependence == False:
        alpha = 1000
        beta = 1 
        
        X_sim = X_sim - (alpha/beta)
        X_fluid = X_fluid - (alpha/beta)
    
    fig_title = "X(t) v Time with initial parameters X(0) = {val1}, Y(0) = {val2}".format(val1 = X0,val2 = Y0 )
    plt.plot(t,(X_sim),'-r', label="sim")
    plt.plot(t,(X_fluid),'-b', label="fluid")
    plt.legend(loc="lower right")
    plt.xlabel("Time (s)")
    plt.xticks(major_ticks)
    
    if t_dependence == False:
        plt.ylabel(r'X(t) - $ \frac{\lambda}{\beta}$')
    else:
        plt.ylabel('X(t)')
    plt.title(fig_title)
    plt.show()
    

def run_experiment(param_set):
    for j in range(int(len(param_set)/4)):
        X0 = param_set[4*j]
        Y0 = param_set[4*j + 1]
        tf = param_set[4*j + 2]
        t_dependence = param_set[4*j + 3]
        
        t,X1,Y1 = schemeB_fluid_with_A(X0,Y0,tf,t_dependence)
        t,X2,Y2 = schemeB_simulation_with_A(X0,Y0,tf,t_dependence)
        
        plot_graphs(X0,Y0,t,X2,Y2,X1,Y1,t_dependence)


#Set of parameters to be used in format (X0,Y0,time duration,time dependence,....)    
param_set = [0,0,50,False,0,1000,50,False,2000,0,50,False,2000,-1000,50,False,0,0,500,True,2000,-1000,500,True]
run_experiment(param_set)
