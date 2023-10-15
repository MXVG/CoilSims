from femm import *
from math import *
import numpy as np
import matplotlib.pyplot as plt
import os

import pandas as pd
from scipy import integrate

mydir="./"

# Projectile sizing:
D_p = 9.525 # Projectile diameter
L_p = 22.225 # Projectile length

M_p = 12.416 # Projectile mass (g)

# Starting velocity
startingVelocity = 0

# System properties
VOLTAGE = 400 # Volts.
N_STAGES = 5 # (unitless)


MAX_TIME = 0.02 #20ms
TIME_STEP = 0.00005 #0.05ms


# Coil dimensions:
W_mm = 0 # Wire diameter
D_ic = 0 # Coil inner diameter
D_oc = 0 # Coil outer diameter
L_c = 0 # Coil length
N_Coil = 0 # Number of turns
R_Coil = 0 # Resistance of coil


# Problem setup:
#Z_proj = -5 # Offset between front face of projectile and reverse face of copper. 
#Z_step = 0 #0.5 # Resolution of the work calculation. 
E_force = 0.1 # If F < E_force, stop iterating.

df = pd.read_excel('StageSetup.xlsx')

ts = np.linspace(0, MAX_TIME, int(MAX_TIME/TIME_STEP))


def I(I, t, k1, k2):
    dIdt = [I[1], -k1*I[1]-k2*I[0]]
    return dIdt

def findClosest(arr, target):
    closest_idx = 0
    min_diff = abs(arr[0] - target)
    for i in range(1, len(arr)):
        diff = abs(arr[i] - target)
        if diff < min_diff:
            closest_idx = i
            min_diff = diff
    return closest_idx




def runSim(simName, experimentName, csv, value, startingVelocity):
    
    
    Z_step = 1
    Z_proj = -1 * float(df['Z_proj'][value])
    
    # Coil dimensions:
    W_mm = float(df['W_mm'][value]) # Wire diameter
    D_ic = float(df['D_ic'][value]) # Coil inner diameter
    D_oc = float(df['D_oc'][value]) # Coil outer diameter
    L_c = float(df['L_c'][value]) # Coil length
    N_Coil = float(df['N_Turns'][value]) # Number of turns
    R_Coil = float(df['R'][value]) # Resistance of coil
    
    L = float(df['L'][value]) # Inductance of coil
    C = float(df['C'][value]) # Capacitance
    V = float(df['V'][value]) # Voltage
    
    current = integrate.odeint(I, [0, V/L], ts, args=(R_Coil/L, 1/(L*C)) )
    
    # Plot RLC Current Response
    plt.cla()
    plt.plot(ts, current[:, 0])
    plt.title("Current with L={0} C={1} R={2} Vi={3}".format(L, C, R_Coil, V))
    plt.ylabel('I [A]')
    plt.xlabel('t [s]')
    plt.grid()
    plt.savefig("{0}/Response Stage{1}.png".format(experimentName, value))
    plt.cla()
    
    
    currI = current[:, 0][0]
    simTime = [TIME_STEP]
    
    acc = []
    acc += [0]
    vel = [startingVelocity]
    pos = [0]
    
    print ("-----------------------")
    print("Start of simulation: STAGE {} of {}".format(value + 1, N_STAGES))
    print ("-----------------------")
    # CONSTANTS ##############
    PGROUP = 1
    CMAT = df['CMAT'][value]
    PMAT = "1010 Steel"
    AMAT = "Air"
    
    print("Coil turns:\t\t{:>10.3f}".format(N_Coil))
    print("Coil resistance:\t{:>10.3f} mOhm".format(R_Coil*1000))
    print("Capacitor Capacitance:\t\t{:>10.3f} uF".format(C / 1E-6))
    print("Starting velocity:\t\t{:>10.3f} m/s".format(startingVelocity))
    
    Z_sweep = L_c - Z_proj 
    def drawRect(x1,y1, x2,y2):
        mi_addnode(x1,y1)
        mi_addnode(x2,y1)
        mi_addnode(x2,y2)
        mi_addnode(x1,y2)
        mi_addsegment(x1,y1, x2,y1)
        mi_addsegment(x2,y1, x2,y2)
        mi_addsegment(x2,y2, x1,y2)
        mi_addsegment(x1,y2, x1,y1)
        
        
    # BUILD THE SIMULATION ##############
    # Make a new magnetostatics document
    newdocument(0)
    # Set up the problem and save our new scratch fille as something:
    mi_probdef(0,"millimeters","axi",1E-8)
    mi_saveas(mydir+"TEMP.fem")
    mi_getmaterial(CMAT)
    mi_getmaterial(PMAT)
    mi_getmaterial(AMAT)
    # Setup coil circuit
    # "COIL", 100A, (0 -> parallel, 1 ->series)
    mi_addcircprop("COIL", currI, 1)
    ## Build the projectile, and assign it to PGROUP
    mi_addnode(D_p/2, Z_proj)
    mi_addnode(D_p/2, Z_proj-L_p)
    mi_addnode(0, Z_proj-L_p)
    mi_addnode(0, Z_proj)
    mi_addsegment(D_p/2, Z_proj, 0, Z_proj)
    mi_addsegment(0, Z_proj, 0, Z_proj-L_p)
    mi_addsegment(0, Z_proj-L_p, D_p/2, Z_proj-L_p)
    mi_addsegment(D_p/2, Z_proj-L_p, D_p/2, Z_proj)
    mi_addblocklabel(D_p/4, Z_proj-L_p/2)
    mi_selectlabel(D_p/4, Z_proj-L_p/2)
    mi_setblockprop(PMAT,1,0,"",0,0,0)
    # Since everything in existence is G0, we can select everything to turn it into projectile
    mi_selectgroup(0)
    mi_setgroup(PGROUP)
    
    
    ## Build the coil
    drawRect(D_ic/2, 0, D_oc/2, L_c)
    mi_addblocklabel((D_ic+D_oc)/4, L_c/2)
    mi_selectlabel((D_ic+D_oc)/4, L_c/2)
    mi_setblockprop(CMAT,1,0,"COIL",1,0,N_Coil)
    

            
    # Create boundary (bc = 0 for Dirichlet, 1 for Neumann)
    # n,R,x,y,bc
    mi_makeABC(7,200,0,0,0)
    
    # Make the air... air.
    mi_clearselected()
    mi_addblocklabel(50, -10)
    mi_selectlabel(50, -10)
    mi_setblockprop(AMAT,1,0,"",0,0,0)
    

    forceMeasured = 1
    distZ = 0
    plot_z = [] # Distance
    plot_f = [] # Force
    plot_I = [] # Current
    plot_t = [] # SimTime 
    
    # Proceed until we've either reached the crossover point or exceeded the simulation sweep distance.distZ < Z_sweep
    while (forceMeasured > -0.1 and currI >= 0):
        
        try :
            mi_analyze()
            mi_loadsolution()
        except:
            input()
            
        mo_hidedensityplot()
        mo_hidecontourplot()

        #mo_showdensityplot(1,0,5,1E-5,"bmag")
        mo_groupselectblock(1)

        plot_t += [simTime]

        forceMeasured=mo_blockintegral(19)

        mi_selectgroup(1)
        mi_movetranslate(0,Z_step)
        distZ = distZ + Z_step

        if (forceMeasured > 0) :
            plot_z += [distZ]
            plot_f += [forceMeasured]

        
            acc += [forceMeasured/(M_p/1000)]
            
            simTime += [ simTime[len(simTime)-1] + (sqrt(vel[len(vel)-1]**2 + 2*acc[len(acc)-1]*1E-3)-vel[len(vel)-1])/acc[len(acc) - 1] ]
            
            if len(vel) > 0:
                vel += [vel[len(vel) - 1] +  acc[len(acc)-1]*(simTime[len(simTime) - 1] - simTime[len(simTime) - 2])]
            
            if len(pos) > 0:
                
                pos += [pos[len(pos) - 1] + Z_step]
                
                #pos += [(vel[len(vel) - 1] - vel[len(vel) - 2])/(simTime[len(simTime) - 1] - simTime[len(simTime) - 2])]
                #Z_step = pos[len(pos) - 1] / 1000

        
        currI = current[:, 0][findClosest(ts, simTime[len(simTime) - 1])]
        
        print("time [ms], pos[mm], I [A]:", simTime[len(simTime) - 1] * 1E3, Z_step, currI)
        plot_I += [currI]      
        mi_addcircprop("COIL", currI, 1)
    
    
    print("Speed gained:\t\t{:>10.3f} m/s \t {:>10.3f} fps".format(vel[len(vel) - 1] - startingVelocity, (vel[len(vel) - 1] - startingVelocity) * 3.28084))
    startingVelocity = vel[len(vel) - 1]    
    
    #mi_close()
    
    plt.cla()
    plt.plot(plot_z, plot_f)
    plt.title("Force by distance (Stage {0})".format(value))
    plt.ylabel('F [N]')
    plt.grid()
    plt.xlabel('Z [mm]')
    plt.savefig("{0}/Force_STAGE{1}.png".format(experimentName, value))
    plt.cla()
    
    plt.cla()
    plt.plot(simTime, vel)
    plt.title("Speed by time (Stage {0})".format(value))
    plt.ylabel('V [m/s]')
    plt.grid()
    plt.xlabel('T [s]')
    plt.savefig("{0}/Speed_STAGE{1}.png".format(experimentName, value))
    plt.cla()
    
    
    simFile = open('{0}/Raw_{1}_{2}.csv'.format(ExperimentName, simName, value), 'w')
    simFile.write('I [A], F [N], Acc [m/s*s], Vel [m/s], Z [mm], pos [mm]\n')
    for i in range(len(plot_z)):
        simFile.write('{},{},{},{},{},{},{}\n'.format(plot_t[i] * 1000,plot_I[i], plot_f[i], acc[i], vel[i], plot_z[i], pos[i]))
    simFile.close()
    # Cleanup
    mi_close()
    mo_close()
    return startingVelocity
    
    
ExperimentName = "{0}mmX{1}mm-FinalReport1".format(D_p, L_p)
os.mkdir(ExperimentName)
f = open('{}/Summary_{}.csv'.format(ExperimentName,ExperimentName), 'w')
openfemm()
hidepointprops()
for i in range(0,N_STAGES) :

    sim = runSim ("SIM_{1}_{0}".format(D_oc, ExperimentName), ExperimentName, f, i, startingVelocity)
    startingVelocity = sim
    f.close()
    
print("Final Velocity [m/s]\t{:>10.3f} \t [fps]\t{:>10.3f} ".format(startingVelocity, startingVelocity * 3.28084))
closefemm() # Goodnight!