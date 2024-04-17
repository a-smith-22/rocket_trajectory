# TITLE: 1-D Vehicle Launch Trajectory
# AUTHOR: Andrew Smith
# DATE: 04/17/2024
# NOTES:
# - (4/16) Added stage mass fraction range to mass calculations
# - (4/16) Included multiple simulations across all MR values, print range of trajectories
# - (4/17) Improved user interface for increased functionality



# I. DEFINE PARAMETERS

# Import libraries
import math
import numpy as np
import time

# Rocket motor characteristics
O50SXL_P = [69.1, 9737000, 140802]             # stage 1, performance values -> [ burn time (s), total impulse (lbf-sec), burn time average thrust (lbf) ]
O50SXL_W = [35656, 33121, 1923, 545, 83, 2408] # stage 1, weight values (lbm) -> [ total loaded, propellent, case, nozzle, other, burnout ]
O50SXL_D = [50, 404]                           # stage 1, dimensions (in) -> [ motor diameter, motor length ]
O50XL_P =  [69.7, 2518000, 36096]              # stage 2, " "
O50XL_W =  [9520, 8650, 551, 240, 79, 824]     # stage 2, " "
O50XL_D =  [50, 122]                           # stage 2, " "
O38_P =    [67.7, 491000, 7246]                # stage 3, " "
O38_W =    [1966, 1699, 133, 91, 46, 243]      # stage 3, " "
O38_D =    [38, 53]                            # stage 3, " "

# Input parameters
print('TITLE: 1-D Vehicle Launch Trajectory (V1.1)') 
print('AUTHOR: Andrew Smith')
print('UPDATED: 04/16/2024')
print('')
print('========================================================')
print('')
print('Rocket motor mass fraction data (theoretical max): ')
print('  Stage 1 (MR): ' + str("{:.2f}".format( O50SXL_W[0]/O50SXL_W[5] )) ) # baseline mass fraction (mo/mf) -> use as minimum for stage mass fraction
print('  Stage 2 (MR): ' + str("{:.2f}".format( O50XL_W[0]/O50XL_W[5] )) )   # " "
print('  Stage 3 (MR): ' + str("{:.2f}".format( O38_W[0]/O38_W[5] )) )
print('')
print('========================================================')
print('')
print('Enter stage mass fraction (mo/mf) ranges as: min max')
MR_1 = [] # set up variables for input
MR_2 = [] # " "
MR_2 = [] # " "                
MR_1 = [float(item) for item in input("  Stage 1: ").split()] # stage 1 motor mass fraction range -> [min, max]
MR_2 = [float(item) for item in input("  Stage 2: ").split()] # stage 2 " "
MR_3 = [float(item) for item in input("  Stage 3: ").split()] # stage 3 " "
if len(MR_1) > 1 or len(MR_2) > 1 or len(MR_3) > 1:    # if multiple values entered, ask for step size
    dMR  = float(input("  Mass fraction step size: ")) # step size for mass fraction values
print('Enter simulation parameters')
dt = float( input("  timestep (s): ") )    # time step for simulation
t_max = float( input("  max time (s): ") ) # maximum time to compute simulation
print('')
print('========================================================')
print('')

# Constant values
lbm_to_kg = 0.453592 # mass
lbf_to_N = 4.44822   # force
in_to_m = 0.0254     # length

# System settings
t = 0.0     # current time (s)
h = 0.0     # current altitude (m)
v = 0.0     # current velocity (m/s)
a = 0.0     # current acceleration (m/s^2)
F = 0.0     # current force (N)
stage = 1   # current rocket stage 

# Vehicle characteristics
pl = 3000 # payload (kg)
cd = 0.02 # drag coefficient (initial value)



# II. GLOBAL FUNCTIONS

def array_conversion(vals, factors):
    # Convert an array of values using a single value or an array of scaling factors.
    new_vals = []
    if( type(factors) == list ): # multiple values provided
        if ( len(factors) == len(vals) ): 
            for i in range(len(vals)):
                    new_vals.append(vals[i] * factors[i])
        else:
            raise Exception("Must enter the same number of factors as input values.")
    else: # single value provided
        for i in range(len(vals)):
            new_vals.append(vals[i] * factors)
    return new_vals
    

def atm_density(alt):
    # Calculate the density of the atmosphere at altitude alt (in km) using linear interpolation. 
    
    # Define known data points -> values taken from textbook (Appendix 2)
    h = [0, 1, 3, 5, 10, 25, 50, 75, 100, 130, 160, 200, 300, 400, 600, 1000] # known altitude data points (in km)
    p = [1.2250, 1.1117, 9.0912e-1, 7.6312e-1, 4.1351e-1, 4.0084e-2, 1.0269e-3, 3.4861e-5, 5.604e-7, 8.152e-9, 1.233e-9, 2.541e-10, 1.916e-11, 2.803e-12, 2.137e-13, 3.561e-15] # density values at above altitudes (in kg/m^3)
    
    # Get index of value before point
    for i in range(len(h)):
        if h[i] <= alt:
            i1 = i # index of largest value below input altitude
    
    # Get index of value after point
    if i1 == 0:
        i2 = 1 # use first two points
    elif i1 == len(h)-1:
        i1 = len(h)-2 # use last two points
        i2 = i1 + 1
    else:
        i2 = i1 + 1 # use point after it if not an edge case
        
    # Compute linear interpolation
    p_lin = p[i1] + ((p[i2]-p[i1])/(h[i2]-h[i1])) * (alt-h[i1]) 
    
    return p_lin


def grav(alt):
    # Calculate the gravitational acceleration at altitude alt (in km) above the surface.
    R0 = 6371.000 # volumetric mean radius (km) -> https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
    g0 = 9.820    # mean surface gravity (m/s^2) -> see above link
    g = g0 * (R0 / (R0 + alt))**2
    return g


def mach(vel, alt):
    # Determine the Mach number based on the current rocket velocity (in m/s) and altitude (in km).
    # Gas properties
    gamma = 1.4     # specific heat ratio (air) -> assume constant
    R0 = 8.314      # ideal gas constant (J/K-mol)
    Mw = 28.96/1000 # average molecular weight of air (kg/mol) -> https://unacademy.com/content/question-answer/chemistry/what-is-the-molecular-weight-of-air/#:~:text=Air%20is%20a%20mixture%20of%20several%20gasses%20where%20the%20two,carbon%20dioxide%20of%20about%200.03%25.&text=We%20get%2028.96%20g%2Fmol,the%20molecular%20weight%20of%20Air.
    R = R0 / Mw     # specific gas constant for air (J/kg-K)
    
    # Define known data points -> values taken from textbook (Appendix 2)
    h = [0, 1, 3, 5, 10, 25, 50, 75, 100, 130, 160, 200, 300, 400, 600, 1000] # known altitude data points (in km)
    T = [288.150, 281.651, 268.650, 255.650, 223.252, 221.552, 270.650, 206.650, 195.08, 469.27, 696.29, 845.56, 976.01, 995.83, 999.85, 1000.00] # temperature data points (in K)
    
    # Get index of value before point
    i1 = 0
    for i in range(len(h)):
        if h[i] <= alt:
            i1 = i # index of largest value below input altitude
    
    # Get index of value after point
    if i1 == 0:
        i2 = 1 # use first two points
    elif i1 == len(h)-1:
        i1 = len(h)-2 # use last two points
        i2 = i1 + 1
    else:
        i2 = i1 + 1 # use point after it if not an edge case
        
    # Compute linear interpolation
    T_lin = T[i1] + ((T[i2]-T[i1])/(h[i2]-h[i1])) * (alt-h[i1])
    
    # Compute Mach number
    a = math.sqrt(gamma * R * T_lin) # speed of sound based on atmospheric temperature
    return (vel / a)
    

def drag_coeff(M):
    # Estimate the drag coefficient on the vehicle based on the Mach number. 
    
    # Values found Cd vs M graph on page 105 of textbook (zero angle of attack)
    Ms = [0, 0.192434211, 0.394736842, 0.587171053, 0.65625, 0.715460526, 0.764802632, 0.809210526, 0.863486842, 0.917763158, 0.967105263, 1.041118421, 1.090460526, 1.144736842, 1.174342105, 1.248355263, 1.292763158, 1.386513158, 1.549342105, 1.623355263, 1.751644737, 1.879934211, 2.077302632, 2.457236842, 2.930921053, 3.5625, 4.154605263, 4.702302632, 5.264802632, 5.496710526] # known Mach number values (approximated from Excel)
    Cd = [0.147040498, 0.147040498, 0.147040498, 0.148286604, 0.153271028, 0.150778816, 0.175700935, 0.190654206, 0.214330218, 0.251713396, 0.291588785, 0.357632399, 0.396261682, 0.416199377, 0.413707165, 0.395015576, 0.380062305, 0.352647975, 0.309034268, 0.290342679, 0.264174455, 0.249221184, 0.236760125, 0.224299065, 0.208099688, 0.186915888, 0.168224299, 0.15576324, 0.145794393, 0.144548287]  # known Cd values (see above)
    
    # Get index of value before point
    i1 = 0
    for i in range(len(Ms)):
        if Ms[i] <= M:
            i1 = i # index of largest value below input Mach number
    
    # Get index of value after point
    if i1 == 0:
        i2 = 1 # use first two points
    elif i1 == len(Ms)-1:
        i1 = len(Ms)-2 # use last two points
        i2 = i1 + 1
    else:
        i2 = i1 + 1 # use point after it if not an edge case
        
    # Compute linear interpolation
    Cd_lin = Cd[i1] + ((Cd[i2]-Cd[i1])/(Ms[i2]-Ms[i1])) * (abs(M)-Ms[i1])
    
    return Cd_lin


def print_progress(b):
    # Option to print summary statistics at each time step of the simulation. 
    if b == True:
        print("(" + str(i) + "/" + str(num_MRs) + ") -> " + "t: " + str("{:.1f}".format( t )) + "s, Stage: " + str(stage) + ", T: " + str("{:.3f}".format( T )) + "N, D: " + str("{:.3f}".format( D )) + "N, G: " + str("{:.0f}".format( G )) + "N, m-dot: " + str("{:.1f}".format( m_dot )) + "kg/s, m: " + str("{:.1f}".format( m )) + "kg, a: " + str("{:.3f}".format( a )) + "m/s^2, v: " + str("{:.1f}".format( v )) + "m/s, M: " + str("{:.3f}".format( M )) + ", h: " + str("{:.3f}".format( h/1000 )) + "km")
    
current_dots = -1 # number of dots displayed last
def progress_bar(n, N, d):
    # Print a progress bar for simulation n of N. Input contains number of previously displayed dots (prevents printing multiple times). 
    num_dots = 20                    # total number of dots to use (width of display)
    dots = round((n / N) * num_dots) # number of filled-in dots to display
    b1 = "□"
    b2 = "■"
    if dots != d: # prevents constant printing
        print(b2*dots + b1*(num_dots-dots))
    return dots



# III. PRE-PROCESSING

# Convert motor characteristics to SI units
O50SXL_P = array_conversion(O50SXL_P, [1.0, lbf_to_N, lbf_to_N]) # stage 1, convert performance values
O50SXL_W = array_conversion(O50SXL_W, lbm_to_kg)                 # stage 1, convert weight values
O50SXL_D = array_conversion(O50SXL_D, in_to_m)                   # stage 1, convert dimensions
O50XL_P =  array_conversion(O50XL_P, [1.0, lbf_to_N, lbf_to_N])  # stage 2, " "
O50XL_W =  array_conversion(O50XL_W, lbm_to_kg)                  # stage 2, " "
O50XL_D =  array_conversion(O50XL_D, in_to_m)                    # stage 2, " "
O38_P =    array_conversion(O38_P, [1.0, lbf_to_N, lbf_to_N])    # stage 3, " "
O38_W =    array_conversion(O38_W, lbm_to_kg)                    # stage 3, " "
O38_D =    array_conversion(O38_D, in_to_m)                      # stage 3, " "

# Mass parameters
# mass flow rates
dm_1 = O50SXL_W[1] / O50SXL_P[0] # stage 1 mass flow rate, assume constant = propellent mass / burn time
dm_2 = O50XL_W[1] / O50XL_P[0]   # stage 2 " "
dm_3 = O38_W[1] / O38_P[0]       # stage 3 " "
m_dot = dm_1                     # initial mass flow rate
# stage mass fractions
if len(MR_1) > 1:
    MR_1 = np.arange(MR_1[0], MR_1[-1] + dMR, dMR) # convert mass fractions to range of values
else:
    MR_1 = [MR_1[0]]                               # if only one value is entered, use just that value
if len(MR_2) > 1:
    MR_2 = np.arange(MR_2[0], MR_2[-1] + dMR, dMR) # " "
else:
    MR_2 = [MR_2[0]]                               # " "
if len(MR_3) > 1:
    MR_3 = np.arange(MR_3[0], MR_3[-1] + dMR, dMR) # " "
else:
    MR_3 = [MR_3[0]]                               # " "

# Thrust parameters
tb_total = O50SXL_P[0] + O50XL_P[0] + O38_P[0] # total burn time of all stages
print("Total burn time:         " + str("{:.0f}".format( tb_total )) + "s")

# Simulation data
i = 0                                                                # current simulation
num_MRs = len(MR_1) * len(MR_2) * len(MR_3)                          # number of simulations to perform
print("Number of simulations:   " + str("{:.0f}".format( num_MRs )) )
expected_time = ( num_MRs * (t_max) / dt ) * 2.3623e-5 + 1.2100e0    # simulations are O(n) -> use computations vs. execution time data for prediction
if expected_time > 1 and expected_time < 60: # in ms
    print("Expected execution time: " + str("{:.1f}".format( expected_time )) + "s")
elif expected_time > 60:                     # in min, sec
    print("Expected execution time: " + str("{:.0f}".format( math.floor(expected_time/60) )) + "min " + str("{:.0f}".format( expected_time%60 )) + "s")
else:                                        # in sec
    print("Expected execution time: " + str("{:.0f}".format( expected_time*1000 )) + "ms")
    
start_sim = input("Start simulation? [Y/N]: ")                       # allow user to halt operation if desired
if start_sim == "Y" or start_sim == "y":
    pass
else:                                                                # halt program
    print("")
    print("[Program terminated]")
    exit()                                      
    

# IV. SIMULATION

print('')
print('========================================================')
print('')
ask_print_progress = input("Display simulation progress? [Y/N]: ")
if ask_print_progress == "Y" or ask_print_progress == "y":
    print_progress_bool = True
else:
    print_progress_bool = False
print("Running simulation...")
start_time = time.time() # simulation start time

# Track trajectory data
h_max = [0.0] * num_MRs     # keep track of maximum height reached
t_at_max = [0.0] * num_MRs  # time to reach above maximum altitude
M_max = [0.0] * num_MRs     # maximum Mach number reached
M_min = [0.0] * num_MRs     # " " in descent
t_at_exo1 = [0.0] * num_MRs # time to reach space (pass Karman line)
t_at_exo2 = [0.0] * num_MRs # time to return to atmosphere
t_crash = [0.0] * num_MRs   # time of rocet crash (return to h=0)
MR_values = []              # list of stage mass fractions used in simulation

# Run simulation
for mr_1 in MR_1:
    for mr_2 in MR_2:
        for mr_3 in MR_3: # iterate over each stage mass fraction
            
            # Print computation progress
            current_dots = progress_bar(i, num_MRs, current_dots)
                       
            # Initialize all variables
            t = 0
            v = 0
            a = 0
            h = 0
            stage = 1
            
            # Compute vehicle initial mass
            m_s1_o = mr_1/(mr_1 - 1) * O50SXL_W[1] # stage one initial (wet) mass -> for MR = mf/mo, mo = mf/MR = MR/(MR-1) * mp
            m_s1_f = 1/(mr_1 - 1) * O50SXL_W[1]    # stage one final (dry) mass   -> for MR = mf/mo, mf = mo*MR = 1/(MR-1) * mp
            m_s2_o = mr_2/(mr_2 - 1) * O50XL_W[1]  # stage two " "
            m_s2_f = 1/(mr_2 - 1) * O50XL_W[1]     # stage two " "
            m_s3_o = mr_3/(mr_3 - 1) * O38_W[1]    # stage three " "
            m_s3_f = 1/(mr_3 - 1) * O38_W[1]       # stage three " "
            m0 = m_s1_o + m_s2_o + m_s3_o          # full rocket initial mass
            
            m = m0                               # initialize mass
            MR_values.append([mr_1, mr_2, mr_3]) # save current stage mass fractions
            
            # Run simulation
            while (t <= t_max):
                
                # Perform vehicle staging
                if t >= O50SXL_P[0]:              # stage 1 has burned out -> engage stage 2
                    stage = 2
                if t >= O50SXL_P[0] + O50XL_P[0]: # stage 2 has burned out -> engage stage 3
                    stage = 3
                
                # Set mass flow rate
                if stage == 1:
                    m_dot = dm_1
                elif stage == 2:
                    m_dot = dm_2
                elif stage == 3:
                    m_dot = dm_3
                if t > tb_total:
                    m_dot = 0 # no mass flow after final burnout time is reached
                    
                # Determine current mass
                if stage == 1:
                    m = m0 - m_dot*t                                                          # subtract expelled propellent
                elif stage == 2:
                    m = m0 - m_s1_f - m_dot*(t - O50SXL_P[0])                                 # " " + subtract above stage burnout mass
                elif stage == 3:
                    m = m0 - m_s1_f - m_s2_f - m_dot*(t - O50SXL_P[0] - O50XL_P[0])           # " "
                if t > tb_total:
                    m = m0 - m_s1_f - m_s2_f - m_s3_f                                         # subtract all stage burnout masses
                
                
                # Find drag force
                # cross sectional area
                if stage == 1:
                    A = math.pi/4 * O50SXL_D[0]**2
                elif stage == 2:
                    A = math.pi/4 * O50XL_D[0]**2
                elif stage == 3:
                    A = math.pi/4 * O38_D[0]**2
                # drag coefficient
                M = mach(v, h/1000) # get mach number
                cd = drag_coeff(M)  # get drag coefficient
                # overall drag
                if t > t_at_max[i]: # return direction -> drag acts upwards
                    D = -0.5 * atm_density(h/1000) * v**2 * A * cd
                else:            # ascension -> drag acts downwards
                    D = 0.5 * atm_density(h/1000) * v**2 * A * cd
               
                
                # Find gravitational force
                G = grav(h/1000)*m # use current altitude and mass
                    
                # Calculate thrust
                if stage == 1:
                    T = O50SXL_P[2] # assume constant thrust
                elif stage == 2:
                    T = O50XL_P[2]  # " "
                elif stage == 3:
                    T = O38_P[2]    # " "
                if t >= tb_total: # all rockets have burned out
                    T = 0         # no thrust produced after this point
               
                # Update physics states
                F = T - G - D # force = thrust - gravity - drag
                a = F / m     # Newton's 2nd law
                v += a*dt     # find velocity based on acceleration
                h += v*dt     # increase height 
                t += dt       # update time
                
                
                # Track state variables
                # maximum Mach numbers
                if M > M_max[i]:
                    M_max[i] = M    # save highest Mach number
                elif M < M_min[i]:
                    M_min[i] = M    # " " in return direction (negative)
                    
                # maximum altitude 
                if h > h_max[i]:
                    h_max[i] = h       # update peak position value (t, h)
                    t_at_max[i] = t    # " "
                elif h <= 0 and t > 0:
                    h = 0          # rocket has crashed
                    v = 0          # " "
                    t_crash[i] = t # save crash time
                    break
                
                # exo-atmospheric flight
                if h < 100e3: # Karman line = 100km
                    if t <= t_at_max[i]: # first time reaching point
                        t_at_exo1[i] = t
                elif h > 100e3:
                        t_at_exo2[i] = t
                        
                
                # Print current data (optional)
                print_progress(print_progress_bool) # toggle on/off to print simulation progress at each timestep (turn off for speed)
                

            # Update simulation count
            i += 1


# V. DISPLAY RESULTS 

print("Done.")
end_time = time.time()               # simulation end time
elapsed_time = end_time - start_time # total computation time (wall time) for simulation
if elapsed_time > 1:
    print("Execution time: " + str("{:.1f}".format( elapsed_time )) + "s ")
else:
    print("Execution time: " + str("{:.0f}".format( elapsed_time*1000 )) + "ms ")

print('')
print('========================================================')
print('')

max_alt = max(h_max)         # highest peak altitude achieved
min_alt = min(h_max)         # lowest peak altitude achieved
index = h_max.index(max_alt) # get index of best performing flight

print("Simulation Data:")
print("  Stage mass fractions: " + "[" + str("{:.1f}".format( MR_1[0] )) + "-" + str("{:.1f}".format( MR_1[-1] )) + "], [" + str("{:.1f}".format( MR_2[0] )) + "-" + str("{:.1f}".format( MR_2[-1] )) + "], [" + str("{:.1f}".format( MR_3[0] )) + "-" + str("{:.1f}".format( MR_3[-1] )) + "]") # display mass fraction ranges
print("  Tracjectory Range:    " + str("{:.3f}".format( min_alt/1000 )) + "-" + str("{:.3f}".format( max_alt/1000 )) + "km" )
print("")

print("Best Flight Data: ")
print("  Stage mass fractions: " + str("{:.2f}".format( MR_values[index][0] )) + ", " + str("{:.2f}".format( MR_values[index][1] )) + ", " + str("{:.2f}".format( MR_values[index][2] ))) # stage mass fractions used
print("  Maximum altitude:     " + str("{:.3f}".format( max_alt/1000 )) + "km")                                                                                                           # primary output of simulation
print("  Apogee time:          " + str("{:.2f}".format( t_at_max[index] )) + "s")                                                                                                         # time to reach above altitude
print("  Mach:                 " + str("{:.3f}".format( M_max[index] )) + " (" + str("{:.3f}".format( M_min[index] )) + ")")                                                              # maximum Mach number on ascension (and descension)
print("  Space flight:         " + str("{:.2f}".format( t_at_exo1[index] )) + "s" + " (" + str("{:.2f}".format( t_at_exo2[index] )) + "s)")                                               # time to escape (and return to) atmosphere
print("  Crash time:           " + str("{:.2f}".format( t_crash[index] )) + "s")                                                                                                          # time to return to starting height