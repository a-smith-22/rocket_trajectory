# TITLE: 1-D Vehicle Launch Trajectory
# AUTHOR: Andrew Smith
# DATE: 04/15/2024
# NOTES:



# I. DEFINE PARAMETERS

# Import libraries
import math

# Input parameters
print('TITLE: 1-D Vehicle Launch Trajectory (V1.0)') 
print('AUTHOR: Andrew Smith')
print('UPDATED: 04/16/2024')
print('')
print('===========================================')
print('')
print('Enter stage mass fraction ranges as [min, max]')
MR_1 = input("  Stage 1: ")                 # stage 1 motor mass fraction range -> [min, max]
MR_2 = input("  Stage 2: ")                 # stage 2 " "
MR_3 = input("  Stage 3: ")                 # stage 3 " "
dMR  = input("  Mass fraction step size: ") # step size for mass fraction values
print('Enter simulation parameters')
dt = float( input("  timestep (s): ") )    # time step for simulation
t_max = float( input("  max time (s): ") ) # maximum time to compute simulation
print('')
print('===========================================')
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
    Cd_lin = Cd[i1] + ((Cd[i2]-Cd[i1])/(Ms[i2]-Ms[i1])) * (M-Ms[i1])
    
    return Cd_lin


def print_progress(b):
    # Option to print summary statistics at each time step of the simulation. 
    if b == True:
        print("t: " + str("{:.1f}".format( t )) + "s, Stage: " + str(stage) + ", T: " + str("{:.3f}".format( T )) + "N, D: " + str("{:.3f}".format( D )) + "N, G: " + str("{:.0f}".format( G )) + "N, m-dot: " + str("{:.1f}".format( m_dot )) + "kg/s, m: " + str("{:.1f}".format( m )) + "kg, a: " + str("{:.3f}".format( a )) + "m/s^2, v: " + str("{:.1f}".format( v )) + "m/s, M: " + str("{:.3f}".format( M )) + ", h: " + str("{:.3f}".format( h/1000 )) + "km")
    



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
# initial mass
m0 = pl + O50SXL_W[0] + O50XL_W[0] + O38_W[0] # initial loaded mass
m = m0                                        # current mass (kg) -> set as initial mass
print("Initial vehicle mass: " + str("{:.0f}".format( m )) + "kg")
# mass flow rates
dm_1 = O50SXL_W[1] / O50SXL_P[0] # stage 1 mass flow rate, assume constant = propellent mass / burn time
dm_2 = O50XL_W[1] / O50XL_P[0]   # stage 2 " "
dm_3 = O38_W[1] / O38_P[0]       # stage 3 " "
m_dot = dm_1                     # initial mass flow rate

# Thrust parameters
tb_total = O50SXL_P[0] + O50XL_P[0] + O38_P[0] # total burn time of all stages
print("Total burn time: " + str("{:.0f}".format( tb_total )) + "s")


# IV. SIMULATION

print('')
print('===========================================')
print('')
print("Running simulation...")

# Track trajectory data
h_max = 0.0     # keep track of maximum height reached
t_at_max = 0.0  # time to reach above maximum altitude
M_max = 0.0     # maximum Mach number reached
M_min = 0.0     # " " in descent
t_at_exo1 = 0.0 # time to reach space (pass Karman line)
t_at_exo2 = 0.0 # time to return to atmosphere

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
        
    # Determine current mass
    if stage == 1:
        m = m0 - m_dot*t                                                          # subtract expelled propellent
    elif stage == 2:
        m = m0 - O50SXL_W[5] - m_dot*(t - O50SXL_P[0])                            # " " + subtract above stage burnout mass
    elif stage == 3:
        m = m0 - O50SXL_W[5] - O50XL_W[5] - m_dot*(t - O50SXL_P[0] - O50XL_P[0])  # " "
    if t > tb_total:
        m = m0 - O50SXL_W[5] - O50XL_W[5] - O38_W[5]                              # subtract all stage burnout masses
    
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
    if t > t_at_max: # return direction -> drag acts upwards
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
    if M > M_max:
        M_max = M    # save highest Mach number
    elif M < M_min:
        M_min = M    # " " in return direction (negative)
        
    # maximum altitude 
    if h > h_max:
        h_max = h    # update peak position value (t, h)
        t_at_max = t # " "
    elif h <= 0 and t > 0:
        h = 0     # rocket has crashed
        v = 0     # " "
        break
    
    # exo-atmospheric flight
    if h < 100e3: # Karman line = 100km
        if t <= t_at_max: # first time reaching point
            t_at_exo1 = t
    elif h > 100e3:
            t_at_exo2 = t
    
    
    # Print current data (optional)
    print_progress(False) # toggle on/off to print simulation progress at each timestep (turn off for speed)


# V. DISPLAY RESULTS

print("Done.")
print("")
print("Maximum altitude: " + str("{:.3f}".format( h_max/1000 )) + "km") # primary output of simulation
print("Apogee time: " + str("{:.2f}".format( t_at_max )) + "s")         # time to reach above altitude
print("Mach: " + str("{:.3f}".format( M_max )) + " (" + str("{:.3f}".format( M_min )) + ")")  # maximum Mach number on ascension (and descension)
print("Space flight: " + str("{:.2f}".format( t_at_exo1 )) + "s" + " (" + str("{:.3f}".format( t_at_exo2 )) + "s)")                 # time to escape (and return to) atmosphere
print("Crash time: " + str("{:.2f}".format( t )) + "s")                 # time to return to starting height