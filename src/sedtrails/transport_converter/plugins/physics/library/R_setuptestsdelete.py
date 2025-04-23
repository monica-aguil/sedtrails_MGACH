import numpy as np

# This is just a quick set-up for the remaining py-test for the R function

# abstau_mean = Fluid_copy.abstau_mean(D) # I need D.tau_mean
# abstau_max = Fluid_copy.abstau_max(D) # I need D.tau_max
# ustar_mean = Fluid_copy.ustar_mean(abstau_mean, S)
# ustar_max = Fluid_copy.ustar_max(abstau_max, S)
# U_mag = Fluid_copy.U_mag(Ux, Uy)

# theta_max_Tracer =  Sand_copy.theta_max_Tracer(abstau_max, S)
# Dstar_Tracer = Sand_copy.Dstar_Tracer(S)
# ratio = Sand_copy.ratio(S)
# theta_cr_Tracer = Sand_copy.theta_cr_Tracer(Dstar_Tracer)
# theta_cr_Tracer_exp = Sand_copy.theta_cr_Tracer_exp(theta_cr_Tracer, ratio)
# w_s_t = Sand_copy.w_s_t(Dstar_Tracer, S)
# rouse = Sand_copy.rouse(w_s_t, ustar_max)
# U_bed = Sand_copy.U_bed(ustar_mean, theta_cr_Tracer_exp, theta_max_Tracer)

rhoFluid = 1027 # kg/m3
rhoParticle = 2650 # kg/m3
g = 9.81 # m/s2
dTracer = 0.0002 # m
dBackground = 0.0002 # m
visc_kin = 1.36e-6 # m2/s

Ux = np.array([[-0.2975, -0.3079], [0.0562, 0.0692]]) # m/s From Natascia's input = timestep 146 rows 939 and 940 and 952 and 953
Uy = np.array([[0.0408, 0.0427], [0.0902, 0.0925]]) # m/s
tau_mean = np.array([[1, 1], [1, 1]]) # N/m^2 To correct as soon as Natascia sends me the file
tau_max = np.array([[1, 1], [1, 1]]) # N/m^2 To correct as soon as Natascia sends me the file

# Fluid properties 
abstau_mean = np.abs(tau_mean)
abstau_max = np.abs(tau_max)
U_mag = np.sqrt(Ux ** 2 + Uy ** 2) # m/s
ustar_mean = np.sqrt(abstau_mean / rhoFluid) # m/s
ustar_max = np.sqrt(abstau_max / rhoFluid) # m/s

# Particle properties
theta_max_Tracer = abstau_max / (g * (rhoParticle - rhoFluid) * dTracer)
Dstar_Tracer = (g * (rhoParticle / rhoFluid - 1) / (visc_kin ** 2)) ** (1/3) * dTracer
ratio = dTracer / dBackground
theta_cr_Tracer = 0.3 / (1 + 1.2 * Dstar_Tracer) + 0.055 * (1 - np.exp(-0.020 * Dstar_Tracer))
theta_cr_Tracer_exp = theta_cr_Tracer * np.sqrt(8 / (3 * (ratio ** 2) + 6 * ratio - 1)) * ((3.2260 * ratio) / 
                                                                                           (4 * ratio - 2 * (ratio + 1 - np.sqrt(ratio ** 2 + 2 * ratio - 1/3))))
w_s_t = (visc_kin / dTracer) * (np.sqrt(10.36 ** 2 + 1.049 * (Dstar_Tracer ** 3)) - 10.36)
rouse = w_s_t / (0.4 * ustar_max)
U_bed = 10 * ustar_mean * (1 - 0.7 * np.sqrt(theta_cr_Tracer_exp / theta_max_Tracer))


# Case 1: theta_max_Tracer < theta_cr_Tracer_exp (done!) R = 0

# Case 2: theta_max_Tracer > theta_cr_Tracer_exp & rouse > 2.5 
# Manually calculating Rb
# Rb = U_bed / U_mag
# Rb[Rb > 1] = 1
# print(Rb)

# Rs = ((Rb*(1-rouse))/(8/7-rouse))*((np.power((8/7*Rb),(8-7*rouse)) - 1)/(np.power((8/7*Rb),(7-7*rouse)) - 1))
# Rs[Rs > 1] = 1
# print(Rs)

# Function to compare matrices
def compare_matrices(matrix1, matrix2):
    return np.array_equal(matrix1, matrix2)

# Test function using pytest
def test_matrices_are_equal():

    Rb, Rs, R = Rfactor()

    # Assert matrix1 and matrix2 are equal
    assert compare_matrices(Rb, R) == True

    # Assert matrix1 and matrix3 are not equal
    assert compare_matrices(Rs, R) == False