import numpy as np
from numpy._typing._array_like import NDArray
import pytest
from sedtrails.transport_converter.plugins.physics.library.R_trial import Rfactor


@pytest.fixture
def sample_data_case_1():
    """
    Defining some random initial values to do test in which R comes as zeros
    """
    rhoFluid = 1027 # kg/m3
    rhoParticle = 2650 # kg/m3
    g = 9.81 # m/s2
    dTracer = 0.0002 # m
    dBackground = 0.0002 # m
    visc_kin = 1.36e-6 # m2/s

    ## Values to change by test
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

    return [theta_max_Tracer, theta_cr_Tracer_exp, rouse, U_mag, U_bed]

@pytest.fixture
def sample_data_case_2():
    """
    Defining some random initial values to do some test in which R = Rb because theta_max_Tracer > theta_cr_Tracer_exp & rouse > 2.5 
    """
    rhoFluid = 1027 # kg/m3
    rhoParticle = 2650 # kg/m3
    g = 9.81 # m/s2
    dTracer = 0.0002 # m
    dBackground = 0.0002 # m
    visc_kin = 1.36e-6 # m2/s

    ## Values to change by test
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

    return [theta_max_Tracer, theta_cr_Tracer_exp, rouse, U_mag, U_bed]

@pytest.fixture
def sample_data_case_3():
    """
    Defining some random initial values to do some test in which R = Rs because theta_max_Tracer > theta_cr_Tracer_exp & rouse < 2.5 
    """
    rhoFluid = 1027 # kg/m3
    rhoParticle = 2650 # kg/m3
    g = 9.81 # m/s2
    dTracer = 0.0002 # m
    dBackground = 0.0002 # m
    visc_kin = 1.36e-6 # m2/s

    ## Values to change by test
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

    return [theta_max_Tracer, theta_cr_Tracer_exp, rouse, U_mag, U_bed]

def test_R_equal_zeros(sample_data_case_1: NDArray):
    """
    When theta_max_A < theta_cr_A, R should be zero
    """
    theta_max_Tracer = sample_data_case_1[0]
    theta_cr_Tracer_exp = sample_data_case_1[1]
    rouse = sample_data_case_1[2]
    U_mag = sample_data_case_1[3]
    U_bed = sample_data_case_1[4]

    def all_zeros(arr_list):
        return all(np.all(arr == 0) for arr in arr_list)

    assert all_zeros(Rfactor(
        theta_max_Tracer, 
        theta_cr_Tracer_exp,
        rouse,
        U_mag,
        U_bed,
        ))
    
def test_R_equal_Rb(sample_data_case_2: NDArray):
    """
    When theta_max_A < theta_cr_A, R should be zero
    """
    theta_max_Tracer = sample_data_case_2[0]
    theta_cr_Tracer_exp = sample_data_case_2[1]
    rouse = sample_data_case_2[2]
    U_mag = sample_data_case_2[3]
    U_bed = sample_data_case_2[4]

    # Function to compare matrices
    def compare_matrices(matrix1, matrix2):
        return np.array_equal(matrix1, matrix2)
    
    Rb, Rs, R = Rfactor(
        theta_max_Tracer, 
        theta_cr_Tracer_exp,
        rouse,
        U_mag,
        U_bed,
        )

    assert compare_matrices(Rb, R) is True
    assert compare_matrices(Rs, R) is False

def test_R_equal_Rs(sample_data_case_3: NDArray):
    """
    When theta_max_A < theta_cr_A, R should be zero
    """
    theta_max_Tracer = sample_data_case_3[0]
    theta_cr_Tracer_exp = sample_data_case_3[1]
    rouse = sample_data_case_3[2]
    U_mag = sample_data_case_3[3]
    U_bed = sample_data_case_3[4]

    # Function to compare matrices
    def compare_matrices(matrix1, matrix2):
        return np.array_equal(matrix1, matrix2)
    
    Rb, Rs, R = Rfactor(
        theta_max_Tracer, 
        theta_cr_Tracer_exp,
        rouse,
        U_mag,
        U_bed,
        )

    assert compare_matrices(Rs, R) is True
    assert compare_matrices(Rb, R) is False
