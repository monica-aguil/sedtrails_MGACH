# For now I will start this as a function and later on I will shift it to a class
# if I see it is possible

import numpy as np
from numpy.typing import NDArray

from sedtrails.transport_converter.plugins.physics.library.Fluid import Fluid
from sedtrails.transport_converter.plugins.physics.library.Particle import Sand

# Create an instance of the class
Fluid_copy = Fluid()
Sand_copy = Sand()


# def Rfactor(
#         fluid: Fluid,
#         theta_max_Tracer
# ):
#     factor = theta_max_Tracer * fluid.p

def get_input(S, D, Ux, Uy):
    abstau_mean = Fluid_copy.abstau_mean(D)
    abstau_max = Fluid_copy.abstau_max(D)
    ustar_mean = Fluid_copy.ustar_mean(abstau_mean, S)
    ustar_max = Fluid_copy.ustar_max(abstau_max, S)
    U_mag = Fluid_copy.U_mag(Ux, Uy)

    theta_max_Tracer =  Sand_copy.theta_max_Tracer(abstau_max, S)
    Dstar_Tracer = Sand_copy.Dstar_Tracer(S)
    ratio = Sand_copy.ratio(S)
    theta_cr_Tracer = Sand_copy.theta_cr_Tracer(Dstar_Tracer)
    theta_cr_Tracer_exp = Sand_copy.theta_cr_Tracer_exp(theta_cr_Tracer, ratio)
    w_s_t = Sand_copy.w_s_t(Dstar_Tracer, S)
    rouse = Sand_copy.rouse(w_s_t, ustar_max)
    U_bed = Sand_copy.U_bed(ustar_mean, theta_cr_Tracer_exp, theta_max_Tracer)

    return abstau_mean, abstau_max, ustar_mean, ustar_max, U_mag, theta_max_Tracer, Dstar_Tracer, ratio, theta_cr_Tracer_exp, w_s_t, rouse, U_bed

# asumption
# require an instance of Fluid


def Rfacto(fluid: Fluid, particle: Particle) -> ndarray:

    """"""

    particle.theta_max_Tracer + 1





def Rfactor(
        theta_max_Tracer, 
        theta_cr_Tracer_exp, 
        rouse,
        U_mag,
        U_bed
        ) -> NDArray: # REVISAR PORQUE NO SE SI LOS OCUPO TODOS
    
    """
    Computes reduction factor R to obtain grain velocities according to Soulsby et al. (2011)

    Returns
    -------
    R (NDArray): Reduction factor
    """

    # theta_max_A: NDArray[np.float64],
    # theta_cr_A: float,
    # w_s_t: float,
    # u_star_m: NDArray[np.float64],
    # u_star_max: NDArray[np.float64],
    # Uc_mag: NDArray[np.float64],
    # Ub: NDArray[np.float64],

    Rb = np.zeros(U_bed.shape)
    Rs = np.zeros(U_bed.shape)
    R = np.zeros(U_bed.shape)

    for ii in range(0, theta_max_Tracer.shape[0]):
        for jj in range(0, theta_max_Tracer.shape[1]):
            if theta_max_Tracer[ii][jj] > theta_cr_Tracer_exp:  # [Equation 8]
                print(U_bed[ii][jj] / U_mag[ii][jj])
                Rb[ii][jj] = U_bed[ii][jj] / U_mag[ii][jj]
                if Rb[ii][jj] > 1:
                    Rb[ii][jj] = 1  # apply velocity limiter
            else:
                Rb[ii][jj] = 0

    # Rouse parameter (rouse < 2.5 indicates complete suspension)

    # Suspended load velocity reduction factor (R_suspload)
    for ii in range(0, theta_max_Tracer.shape[0]):
        for jj in range(0, theta_max_Tracer.shape[1]):
            if Rb[ii][jj] == 0:
                Rs[ii][jj] = 0
            else:
                # didn't find a good way of presenting this huge equation readable in python!!!!
                Rs[ii][jj] = np.multiply(
                    np.divide(
                        np.multiply(Rb[ii][jj], (1 - rouse[ii][jj])), (8 / 7 - rouse[ii][jj]) 
                    ),
                    np.divide(
                        (np.power((8 / 7 * Rb[ii][jj]), (8 - 7 * rouse[ii][jj])) - 1),
                        (np.power((8 / 7 * Rb[ii][jj]), (7 - 7 * rouse[ii][jj])) - 1),
                    ),
                )

            # Apply velocity limiter (grain velocity cannot exceed flow velocity)
            if Rs[ii][jj] > 1:
                Rs[ii][jj] = 1
            elif np.isnan(Rs[ii][jj]):  # I am not sure about this computation
                Rs[ii][jj] = 0

    # Reduction factor R
    for ii in range(0, theta_max_Tracer.shape[0]):
        for jj in range(0, theta_max_Tracer.shape[1]):
            if (
                rouse[ii][jj] < 2.5
            ):  # if all material is in suspension use Rs (for suspended load)
                R[ii][jj] = Rs[ii][jj]
            else:  # otherwise use Rb (for bed load)
                R[ii][jj] = Rb[ii][jj]
    return Rb, Rs, R 
