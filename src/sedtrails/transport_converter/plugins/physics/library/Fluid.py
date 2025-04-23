"""
Class representing Fluid physical properties used by the transport converter.
"""

from dataclasses import dataclass, field
import numpy as np
from numpy import ndarray
from typing import Optional 



@dataclass
class D:
    # TODO: find a good name for this
    # TODO: define this class

    tau_mean : ndarray




@dataclass
class Fluid:
    """
    For example, water and wind.
    
    Attributes
    ----------
    abstau_mean : float
        absolute mean shear stress of the fluid [N/m^2].
    abstau_max : float
        absolute maximum shear stress of the fluid [N/m^2].
    U_mag : float
        total velocity magnitude of the fluid [m/s].
    ustar_mean : float
        mean friction velocity of the fluid [m/s].
    ustar_max : float
        maximum friction velocity of the fluid [m/s].
    """

    type: str  # for example water or wind
    d: D # object of type D
    abs_tau_mean: Optional[ndarray] = field(default = None, init = False)
    

    def __post_init__(self):
        self.compute_absolute_tau_mean() 


    # concrete method
    def compute_absolute_tau_mean(self) -> ndarray:
        """
        Calculates the absolute mean shear stress of the fluid [N/m^2].
        """
        self.abs_tau_mean = np.abs(self.d.tau_mean)
    
    
    def abstau_max(self, D) -> float:
        """
        Calculates the absolute maximum shear stress of the fluid [N/m^2].
        """
        return np.abs(D.tau_max)
    
    def U_mag(self, Ux, Uy) -> float:
        """
        Calculates the total velocity magnitude of the fluid [m/s].
        """
        return np.sqrt(Ux ** 2 + Uy ** 2)
    
    def ustar_mean(self, abstau_mean, S) -> float:
        """
        Calculates the mean friction velocity of the fluid [m/s].
        """
        return np.sqrt(abstau_mean / S.rhoFluid)

    def ustar_max(self, abstau_max, S) -> float:
        """
        Calculates the maximum friction velocity of the fluid [m/s].
        """
        return np.sqrt(abstau_max / S.rhoFluid)
    

if __name__ == '__main__':

    my_array = np.array([1,2,3])

    d = D(tau_mean=my_array)
    print(d.tau_mean)

    fluid = Fluid("water", d)

    print(fluid.abs_tau_mean)