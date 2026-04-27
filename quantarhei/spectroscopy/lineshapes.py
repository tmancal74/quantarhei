#
# This module defines several normalized lineshapes for 1D and 2D spectroscopy
#
from __future__ import annotations

import numpy
from scipy import special


class Storage:
    def __init__(self, rel_tol: float = 1e-8) -> None:
        self.params: list[tuple[float, ...]] = []
        self.results: list[numpy.ndarray] = []
        self.positions: list[tuple[float, float]] = []
        self.rel_tol = rel_tol  # relative tolerance for parameter comparison
        self.lookup_count = 0
        self.ok_count = 0
        self.move_count = 0
        self.flip_count = 0
        self.lookup_on = False

    def _is_close(self, a: tuple[float, ...], b: tuple[float, ...]) -> bool:
        """Check if two 4-tuples are close within relative tolerance."""
        return all(numpy.isclose(a_i, b_i, rtol=self.rel_tol, atol=0.0) for a_i, b_i in zip(a, b))

    def lookup(self, param_tuple: tuple[float, ...]) -> int:
        """Return the index of param_tuple if found, else -1."""
        self.lookup_count += 1
        for idx, existing in enumerate(self.params):
            if self._is_close(param_tuple, existing):
                return idx
        return -1

    def store(self, param_tuple: tuple[float, ...], pos: tuple[float, float], matrix: numpy.ndarray) -> None:
        """Store new parameters and corresponding matrix."""
        self.params.append(param_tuple)
        self.results.append(matrix)
        self.positions.append(pos)

# Singleton instance
storage = Storage()





###############################################################################
#
#    1D absorptive lineshapes
#
###############################################################################

def gaussian(omega: numpy.ndarray, cent: float, delta: float) -> numpy.ndarray:
    """Normalized Gaussian line shape

    """
    return numpy.sqrt(numpy.log(2.0)/numpy.pi)\
                     *numpy.exp(-numpy.log(2.0)*((omega-cent)/delta)**2) \
                     /delta


def lorentzian(omega: numpy.ndarray, cent: float, gamma: float) -> numpy.ndarray:
    """Normalized Lorenzian line shape

    """
    return (gamma/numpy.pi)/((omega-cent)**2 + gamma**2)


def lorentzian_im(omega: numpy.ndarray, cent: float, gamma: float) -> numpy.ndarray:
    """Imaginary part of a normalized Lorenzian line shape

    """
    return 1j*((omega-cent)/numpy.pi)/((omega-cent)**2 + gamma**2)


def voigt(omega: numpy.ndarray, cent: float, delta: float, gamma: float = 0.0) -> numpy.ndarray:
    """Normalized Voigt line shape for absorption

    """
    z = (omega - cent + 1j*gamma)*numpy.sqrt(numpy.log(2.0))/delta

    return numpy.sqrt(numpy.log(2.0))*\
                      numpy.real(special.wofz(z)) \
                      /(numpy.sqrt(numpy.pi)*delta)


def cvoigt(omega: numpy.ndarray, cent: float, delta: float, gamma: float = 0.0) -> numpy.ndarray:
    """Complex normalized Voigt line shape

    """
    a = (delta**2)/(4.0*numpy.log(2))
    z = (gamma - 1j*(omega - cent))/(2.0*numpy.sqrt(a))


    return numpy.real(special.erfcx(z))*numpy.sqrt(numpy.pi/a)/2.0


###############################################################################
#
#    2D lineshapes
#
###############################################################################


def gaussian2D(omega1: numpy.ndarray, cent1: float, delta1: float, omega2: numpy.ndarray, cent2: float, delta2: float, corr: float = 0.0) -> numpy.ndarray:
    """Two-dimensional complex Gaussian lineshape

    """
    gamma1 = 0.0
    gamma2 = 0.0
    return voigt2D(omega1, cent1, delta1, gamma1,
                   omega2, cent2, delta2, gamma2, corr=corr)


def voigt2D(omega1: numpy.ndarray, cent1: float, delta1: float, gamma1: float,
            omega2: numpy.ndarray, cent2: float, delta2: float, gamma2: float, corr: float = 0.0) -> numpy.ndarray:
    """Two-dimensional complex Voigt lineshape

    Parameters
    ----------
    omega1, omega2 : real arrays
        Arrays of frequencies

    cent1, cent2 : real
        Center of the lineshape coordinates

    delta1, delta2 : real
        Gaussian linshape widths

    gamma1, gamma2 : real
        Exponential decays of the signal in t1 and t3 times (here denoted as 2)

    corr : real
        Correlation in the lineshape

    """
    N1 = omega1.shape[0]
    N2 = omega2.shape[0]

    if corr == 0.0:

        dat1 = cvoigt(omega1, cent1, delta1, gamma1)
        dat2 = cvoigt(omega2, cent2, delta2, gamma2)

        #data = numpy.zeros((N1, N2), dtype=COMPLEX)
        #
        #for k in range(N1):
        #    data[:, k] = dat1[k]*dat2[:]

        if storage.lookup_on:

            params = (delta1, delta2, gamma1, gamma2)

            indx = storage.lookup(params)
            if indx < 0:
                data = numpy.outer(dat2, dat1)
                storage.store(params, (cent1, cent2), data)
            else:
                data = storage.results[indx]
                pos = storage.positions[indx]
                if storage._is_close((cent1, cent2), pos):
                    print("Ok")
                    storage.ok_count += 1
                else:
                    if (numpy.abs(cent1) == numpy.abs(pos[0])) and (cent2 == pos[1]):
                        print("Just a flip", pos[0], cent1, pos[1], cent2)
                        storage.flip_count += 1
                    else:
                        print("Has to be moved by", pos[0], cent1, pos[1], cent2)
                        storage.move_count += 1

        else:
            data = numpy.outer(dat2, dat1)

        #print(indx, len(storage.params), storage.lookup_count)

    else:

        raise Exception("Not implemented yet")

    return data


def lorentzian2D(omega1: numpy.ndarray, cent1: float, gamma1: float, omega2: numpy.ndarray, cent2: float, gamma2: float, corr: float = 0.0) -> numpy.ndarray:
    """Two-dimensional complex Lorentzian lineshape

    """
    N1 = omega1.shape[0]
    N2 = omega2.shape[0]

    if corr == 0.0:

        dat1 = lorentzian(omega1, cent1, gamma1) + \
               lorentzian_im(omega1, cent1, gamma1)
        dat2 = lorentzian(omega2, cent2, gamma2) + \
               lorentzian_im(omega2, cent2, gamma2)

        #data = numpy.zeros((N1, N2), dtype=COMPLEX)
        #
        #for k in range(N1):
        #    data[k, :] = dat1[k]*dat2[:]
        data = numpy.outer(dat1,dat2)

    else:

        raise Exception("Not implemented yet")

    return data
