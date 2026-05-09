import numpy

from .types import check_numpy_array

X = (1.0, 0.0, 0.0)
Y = (0.0, 1.0, 0.0)
Z = (0.0, 0.0, 1.0)
Dxy = (1.0, 1.0, 0.0) / numpy.sqrt(2)
Axy = (-1.0, 1.0, 0.0) / numpy.sqrt(2)


def normalize2(vec: numpy.ndarray, norm: float = 1.0) -> numpy.ndarray:  # type: ignore[explicit-any]
    """Normalize a vector to a specified magnitude.

    Parameters
    ----------
    vec : numpy.ndarray
        Input vector to normalize.
    norm : float, optional
        Target magnitude. Default is ``1.0``.

    Returns
    -------
    numpy.ndarray
        Vector rescaled so that its Euclidean norm equals ``norm``.
    """
    vec = check_numpy_array(vec)
    vel = numpy.sqrt(numpy.dot(vec, vec))
    out = (vec / vel) * norm
    return out


def norm(vec: numpy.ndarray) -> float:  # type: ignore[explicit-any]
    """Return the Euclidean norm of a vector.

    Parameters
    ----------
    vec : numpy.ndarray
        Input vector.

    Returns
    -------
    float
        Euclidean norm (square root of the dot product of the vector with
        itself).
    """
    vel = numpy.sqrt(numpy.dot(vec, vec))
    return vel
