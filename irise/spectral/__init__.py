import numpy as np


def reshape(spectra):
    """Rearrange a compressed 1d array of spherical harmonics to 2d

    Args:
        spectra (np.ndarray): 1 dimensional storage of 2d spherical modes

    Returns:
        np.ndarray:
            2-dimensional array of the reshaped input with zonal and meridional
            wavenumber coordinates

    """
    # Account for complex inputs as two dimensions
    if spectra.ndim == 2:
        spectra = spectra[:, 0] + spectra[:, 1]*1j
    if spectra.ndim != 1:
        raise ValueError('Spectra must be a 1-dimensional array')

    # Deduce truncation from shape
    trunc = find_trunc(len(spectra))

    # Zeros for output
    spectra_2d = np.zeros((trunc, trunc))
    idx0 = 0
    idx1 = trunc
    for i in range(trunc):
        spectra_2d[i, i:trunc] = spectra[idx0:idx1]
        idx0 += trunc - i
        idx1 += trunc - i - 1
    return spectra_2d


def find_trunc(n):
    """Solve quadratic: (T+1)(T+2)/2 = n

    Args:
        n (int): Number of elements of the flattened array of modes

    Returns:
         int: Inferred spectral truncation, T, of the cube
    """
    # Solve the quadratic equation
    trunc = (-3 + (np.sqrt(1 + 8*n)))/2

    # Round the result to nearest integer to compensate any floating-point error
    trunc = int(round(trunc))

    return trunc
