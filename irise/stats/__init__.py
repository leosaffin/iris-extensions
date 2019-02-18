"""
"""
import numpy as np
from properscoring import crps_ensemble


def crps(forecast, truth):
    r"""Continous ranked probability score

    .. math:
        CRPS = \int_{-\infty}^{\infty} [P(x) - P_a(x)]^2 dx

    Where :math:`P(x)` is the predicted cumulative probability distribution
    and :math:`P_a(x)` is the true cumulative probability distribution,
    which is typically represented by the Heaviside function.

    Args:
        forecast (iris.cube.Cube):

        truth (iris.cube.Cube): Same as forecast but without the
            `ensemble_member` coordinate

    Returns:

    """

    crps_ensemble(observations=truth, forecasts=forecast,
                  weights=None, issorted=False, axis=-1)
    return


def hellinger_distance(p, q):
    r"""Calculate the similarity between probability distributions p and q

    :math:`H(P,Q) = \frac{1}{\sqrt{2}} ||\sqrt{P} - \sqrt{Q}||_2`

    Args:
        p,q (iris.cube.Cube): Two probability distributions

    Returns:
        h: The Hellinger distance

    """

    h = ((p**0.5 - q**0.5)**2).mean()**0.5 / 2**0.5

    return h


def overlapping_coefficient(pdf1, pdf2):
    r"""Calculate the overlap between two probability density functions

    :math:`OVL = \sum_X (f_1(X), f_2(X)) dX`

    Args:
        pdf1:
        pdf2:

    Returns:

    """
    ovl = np.minimum(pdf1, pdf2, axis=0).sum()

    return ovl


def acc(forecast, analysis, climatology, **kwargs):
    """Calculate the anomaly correlation coefficient

    Args:
        forecast (iris.cube.Cube):

        analysis (iris.cube.Cube):

        climatology (iris.cube.Cube):

        kwargs:
            Extra arguments to pass to pearsonr (e.g. weights)
    """
    correlation = pearsonr(
        forecast - climatology, analysis - climatology, **kwargs)

    return correlation
