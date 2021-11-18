import numpy as np

import iris.plot as iplt


def get_contours(evaporation_tracer, threshold=1e-4):
    cs = iplt.contour(evaporation_tracer, [threshold])

    return cs.allsegs[0]


def is_closed(contour, threshold):
    """Checks that a contour is closed

    Checks that the final point along a contour is sufficiently close to the
    initial point on a countour to determine if it is closed

    Args:
        contour_section (np.Array):
            An array of coordinates for each point along the contour of shape
            (N, 2).

        threshold (scalar):

    Returns:
        True: If contour is closed

        False: If contour is open
    """
    return haversine(contour[0], contour[-1]) < threshold


def haversine(x1, x2):
    """ Calculate the great circle distance between two points on the earth
    (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(np.deg2rad, [x1[0], x1[1], x2[0], x2[1]])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    r = 6371  # Radius of earth in kilometers
    return c * r
