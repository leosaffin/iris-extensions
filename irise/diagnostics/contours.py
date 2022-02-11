from shapely.geometry.polygon import Polygon

import numpy as np

import iris.plot as iplt


def get_contours(evaporation_tracer, threshold=1e-4):
    cs = iplt.contour(evaporation_tracer, [threshold])

    return cs.allsegs[0]


def as_array(polygons):
    return np.array(polygons[0].exterior.coords)


def filter_invalid(contours, maximum_distance):
    """Remove contours that are not closed or contained within other contours

    Args:
        contours:
        maximum_distance:

    Returns:

    """
    # Sort contours by distance
    contours = sorted([contour for contour in contours if len(contour) > 2], key=contour_length)
    polygons = [Polygon(contour) for contour in contours]

    # Check whether the shorter contours are contained by the longer contours
    to_be_removed = []
    for n, poly1 in enumerate(polygons):
        for m, poly2 in enumerate(polygons):
            if n != m and poly2.contains(poly1):
                to_be_removed.append(n)

    # Remove open contours
    for n, contour in enumerate(contours):
        if not is_closed(contour, threshold=maximum_distance):
            to_be_removed.append(n)

    for n in reversed(sorted(set(to_be_removed))):
        del contours[n]

    return contours


def contour_length(points):
    """Contour length in kilometres

    Args:
        points: Nx2 array of longitude and latitude points around a contour (degrees)

    Returns:
        float: The total length of the contour in kilometres
    """
    conlen = haversine(points[-1], points[0])
    for n in range(len(points) - 1):
        conlen += haversine(points[n], points[n+1])

    return conlen


def is_closed(contour, threshold):
    """Checks that a contour is closed

    Checks that the final point along a contour is sufficiently close to the
    initial point on a countour to determine if it is closed

    Args:
        contour (np.Array):
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
