from math import radians, cos, sin, asin, sqrt
import numpy as np


def get_contour_verts(cn):
    """
    Args:
        cn: A set of contours returned from matplotlib.pyplot.contour

    Returns:
        contours (list): A list of contour lines corresponding to each value
            contoured. For each values the list will contain a list of
            contours. Each individual contour is an array of coordinates for
            each point along the contour of shape (N, 2).
    """
    contours = []
    # Loop over each differently valued contour in cn
    for cc in cn.collections:
        paths = []
        # Loop over each separated contour of the individual value
        for pp in cc.get_paths():
            xy = []
            # for each segment of that section
            for vv in pp.iter_segments():
                xy.append(vv[0])
            paths.append(np.vstack(xy))
        contours.append(paths)

    return contours


def find_longest_contour(contours):
    """Finds the longest contour from a set of contours

    Args:
        contours (list):
            A list of contours. Each individual contour is an array of
            coordinates for each point along the contour of shape (N, 2).

    Returns:
        imax (int): The index of the longest contour

        len_max (int): The length of the longest contour
    """
    # Get a list of the lengths of each contour
    lengths = np.array([contour_length(x) for x in contours])

    # Index of longest contour
    imax = lengths.argmax()

    # Value of longest contour
    len_max = lengths[imax]

    return imax, len_max


def is_closed_contour(contour_section, threshold=100):
    """Checks that a contour is closed

    Checks that the final point along a contour is sufficiently close to the
    initial point on a countour to determine if it is closed

    Args:
        contour_section (np.Array):
            An array of coordinates for each point along the contour of shape
            (N, 2).

    Returns:
        True: If contour is closed

        False: If contour is open
    """
    return haversine(contour_section[0], contour_section[-1]) < threshold


def get_tropopause_contour(contours, max_iterations=10):
    """Find the tropopause contour that wraps around the globe

    Args:
        contours (list):
            A list of contours. Each individual contour is an array of
            coordinates for each point along the contour of shape (N, 2).

    """
    # Find the longest contour
    imax = find_longest_contour(contours)[0]
    tropopause_contour = contours[imax]

    joined_contours = [imax]

    iter_no = 1
    while (not(is_closed_contour(tropopause_contour)) and
           iter_no < max_iterations):
        # If the original longest contour is not closed,
        # search for continuation contours and concatenate them
        for ii in np.delete(np.arange(len(contours)), joined_contours):
            tropopause_contour, joined = join_contours(
                tropopause_contour, contours[ii])
            if joined:
                joined_contours.append(ii)

        iter_no = iter_no + 1

    return tropopause_contour


def join_contours(long_contour, short_contour, threshold=100):
    # End of large contour joins start of small contour
    if haversine(long_contour[-1], short_contour[0]) < threshold:
        long_contour = np.append(long_contour, short_contour, axis=0)

    # Start of large contour joins end of small contour
    elif haversine(long_contour[0], short_contour[-1]) < threshold:
        long_contour = np.append(short_contour, long_contour, axis=0)

    # Start of large contour joins start of small contour
    elif haversine(long_contour[0], short_contour[0]) < threshold:
        long_contour = np.append(
            np.flipud(long_contour), short_contour, axis=0)

    # End of large contour joins end of small contour
    elif haversine(long_contour[-1], short_contour[-1]) < threshold:
        long_contour = np.insert(
            short_contour, np.flipud(long_contour), axis=0)

    else:
        long_contour, False

    return long_contour, True


def increase_circuit_resolution(points, resolution):
    """Add points around the circuit loop
    """
    # Loop over all points around the circuit.
    n = 0
    while n < len(points):
        # Allow the loop to connect the final point with the first point
        # i.e np1 = 0
        np1 = (n+1) % len(points)

        # Calculate distances from current point to next point and
        # Insert a point if the next point is further away than the required
        # resolution
        distance = haversine(points[n], points[np1])
        if distance > resolution:
            # Put the next point along the same line
            dlon = points[np1, 0] - points[n, 0]
            dlat = points[np1, 1] - points[n, 1]

            new_lon = points[n, 0] + (resolution / distance) * dlon
            new_lat = points[n, 1] + (resolution / distance) * dlat

            points = np.insert(points, np1, [new_lon, new_lat], axis=0)

        # Always jump to the next point. This will either be the inserted point
        # or the next point along that is within the required resolution
        n += 1

    return points


def contour_length(points):
    """Contour length in kilometres

    Args:
        points: Nx2 array of longitude and latitude points around a contour (degrees)

    Returns:
        float: The total length of the contour in kilometres
    """
    conlen = 0
    for n in range(len(points) - 1):
        conlen += haversine(points[n], points[n+1])

    return conlen


def haversine(x1, x2):
    """ Calculate the great circle distance between two points on the earth
    (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [x1[0], x1[1], x2[0], x2[1]])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * asin(sqrt(a))
    r = 6371  # Radius of earth in kilometers
    return c * r
