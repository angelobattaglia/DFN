"""
Module containing functions that replicate some useful matlab routines

Created on 8 March 2019
@authors:
    * Francesco Della Santa (FDS), Politecnico di Torino, ITALY

Updates:
dd Mon YYYY:
    * ...

"""

import numpy as np
import scipy as sp  # Still not used


def cart2pol(x, y):
    """
    Transform Cartesian to polar coordinates (look at car2pol help on matlab)
    :param x:
    :param y:
    :return theta: radians of the angle between [x,y] and x axis (anti-clockwise with x=east, y=north)
    :return rho: radius, norm of [x,y]
    """
    theta = np.arctan2(y, x)
    rho = np.hypot(x, y)
    return theta, rho


def pol2cart(theta, rho):
    """
    Transform polar to Cartesian coordinates (look at the pol2cart help on matlab)
    :param theta: radians of the angle between [x,y] and x axis (anti-clockwise with x=east, y=north)
    :param rho: radius, norm of [x,y]
    :return: x, y
    """
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)
    return x, y


def cart2sph(x, y, z):
    """
    Transform Cartesian to spherical coordinates (look at the cart2sph help on matlab)
    :param x:
    :param y:
    :param z:
    :return az: radians of the angle between projection on (x,y) plane of [x,y,z] and x axis (anti-clockwise with
    respect to z)
    :return el: radians of the angle between [x,y,z] and (x,y) plane (el >/=/< 0 if z >/=/< 0). then -pi/2 <= el <= pi/2
    :return r: radius, norm of [x,y,z]
    """
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    el = np.arctan2(z, hxy)
    az = np.arctan2(y, x)
    return az, el, r


def sph2cart(az, el, r):
    """
    Transform spherical to Cartesian coordinates (look at the sph2cart help on matlab)
    :param az: radians of the angle between projection on (x,y) plane of [x,y,z] and x axis (anti-clockwise with
    respect to z)
    :param el: radians of the angle between [x,y,z] and (x,y) plane (el >/=/< 0 if z >/=/< 0). then -pi/2 <= el <= pi/2
    :param r: radius, norm of [x,y,z]
    :return: x, y, z
    """
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    return x, y, z


def rotx(ang):
    """
    Create a rotation matrix (3 x 3 np.array) in the 3D space with respect to x axis
    :param ang: radians. If ang <(>) 0, rotation is (anti-)clockwise with respect to x.
    :return: Rx: np.array of shape (3, 3) representing the rotation matrix
    """
    Rx = np.array(
        [
            [1, 0, 0],
            [0, np.cos(ang), -np.sin(ang)],
            [0, np.sin(ang), np.cos(ang)]
        ]
    )

    return Rx


def roty(ang):
    """
    Create a rotation matrix (3 x 3 np.array) in the 3D space with respect to y axis
    :param ang: radians. If ang <(>) 0, rotation is (anti-)clockwise with respect to y.
    :return: Rx: np.array of shape (3, 3) representing the rotation matrix
    """
    Ry = np.array(
        [
            [np.cos(ang), 0, np.sin(ang)],
            [0, 1, 0],
            [-np.sin(ang), 0, np.cos(ang)]
        ]
    )

    return Ry


def rotz(ang):
    """
    Create a rotation matrix (3 x 3 np.array) in the 3D space with respect to z axis
    :param ang: radians. If ang <(>) 0, rotation is (anti-)clockwise with respect to z.
    :return: Rz: np.array of shape (3, 3) representing the rotation matrix
    """
    Rz = np.array(
        [
            [np.cos(ang), -np.sin(ang), 0],
            [np.sin(ang), np.cos(ang), 0],
            [0, 0, 1]
        ]
    )

    return Rz





