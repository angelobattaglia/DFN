"""
Module useful if you need some probabilities distributions not implemented yet in other modules or packages

Created on 15 March 2019
@authors:
    * Francesco Della Santa (FDS), Politecnico di Torino, ITALY

Updates:
dd Mon YYYY:
    * ...

"""

import numpy as np
import scipy as sp
import scipy.integrate as sp_int
import warnings as warn
import matlab_clones as mc


def create_uniform1D_pdf(a, b):
    def uniform_pdf(inputs):
        """
        Uniform pdf in [a, b]
        :param inputs: list or 1D np array
        :return: 1D np array of pdf values corresponding to inputs
        """
        return ((np.array(inputs) >= a) * (np.array(inputs) <= b)) * (1 / (b - a))

    return uniform_pdf


def create_uniform1D_cdf(a, b):
    def uniform_cdf(inputs):
        """
        Uniform cdf in [a, b]
        :param inputs: list or 1D np array
        :return: 1D np array of cdf values (that is P(X <= inputs[i])) corresponding to inputs
        """

        numerator_array = (np.zeros(np.array(inputs).shape) +
                           ((np.array(inputs) >= a) * (np.array(inputs) <= b)) * np.array(inputs) +
                           b * (np.array(inputs) > b))

        return abs(((numerator_array >= a) * (numerator_array - a)) / (b - a))

    return uniform_cdf


class ContinuousDistribution1D:
    """
    Class of objects that describes continuous distributions in 1D.
    """
    def __init__(self, name='noname_dist',
                 pdf=create_uniform1D_pdf(0, 1), cdf=create_uniform1D_cdf(0, 1),
                 **options):
        """
        Initialization method for the Distribution class
        :param name: string (default noname_dist)
        :param cdf: a python function that describes the cumulative density function (default uniform in [0,1])
        :param pdf: a python function that describes the probability density function (default uniform in [0,1])
        :param options:
            uniform_default: boolean.
            If False, the object is not the unif. dist. in [0,1] even if cdf and pdf are None. (default True)
            mean: mean/expected value of the distribution
            median: median value of the distribution
            inv_cdf: a python function that describes the inverse of the probability density function (default None).
        """
        self.name = name

        if (cdf is None) and (pdf is None) and (options.get('uniform_default', True)):
            self.cdf = create_uniform1D_cdf(0, 1)
            self.pdf = create_uniform1D_pdf(0, 1)
            self.mean = 0.5
            self.median = 0.5
            self.mode = None
            self.var = 1 / 12

            def identity(x):
                return x

            self.inv_cdf = identity

        elif cdf is None:
            self.pdf = pdf

            def built_cdf(x):
                return sp_int.quad(self.pdf, -sp.inf, x)[0]

            self.cdf = built_cdf
            self.mean = options.get('mean')
            self.median = options.get('median')
            self.mode = options.get('mode')
            self.var = options.get('var')
            self.inv_cdf = options.get('inv_cdf')

        else:
            self.pdf = pdf
            self.cdf = cdf
            self.mean = options.get('mean')
            self.median = options.get('median')
            self.mode = options.get('mode')
            self.var = options.get('var')
            self.inv_cdf = options.get('inv_cdf')

    def sample(self, nsamples=1):
        """
        Function that extract samples of the distribution
        :param nsamples: integer value (default 1)
        :return: array of samples
        """
        if self.inv_cdf is None:
            warn_msg = """
            The inverse of the CDF has not been defined. 
            Please define this function and assign it to the attribute "inv_cdf" of the object.
            """
            return warn.warn(warn_msg)
        else:
            return self.inv_cdf(np.random.rand(nsamples))


class VonMisesFisher:
    """
    Class of objects that describes the Von Mises - Fisher distribution on the sphere (radius 1),
    varying its concentration parameter k and its mode_vector
    """
    def __init__(self, k=10, mode_vector=np.array([[0.], [0.], [1.]])):
        """
        Initialization method
        :param k: concentration parameter (integer; default 10)
        :param mode_vector: mode vector (np.array of shape (3, 1); default [0,0,1])
        """
        self.k = k
        self._ck = sp.sinh(k) * (2 / k)
        self.mode_vector = mode_vector

        tuple_mode_vector = tuple(self.mode_vector[:, 0].tolist())
        az, el, r = mc.cart2sph(*tuple_mode_vector)
        self.mode_vector_sph = np.array([[az], [el], [r]])

        if (mode_vector == np.array([[0.], [0.], [1.]])).all():
            self.rot = np.identity(3)
        else:
            rot1 = mc.roty((np.pi / 2) - el)
            rot2 = mc.rotz(az)
            self.rot = rot2 @ rot1

        def w_mvf_pdf(x):
            return sp.exp(np.array(x, ndmin=1) * self.k) / self._ck

        def w_mvf_cdf(x):
            return (sp.exp(self.k * np.array(x, ndmin=1)) - sp.exp(-self.k)) / (self.k * self._ck)

        def w_mvf_inv_cdf(x):
            return sp.log(sp.exp(-self.k) + self.k * self._ck * np.array(x, ndmin=1)) / self.k

        self.w = ContinuousDistribution1D(name='w_vmf', pdf=w_mvf_pdf, cdf=w_mvf_cdf, inv_cdf=w_mvf_inv_cdf)

    def __v_mvf_sampling(self, nsamples=1):
        theta_unif = np.random.rand(nsamples) * 2 * np.pi
        return np.array([np.cos(theta_unif), np.sin(theta_unif)])

    def sample(self, nsamples=1):
        """
        method that returns a matrix 3 x nsamples of vectors extracted from the distribution
        :param nsamples: number of samples to extract (integer; at least 1)
        :return: vmf_samples
        """
        w_samples = self.w.sample(nsamples)
        v_samples = self.__v_mvf_sampling(nsamples)

        vmf_samples = np.vstack([np.sqrt(1 - w_samples ** 2) * v_samples,
                                 w_samples])

        if (self.mode_vector != np.array([[0.], [0.], [1.]])).any():
            vmf_samples = self.rot @ vmf_samples

        return vmf_samples


# ------------ RANDOM VON MISES - FISHER DISTRIBUTION (START) ------------
mi_rand = np.random.rand(3, 1)
mi_rand = mi_rand / np.linalg.norm(mi_rand)
k_rand = 10 + 90 * np.random.rand()

vmf_rand = VonMisesFisher(k=k_rand, mode_vector=mi_rand)
# ------------ RANDOM VON MISES - FISHER DISTRIBUTION (END) ------------


class PowerLawBounded(ContinuousDistribution1D):
    """
    Class of objects that describes bounded power law distributions in 1D.
    """
    def __init__(self, name='noname_powlaw',
                 alpha=5,
                 radius_l=5, radius_u=30,
                 **options):
        """

        :param name:
        :param alpha: float number greater than 1 (default 2)
        :param radius_l: float number representing the minimum value/radius (default 1)
        :param radius_u: float number representing the maximum value/radius (default 10)
        :param options:
            mean: mean/expected value of the distribution
            median: median value of the distribution
            inv_cdf: a python function that describes the inverse of the probability density function (default None).
        """
        self.alpha = alpha
        self.radius_l = radius_l
        self.radius_u = radius_u

        alpha_exp = 1 - self.alpha

        def pl_pdf(r):
            pdf_expression = (alpha_exp /
                              (self.radius_u ** alpha_exp - self.radius_l ** alpha_exp)) * \
                             (np.array(r, ndmin=1) ** (-self.alpha))
            pdf_ret = pdf_expression * (np.array(r, ndmin=1) >= self.radius_l) * (np.array(r, ndmin=1) <= self.radius_u)
            pdf_ret = np.array(pdf_ret, ndmin=1)
            return pdf_ret

        def pl_cdf(r):
            cdf_expression = (np.array(r, ndmin=1) ** alpha_exp - self.radius_l ** alpha_exp) / \
                             (self.radius_u ** alpha_exp - self.radius_l ** alpha_exp)
            cdf_ret = cdf_expression
            cdf_ret[np.array(r, ndmin=1) <= self.radius_l] = 0.
            cdf_ret[np.array(r, ndmin=1) >= self.radius_u] = 1.
            return cdf_expression

        def pl_inv_cdf(x):
            return self.radius_l * (1 + np.array(x) *
                                    ((self.radius_u / self.radius_l) ** alpha_exp - 1)) ** (1 / alpha_exp)

        ContinuousDistribution1D.__init__(self, name=name, pdf=pl_pdf, cdf=pl_cdf, inv_cdf=pl_inv_cdf)
        


# ------------ DEFAULT BOUNDED POWER LAW DISTRIBUTION (START) ------------
pl_default = PowerLawBounded()
# ------------ DEFAULT BOUNDED POWER LAW DISTRIBUTION (END) ------------






