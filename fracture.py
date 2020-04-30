import numpy as np
import matlab_clones as mc


class Fracture:
    """
    Classe per la descrizione delle fratture: riceve i vari parametri in ingresso ed effettua le dovute rotazioni
    """

    def __init__(self, n, rx, alpha, ar, v, bar):
        """
        Costruttore della classe
        :param n: numero di vertici
        :param rx: lunghezza semiasse ascisse
        :param alpha: angolo di rotazione attorno all'asse z (prima rotazione)
        :param ar: aspect ratio: rapporto tra semiasse x e semiasse y dell'ellisse circoscritta al poligono
        :param v: versore normale
        :param bar: baricentro
        """

        self.n = n
        self.rx = rx
        self.ry = rx / ar
        self.alpha = alpha
        self.phi = mc.cart2sph(v[0], v[1], v[2])[1]
        self.teta = mc.cart2sph(v[0], v[1], v[2])[0]
        self.vn = v
        self.bar = bar
        # vertici = matrice 3xn contenente i vertici come vettore colonna ordinati in senso antiorario
        self.vertici = np.zeros((3, n))

        # Si generano i vertici del poligono
        for i in range(n):
            self.vertici[0, i] = self.rx * np.cos(i * 2 * np.pi / n)
            self.vertici[1, i] = self.ry * np.sin(i * 2 * np.pi / n)

        # Si ruotano i vertici
        rotMatrix = np.dot(np.dot(mc.rotz(self.teta), mc.roty(np.pi / 2 - self.phi)), mc.rotz(alpha))
        self.vertici = np.dot(rotMatrix, self.vertici)

        # Si traslano i vertici
        for i in range(n):
            self.vertici[:, i] = self.vertici[:, i] + bar

        self.xmin = np.min(self.vertici[0, :])
        self.xmax = np.max(self.vertici[0, :])
        self.ymin = np.min(self.vertici[1, :])
        self.ymax = np.max(self.vertici[1, :])
        self.zmin = np.min(self.vertici[2, :])
        self.zmax = np.max(self.vertici[2, :])