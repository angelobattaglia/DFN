import distributions as dist
import fracture as frac
import numpy as np
import plotly.graph_objs as go
import plotly.offline as poff
import dill
import matlab_clones as mc
import trace_1


class DiscreteFractureNetwork:
    """
    Classe per la descrizione del DiscreteFractureNetwork
    """

    def __init__(self, N, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax,
                 alpha_pl, radius_l, radius_u, k, mode_vector, tol=1.0e-10, fixed_n_edges=0):
        """
        Costruttore della classe
        :param N: numero di fratture
        :param Xmin, Xmax, Ymin, Ymax, Zmin, Zmax: estremi del dominio del DFN
        :param alpha_pl: parametro (> 1) della distribuzione Power LawBounded per i semiassi delle ascisse
        :param radius_l: limite inferiore dei semiassi delle ascisse
        :param radius_u: limite superiore dei semiassi delle ascisse
        :param k: parametro di concentrazione della legge Von Mises-Fisher
        :param mode_vector: punto medio della legge Von Mises-Fisher
        :param fixed_n_edges: numero fissato dei lati dei poligoni (>= 8); se non specificato ogni poligono ha
                              un numero casuale di lati da 8 a 16 (distribuiti uniformemente)
        :param tol: tolleranza
        """

        self.N = 0  # il valore verrà poi aggiornato nel metodo genfrac
        self.tol = tol
        self.Xmin = Xmin
        self.Xmax = Xmax
        self.Ymin = Ymin
        self.Ymax = Ymax
        self.Zmin = Zmin
        self.Zmax = Zmax

        self.fixed_n_edges = fixed_n_edges

        self.alpha_pl = alpha_pl
        self.radius_l = radius_l
        self.radius_u = radius_u
        # Si definisce la distribuzione power law bounded per i semiassi delle x
        self.pl_dist = dist.PowerLawBounded(alpha=alpha_pl, radius_l=radius_l, radius_u=radius_u)

        self.k = k
        self.mode_vector = mode_vector / np.linalg.norm(mode_vector) # nel caso il modulo non fosse unitario
        # Si definisce la distribuzione Von Mises-Fisher
        self.vmf_dist = dist.VonMisesFisher(k=k, mode_vector=self.mode_vector)

        self.fractures = []             # Lista delle fratture
        self.poss_intersezioni = []     # Lista di liste con le possibili intersezioni
        self.intersezioni = []          # Lista di liste con le effettive intersezioni
        self.frac_traces = []           # Lista di liste con le tracce per ogni poligono
        self.traces = []                # Lista con tutte le tracce del DFN

        self.genfrac(N)

    def genfrac(self, n_to_gen):
        """
        Genera n_to_gen fratture, aggiornando opportunamente tutte le strutture dati del DFN
        :param n_to_gen: numero di fratture da generare
        """

        semiaxes_x = self.pl_dist.sample(n_to_gen)  # Vettore di n_to_gen semiassi x
        ars = np.random.uniform(1, 3, n_to_gen)     # Vettore di n_to_gen aspect ratio
        normals = self.vmf_dist.sample(n_to_gen)    # Matrice 3 x n_to_gen contenente le normali

        if self.fixed_n_edges == 0:
            n_edges = np.random.random_integers(8, 16, n_to_gen)    # Vettore di n_to_gen numero di lati
        else:
            n_edges = [self.fixed_n_edges] * n_to_gen

        alpha_angles = np.random.uniform(0, 2 * np.pi, n_to_gen) # Vettore di n_to_gen angoli di rotazione
                                                                 # attorno alla normale

        # Matrice n_to_gen x 3 contenente i baricentri
        centers = np.random.uniform(np.array([self.Xmin, self.Ymin, self.Zmin]),
                                    np.array([self.Xmax, self.Ymax, self.Zmax]),
                                    (n_to_gen, 3))

        for i in range(n_to_gen):
            tmp = frac.Fracture(n_edges[i], semiaxes_x[i], alpha_angles[i], ars[i], normals[:, i], centers[i, :])
            self.fractures.append(tmp)
            self.poss_intersezioni.append([])
            self.intersezioni.append([])
            self.frac_traces.append([])

        if self.N == 0: # Si entra in questo blocco solo al momento della creazione del DFN
            self.N += n_to_gen
            for i in range(self.N - 1):
                for j in range(i + 1, self.N):
                    self.aggiorna_int(i, j)

        else: # Si entra in questo blocco nel caso di aggiunte di poligoni al DFN già esistente
            # Confrontiamo le vecchie fratture con le nuove
            for i in range(self.N):
                for j in range(n_to_gen):
                    self.aggiorna_int(i, self.N + j)

            # Confrontiamo le nuove fratture tra loro
            for i in range(self.N, self.N + n_to_gen - 1):
                for j in range(i + 1, self.N + n_to_gen):
                    self.aggiorna_int(i, j)

            self.N += n_to_gen

    def aggiungi(self, v):
        """
        Aggiunge i poligoni, già esistenti, presenti nella lista v
        :param v: lista di oggetti di classe Fracture
        """
        for i in range(len(v)):
            self.fractures.append(v[i])
            self.poss_intersezioni.append([])
            self.intersezioni.append([])
            self.frac_traces.append([])

        # Confrontiamo le vecchie fratture con le nuove
        for i in range(self.N):
            for j in range(len(v)):
                self.aggiorna_int(i, self.N + j)

        # Confrontiamo le nuove fratture tra loro
        for i in range(self.N, self.N + len(v) - 1):
            for j in range(i + 1, self.N + len(v)):
                self.aggiorna_int(i, j)

        self.N += len(v)

    def rimuovi(self, v):
        """
        Rimuove le fratture nelle posizioni indicate dalla lista v, e aggiorna opportunamente le strutture dati del DFN
        :param v: lista contenente gli indici dei poligoni da rimuovere
        """
        v = list(set(v))
        v = sorted(v, reverse=True) # Per evitare errori che possono sorgere con la rinumerazione
        for i in v:
            self.fractures.pop(i)
            self.poss_intersezioni.pop(i)
            self.intersezioni.pop(i)
            for j in range(len(self.frac_traces[i])):
                # Presa la j-esima traccia generata dall'i-esimo poligono, la si rimuove dalla lista delle tracce
                # generate dall'altro genitore
                tr = self.frac_traces[i][j]
                if tr.i1 == i:
                    self.frac_traces[tr.i2].remove(tr)
                else:
                    self.frac_traces[tr.i1].remove(tr)
                self.traces.remove(tr)
            self.frac_traces.pop(i)

            for j in range(len(self.poss_intersezioni)):
                if i in self.poss_intersezioni[j]:
                    if i in self.intersezioni[j]:
                        self.intersezioni[j].remove(i)

                    self.poss_intersezioni[j].remove(i)

                # Si rinumerano gli indici relativi ai poligoni
                for k in range(len(self.poss_intersezioni[j])):
                    if self.poss_intersezioni[j][k] > i:
                        self.poss_intersezioni[j][k] -= 1
                for k in range(len(self.intersezioni[j])):
                    if self.intersezioni[j][k] > i:
                        self.intersezioni[j][k] -= 1
            # Si rinumerano gli indici dei genitori delle tracce
            for tr in self.traces:
                if tr.i1 > i:
                    tr.i1 -= 1
                if tr.i2 > i:
                    tr.i2 -= 1

        self.N = self.N - len(v)

    def inters_BB(self, fr1, fr2):
        """
        Metodo per controllare se i bounding box di due fratture si intersecano
        :param fr1: oggetto della classe Fracture
        :param fr2: oggetto della classe Fracture
        :return: booleano
        """
        max_x_min = max(fr1.xmin, fr2.xmin)
        min_x_max = min(fr1.xmax, fr2.xmax)
        max_y_min = max(fr1.ymin, fr2.ymin)
        min_y_max = min(fr1.ymax, fr2.ymax)
        max_z_min = max(fr1.zmin, fr2.zmin)
        min_z_max = min(fr1.zmax, fr2.zmax)
        if max_x_min <= min_x_max and max_y_min <= min_y_max and max_z_min <= min_z_max:
            return True
        else:
            return False

    def gen_trace(self, i1, i2):
        """
        Metodo per calcolare la traccia tra due poligoni
        :param i1, i2: indici dei poligoni
        :return: oggetto della classe Trace (se la traccia esiste); None se la traccia non esiste
        """

        p1 = self.fractures[i1]
        p2 = self.fractures[i2]

        # Calcoliamo la retta d'intersezione tra i piani contenenti i poligoni
        t = np.cross(p1.vn, p2.vn)
        A = np.array([p1.vn, p2.vn, t])
        b = np.array([np.dot(p1.vertici[:, 0], p1.vn), np.dot(p2.vertici[:, 0], p2.vn), 0])

        r0 = np.linalg.solve(A, b)
        s = []
        s.extend(self.inters_2D(p1, t, r0))
        if len(s) == 2:
            s.extend(self.inters_2D(p2, t, r0))
            if len(s) == 4:
                q = np.argsort(s)
                if (q[0] == 0 and q[1] == 1) or (q[0] == 2 and q[1] == 3):
                    return None
                else:
                    s.sort()
                    x1 = r0 + s[1] * t
                    x2 = r0 + s[2] * t
                    tr = trace_1.Trace(p1, p2, i1, i2, np.array([x1, x2]).T)
                    return tr

        else:
            return None

    def inters_2D(self, p, t, r0):
        """
        Metodo che, data una retta in forma parametrica (X(s) = r0 + s*t) e una frattura, determina se si intersecano
        in un segemento e ne calcola gli estremi
        :param p: oggetto della classe Fracture
        :param t: direzione della retta
        :param r0: punto della retta
        :return: lista s con i due nuemeri reali [s1, s2] che parametrizzano il segmento se c'e' intersezione;
                 lista vuota se non c'e' intersezione
        """

        # Si ruota il poligono per ricondursi sul piano xy
        rotMatrix = np.dot(np.dot(mc.rotz(- p.alpha), mc.roty(p.phi - np.pi / 2)), mc.rotz(- p.teta))
        vertici = p.vertici.copy()

        for i in range(p.n):
            vertici[:, i] -= p.bar
        vertici = np.dot(rotMatrix, vertici)
        t = np.dot(rotMatrix, t)
        r0 = np.dot(rotMatrix, r0 - p.bar)

        conta = 0   # Contatore dei lati intersecati
        i = 0
        s = []
        while conta < 2 and i < p.n:
            j = (i + 1) % p.n
            A = np.array([[t[0], vertici[0, i] - vertici[0, j]], [t[1], vertici[1, i] - vertici[1, j]]])
            b = np.array([vertici[0, i] - r0[0], vertici[1, i] - r0[1]])
            if - self.tol <= np.linalg.det(A) <= self.tol: # Si controlla il parallelismo
                s1 = (vertici[0, i] - r0[0]) / t[0]
                s2 = (vertici[1, i] - r0[1]) / t[1]
                if abs(s1 - s2) < self.tol:
                    s.append(s1)
                    s2 = (vertici[0, j] - r0[0]) / t[0]
                    s.append(s2)
                    return s
                i += 1
            else:
                x = np.linalg.solve(A, b)
                if - self.tol <= x[1] <= 1 + self.tol:
                    conta += 1
                    s.append(x[0])
                i += 1

        s.sort()
        return s

    def aggiorna_int(self, i, j):

        fr1 = self.fractures[i]
        fr2 = self.fractures[j]
        if self.inters_BB(fr1, fr2) is True:
            self.poss_intersezioni[i].append(j)
            self.poss_intersezioni[j].append(i)
            tr = self.gen_trace(i, j)
            if tr is not None:
                self.traces.append(tr)
                self.frac_traces[i].append(tr)
                self.frac_traces[j].append(tr)
                self.intersezioni[i].append(j)
                self.intersezioni[j].append(i)

    def visual3D(self, filename='tmp-plot'):
        """
        Metodo per la visualizzazione grafica delle fratture e delle tracce
        :param filename: nome del file che verrà creato
        """

        all_polygons = []
        for i in range(self.N):
            x = self.fractures[i].vertici[0, :].tolist()
            y = self.fractures[i].vertici[1, :].tolist()
            z = self.fractures[i].vertici[2, :].tolist()
            x.append(x[0])
            y.append(y[0])
            z.append(z[0])
            perimetro = go.Scatter3d(x=x, y=y, z=z, mode='lines', marker=dict(color='red'))
            area = go.Mesh3d(x=x, y=y, z=z, color='#FFB6C1', opacity=0.60)

            # Caso in cui il poligono sia parallelo ad uno dei piani coordinati
            if np.linalg.norm(self.fractures[i].vn - np.array([1, 0, 0])) < self.tol:
                area.delaunayaxis = 'x'
            elif np.linalg.norm(self.fractures[i].vn - np.array([0, 1, 0])) < self.tol:
                area.delaunayaxis = 'y'
            elif np.linalg.norm(self.fractures[i].vn - np.array([0, 0, 1])) < self.tol:
                area.delaunayaxis = 'z'

            all_polygons.append(perimetro)
            all_polygons.append(area)

        for i in range(len(self.traces)):
            x = self.traces[i].estremi[0, :].tolist()
            y = self.traces[i].estremi[1, :].tolist()
            z = self.traces[i].estremi[2, :].tolist()
            segmento = go.Scatter3d(x=x, y=y, z=z, mode='lines', marker=dict(color='black'))
            all_polygons.append(segmento)

        fig_3d_alltogether = go.Figure(data=all_polygons)
        poff.plot(fig_3d_alltogether, filename=filename+'.html')

    def scrittura1(self, filename='file1.txt'):
        """
        Metodo per scrivere su file come richiesto al punto 7
        """

        with open(filename, 'w') as f1:
            print(self.N, file=f1)
            for i in range(self.N):
                print(i, self.fractures[i].n, file=f1)
                for j in range(self.fractures[i].n):
                    print(self.fractures[i].vertici[0, j],
                          self.fractures[i].vertici[1, j],
                          self.fractures[i].vertici[2, j], file=f1)

    def scrittura2(self, filename='file2.txt'):
        """
        Metodo per scrivere su file come richiesto al punto 8
        """

        with open(filename, 'w') as f2:
            print(self.N, file=f2)

            somma_vertici = 0
            for i in range(self.N):
                somma_vertici = somma_vertici + self.fractures[i].n
            print(somma_vertici, file=f2)

            indice = 0
            for i in range(self.N):
                print(i, self.fractures[i].n, indice, file=f2)
                indice = indice + self.fractures[i].n

            for i in range(self.N):
                for j in range(self.fractures[i].n):
                    print(self.fractures[i].vertici[0, j],
                          self.fractures[i].vertici[1, j],
                          self.fractures[i].vertici[2, j], file=f2)

    def save(self):
        """
        Salva l'oggetto DiscreteFractureNetwork come file .pkl
        """

        with open('DFN.pkl', 'wb') as f3:
            dill.dump(self, f3)
