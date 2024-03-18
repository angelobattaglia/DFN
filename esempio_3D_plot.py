"""
Questo modulo contiente degli esempi per la visualizzazione 3D di poligoni tramite il pacchetto plotly.
Per approfondimenti visitare: https://plot.ly/python/
"""

import numpy as np
import plotly.graph_objs as go
import plotly.offline as poff
import DFN

ottagono_3d = np.array([[-0.33095549,  2.80413031,  5.46457847,  6.63418824,  5.86620843,
                         3.45398123,  0.31889543, -2.34155273, -3.5111625 , -2.74318269],
                        [-0.24405635,  0.3613702 ,  3.7662977 ,  8.67015958, 13.19984728,
                         15.62517404, 15.01974749, 11.61481999,  6.71095811,  2.18127042],
                        [10.55212849,  9.10664097,  6.003798  ,  2.42878012, -0.25287733,
                         -1.01687236,  0.42861517,  3.53145814,  7.10647601,  9.78813346]])

x = ottagono_3d[0, :].tolist()
y = ottagono_3d[1, :].tolist()
z = ottagono_3d[2, :].tolist()

# Ripetizione primo punto per "chiudere" il giro intorno al poligono
x.append(x[0])
y.append(y[0])
z.append(z[0])

perimetro = go.Scatter3d(x=x, y=y, z=z,
                         mode='lines',
                         marker=dict(
                             color='red'
                         )
                         )
area = go.Mesh3d(x=x, y=y, z=z, color='#FFB6C1', opacity=0.60)

fig_3d = go.Figure(data=[perimetro, area])
poff.plot(fig_3d)

# ATTENZIONE: Se il poligono risultasse parallelo al piano (y,z), cioè avesse x costante per ogni suo vertice,
# si deve aggiungere il seguente attributo: area.delaunayaxis = 'x'. Analogamente funziona per y e z costante.
# Se non aggiunto, si ottiene l'effetto seguente:

x_costante = [0] * len(x)

perimetro_cost = go.Scatter3d(x=x_costante, y=y, z=z,
                              mode='lines',
                              marker=dict(
                                  color='red'
                              )
                              )
area_cost = go.Mesh3d(x=x_costante, y=y, z=z, color='#FFB6C1', opacity=0.60)

fig_3d_cost = go.Figure(data=[perimetro_cost, area_cost])
poff.plot(fig_3d_cost)


# Per realizzare più plot nella stessa figura, basta aggiungerli alla lista passanta in argomento come "data"
# all'oggetto go.Figure, cioè:

all_polygons = [perimetro, area, perimetro_cost, area_cost]
fig_3d_alltogether = go.Figure(data=all_polygons)
poff.plot(fig_3d_alltogether)




