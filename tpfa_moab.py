#!/usr/bin/env python3

    # TODO:
    # - Rever o uso do scd (para elementos de tamanho variável talvez seja melhor montar
    #   a malha como não-estruturada ou ler do arquivo)

import sys
import numpy as np
from pymoab import core, types, rng
from pymoab.scd import ScdInterface, ScdBox
from pymoab.hcoord import HomCoord
from pymoab.tag import Tag
from math import ceil

num_elem_x = 10
num_elem_y = 3
num_elem_z = 0
dx = 2.0
dy = 5/3
dz = 0
dim = 1

def get_centroid_coords(v):
    global dx, dy, dz
    centroid_x = v[0] + (dx/2)
    centroid_y = v[1] + (dy/2)
    centroid_z = v[2] + (dz/2)
    return np.array([centroid_x, centroid_y, centroid_z])

def main():
    global num_elem_x, num_elem_y, num_elem_z, dx, dy, dz, dim

    if len(sys.argv) == 3:
        num_elem_x = int(sys.argv[1])
        dx = float(sys.argv[2])
        dim = 0
    elif len(sys.argv) == 5:
        num_elem_x = int(sys.argv[1])
        num_elem_y = int(sys.argv[3])
        dx = float(sys.argv[2])
        dy = float(sys.argv[4])
        dim = 1
    elif len(sys.argv) == 7:
        num_elem_x = int(sys.argv[1])
        num_elem_y = int(sys.argv[3])
        num_elem_z = int(sys.argv[5])
        dx = float(sys.argv[2])
        dy = float(sys.argv[4])
        dz = float(sys.argv[6])
        dim = 2

    mbcore = core.Core()

    vertex_coords = np.array([])
    for k in range(num_elem_z + 1)
        for j in range(num_elem_y + 1):
            for i in range(num_elem_x + 1):
                vertex_coords = np.append(vertex_coords, [i*dx, j*dy, k*dz])

    scdint = ScdInterface(mbcore)
    scdbox = scdint.construct_box(HomCoord(0,0,0), \
                                  HomCoord(ceil(num_elem_x*dx), ceil(num_elem_y*dy), ceil(num_elem_z*dz)), \
                                  vertex_coords)

    # Incialização da matriz de permeabilidade e de conectividade
    K_perm = np.repeat(1, scdbox.num_elements())
    connectivity = np.zeros((scdbox.num_elements(), scdbox.num_elements()), dtype=np.bool_)

    # Pegar adjacências e montar matriz de conectividade
    entities = mbcore.get_entities_by_handle(0)
    adjacencies = [mbcore.get_adjacencies(e, dim, True) for e in entities if mbcore.type_from_handle(e) == types.MBHEX]
    i, j = 0, 0
    for a in adjacencies:
        for b in adjacencies:
            if b != a:
                intersection = rng.intersect(a, b)
                if not intersection.empty():
                    connectivity[i][j] = 1
                    connectivity[j][i] = 1
            j += 1
        i += 1

    # Determinando as coordenadas do centroide de cada elemento e aramazenando-as em tags.
    centroid_tag = mbcore.tag_get_handle('centroid', 3, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
    for e in entities:
        if mbcore.type_from_handle(e) == types.MBHEX:
            elem_vertex = mbcore.get_connectivity(e)
            centroid_coord = get_centroid_coords(mbcore.get_coords(elem_vertex[0]))
            mbcore.tag_set_data(centroid_tag, e, centroid_coord)

    # Montagem da matriz de coeficientes do sistema




if __name__ == '__main__':
    main()
