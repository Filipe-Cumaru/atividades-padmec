#!/usr/bin/env python3
import sys
import numpy as np
from pymoab import core, types, rng
from pymoab.tag import Tag
from math import ceil

# Variáveis globais
# ni -> Número de elementos ao longo do eixo i
# di -> Tamanho do elemento ao longo do eixo i
# dim -> Número de dimensões das adjacências dos elementos (ex.: numa malha 2D,
#        os elementos compartilham arestas que tem dim = 1)
nx = 10
ny = 3
nz = 1
dx = 2.0
dy = 5/3
dz = 1
dim = 2
num_elements = 30

def next_index(i):
    global nx, ny, nz

    if ((i+1) % (nx+1)) == 0:
        i += 1
        if (i+1) % ((nx+1)*ny+1) == 0:
            i += (nx+1)

    return i

# create_mesh_connectivity: Função para montagem da conectividade das elementos.
# Parâmetros:
#   - vertex_handles: vetor de EntityHandles representando vértices do MOAB.
#   - vertex_coords: coordenadas dos vértices
def create_mesh_connectivity(vertex_handles, vertex_coords):
    global nx, ny, nz, num_elements

    k = 0
    mesh_connectivity = np.zeros((num_elements, 8), dtype=np.uint64)

    for i in range(num_elements):
        mesh_connectivity[i] = [vertex_handles[k], vertex_handles[k+1],\
                                vertex_handles[k+nx+1], vertex_handles[k+nx+2],\
                                vertex_handles[k+(nx+1)*(ny+1)], vertex_handles[k+(nx+1)*(ny+1)+1],\
                                vertex_handles[k+(nx+1)*(ny+2)], vertex_handles[k+(nx+1)*(ny+2)+1]]
        k = next_index(k+1)

    return mesh_connectivity

def get_centroid_coords(v):
    global dx, dy, dz
    centroid_x = v[0] + (dx/2)
    centroid_y = v[1] + (dy/2)
    centroid_z = v[2] + (dz/2)
    return np.array([centroid_x, centroid_y, centroid_z])

def equiv_perm(k1, k2):
    return (2*k1*k2)/(k1 + k2)

def centroid_dist(c1, c2):
    return (c1[0] + c2[0])**2 + (c1[1] + c2[1])**2 + (c1[2] + c2[2])**2


def main():
    global nx, ny, nz, dx, dy, dz, dim, num_elements

    # Tratamento da entrada. O número de dimensões da malha é determinado a partir
    # da quantidade de argumentos.
    if len(sys.argv) == 7:
        nx = int(sys.argv[1])
        ny = int(sys.argv[3])
        nz = int(sys.argv[5])
        dx = float(sys.argv[2])
        dy = float(sys.argv[4])
        dz = float(sys.argv[6])
        dim = 2
        num_elements = nx*ny*nz
        num_vertex = (nx+1)*(ny+1)*(nz+1)
    else:
        print("Not enough arguments")
        return

    # Criando instância da classe Core que gerencias as operações na malha.
    mbcore = core.Core()

    # Inicializando o vetor de coordenadas dos vértices.
    vertex_coords = np.zeros(num_vertex*3)
    for i in range(num_vertex):
        vertex_coords[3*i] = (i % (nx+1))*dx
        vertex_coords[3*i+1] = ((i // (nx+1)) % (ny+1))*dy
        vertex_coords[3*i+2] = ((i // ((nx+1)*(ny+1))) % (nz+1))*dz

    # O método create_vertices cria os handles associados a cada coordenada em vertex_coords
    vertex_handles = mbcore.create_vertices(vertex_coords)

    # Em mesh_connectivity são aramazenados os conjuntos de vértices que compõem um elemento,
    # ou seja, determina a conectividade dos vértices na malha.
    print("Creating connectivity")
    mesh_connectivity = create_mesh_connectivity(vertex_handles, vertex_coords)
    print("Done")
    for m in mesh_connectivity:
        print(m)
    return

    # De posse da conectividade da malha, criamos os elementos um a um. A troca de valores
    # nas duas primeiras linhas do laço são necessárias devido a forma como o MOAB interpreta
    # a malha e as adjacências dos elementos.
    elem_handles = rng.Range()
    for c in mesh_connectivity:
        c[2], c[3] = c[3], c[2]
        c[6], c[7] = c[7], c[6]
        temp = mbcore.create_element(types.MBHEX, c)
        elem_handles.insert(temp)

    # Incialização da matriz de conectividade. (Neste caso, a conectividade é em
    # relação aos elementos, ou seja, quais elementos são vizinhos.)
    connectivity = np.zeros((num_elements, num_elements), dtype=np.bool_)

    # Encontrando adjacências para preencher a matriz de conectividade
    adjacencies = [mbcore.get_adjacencies(e, dim, True) for e in elem_handles]

    # Para cada adjacência diferente, verifica-se se existem uma fronteira compartilhada.
    # Caso positivo, os dois elementos são vizinhos e isto é indicado em connectivity.
    i, j = 0, 0
    for a in adjacencies:
        for b in adjacencies:
            if b != a:
                intersection = rng.intersect(a, b)
                if not intersection.empty():
                    connectivity[i][j] = 1
                    connectivity[j][i] = 1
            j += 1
        j = 0
        i += 1

    # Determinando as coordenadas do centroide de cada elemento e aramazenando-as em tags.
    # Uma tag é um valor associado a cada elemento. Aqui, cada elemento possui duas tags: uma
    # que armazena o valor das coordenadas do centroide e outra que armazena a permeabilidade.
    centroid_tag = mbcore.tag_get_handle('centroid', 3, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
    permeability_tag = mbcore.tag_get_handle('permeability', 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
    for e in elem_handles:
        elem_vertex = mbcore.get_connectivity(e)
        centroid_coord = get_centroid_coords(mbcore.get_coords([elem_vertex[0]]))
        mbcore.tag_set_data(centroid_tag, e, centroid_coord)
        mbcore.tag_set_data(permeability_tag, e, np.array([1], dtype=np.float_))

    mbcore.write_file("tpfa_mesh.h5m")
    print("New h5m file created")


if __name__ == '__main__':
    main()
