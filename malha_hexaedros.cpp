/*
    Compilar: make malha_hexaedros
    Executar: ./malha_hexaedros <número de hexaedros em x> <tamanho do hexaedro em x> \
                                <número de hexaedros em y> <tamanho do hexaedro em y> \
                                <número de hexaedros em z> <tamanho do hexaedro em z>
*/

#include <iostream>      /* std::cout */
#include <cstdlib>       /* atoi, atof */
#include <cmath>         /* ceil */
#include "moab/Core.hpp"    /* moab::Core, moab::ErrorCode */

using namespace std;
using namespace moab;

int main(int argc, char *argv[]) {
    ErrorCode rval;
    Core mbcore;
    int num_hex_x_axis, num_hex_y_axis, num_hex_z_axis;
    double dx, dy, dz;

    // Verificação da entrada e atribuição das variáveis.
    if (argc < 7) {
        num_hex_x_axis = 2;
        dx = 1;
        num_hex_y_axis = 2;
        dy = 1;
        num_hex_z_axis = 2;
        dz = 1;
    }
    else {
        num_hex_x_axis = atoi(argv[1]);
        dx = atof(argv[2]);
        num_hex_y_axis = atoi(argv[3]);
        dy = atof(argv[4]);
        num_hex_z_axis = atoi(argv[5]);
        dz = atof(argv[6]);
    }

    // Inicialização do vetor contendo as coordenadas dos vértices na malha.
    const int NUM_VERTEX = (num_hex_x_axis + 1)*(num_hex_y_axis + 1)*(num_hex_z_axis + 1);
    const int NUM_HEXAHEDRA = num_hex_x_axis*num_hex_y_axis*num_hex_z_axis;
    double vertex_coords[NUM_VERTEX*3] = {};
    int it = 0;

    for (int k = 0; k < (num_hex_z_axis + 1); k++) {
        for (int j = 0; j < (num_hex_y_axis + 1); j++) {
            for (int i = 0; i < (num_hex_x_axis + 1); i++) {
                vertex_coords[it] = i*dx;
                vertex_coords[it+1] = j*dy;
                vertex_coords[it+2] = k*dz;
                it += 3;
            }
        }
    }

    // Criando os vértices da malha
    Range vertex_handles;
    rval = mbcore.create_vertices(vertex_coords, NUM_VERTEX, vertex_handles);
    MB_CHK_SET_ERR(rval, "create_vertices failed");

    // A matriz connectivity aramazena em cada linha os 8 vértices que compõem o hexaedro.
    // Os oito vértices são determinados escolhendo um intervalo de cada coordenada (x, y, z)
    // e verificando quais vértices estão contidos na região do intervalo. Se um vértice pertence
    // ao intervalo, então ele faz parte do hexaedro.
    double x_inf = 0.0, x_sup = dx, y_inf = 0.0, y_sup = dy, z_inf = 0.0, z_sup = dz;
    int m = 0, n = 0;
    EntityHandle connectivity[NUM_HEXAHEDRA][8];

    while (x_sup <= num_hex_x_axis*dx) {
        while (y_sup <= num_hex_y_axis*dy) {
            while (z_sup <= num_hex_z_axis*dz) {
                n = 0;
                for (int i = 0; i < vertex_handles.size(); i++) {
                    double x_i = vertex_coords[3*i], y_i = vertex_coords[3*i+1], z_i = vertex_coords[3*i+2];
                    if ((x_i >= x_inf && x_i <= x_sup) &&
                        (y_i >= y_inf && y_i <= y_sup) &&
                        (z_i >= z_inf && z_i <= z_sup)) {
                            connectivity[m][n] = vertex_handles[i];
                            n += 1;
                    }
                }
                z_inf = z_sup;
                z_sup += dz;
                if (n == 8) m += 1;
            }
            y_inf = y_sup;
            y_sup += dy;
            z_inf = 0.0;
            z_sup = dz;
        }
        x_inf = x_sup;
        x_sup += dx;
        y_inf = 0.0;
        y_sup = dy;
    }

    // Criando os elementos da malha, i.e., os hexaedros.
    Range hex_handles;
    EntityHandle elem;

    for (int i = 0; i < NUM_HEXAHEDRA; i++) {
        rval = mbcore.create_element(MBHEX, connectivity[i], 8, elem);
        MB_CHK_SET_ERR(rval, "create_element failed");
        hex_handles.insert(elem);
    }
    cout << "Entities created: " << hex_handles;


    // Escrevendo os elementos da malha num arquivo.
    rval = mbcore.write_file("my_first_mesh.vtk");
    MB_CHK_SET_ERR(rval, "write_file failed");

    return 0;
}
