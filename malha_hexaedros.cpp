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
#include "moab/ScdInterface.hpp"    /* moab::ScdInterface, moab::ScdBox */

using namespace std;
using namespace moab;

int main(int argc, char *argv[]) {
    ErrorCode ret;
    Core mbcore;
    int num_hex_x_axis, num_hex_y_axis, num_hex_z_axis;
    double dx, dy, dz;

    // Verificação da entrada e atribuição das variáveis.
    if (argc < 7) {
        num_hex_x_axis = 2;
        dx = 0.1;
        num_hex_y_axis = 2;
        dy = 0.1;
        num_hex_z_axis = 2;
        dz = 0.1;
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

    // A malha é criada usando a classe ScdInterface para malhas estruturadas,
    // já que a malha é constituída por hexaedros regulares.
    ScdInterface *structured_interface;
    ret = mbcore.query_interface(structured_interface);
    MB_CHK_SET_ERR(ret, "query_interface failed");

    ScdBox *structured_box = NULL;
    ret = structured_interface->construct_box(HomCoord(0.0, 0.0, 0.0),
            HomCoord(ceil(num_hex_x_axis*dx), ceil(num_hex_y_axis*dy), ceil(num_hex_z_axis*dz)),
            vertex_coords,
            NUM_VERTEX,
            structured_box);
    MB_CHK_SET_ERR(ret, "construct_box failed");

    // Iterando sobre os elementos da malha.
    for (int i = 0; i < (num_hex_x_axis + 1); i++) {
        for (int j = 0; j < (num_hex_y_axis + 1); j++) {
            for (int k = 0; k < (num_hex_z_axis + 1); k++) {
                cout << "Vertex (" << i*dx << "," << j*dy << "," << k*dz << ") "
            		  << "has handle: " << structured_box->get_vertex(i,j,k) << endl;
            }
        }
    }

    return 0;
}
