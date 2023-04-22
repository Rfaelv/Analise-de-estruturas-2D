import numpy as np


def rotation2D(x_coord, y_coord, node_number, elem_number):
    R = np.array(np.zeros((6,6)))     
    dx = x_coord[node_number[1][elem_number]-1] - x_coord[node_number[0][elem_number]-1]
    dy = y_coord[node_number[1][elem_number]-1] - y_coord[node_number[0][elem_number]-1]
    L = (dx**2 + dy**2)**(1/2)
    cx = dx/L
    cy = dy/L

    R[0][0] = cx
    R[0][1] = cy
    R[1][0] = -cy
    R[1][1] = cx
    R[2][2] = 1
    R[3][3] = cx
    R[3][4] = cy
    R[4][3] = -cy
    R[4][4] = cx
    R[5][5] = 1

    return [R, L]

def gauss(NGL, SG, F):
    U = np.array(np.zeros(NGL))

    for k in range(NGL):
        for i in range(k, NGL):
            F[i] = F[i] - SG[i][k-1] / SG[k-1][k-1]*F[k-1]
            for j in range(k, NGL):
                SG[i][j] = SG[i][j] - SG[i][k-1] / SG[k-1][k-1] * SG[k-1][j]


    U[NGL-1] = F[NGL-1] / SG[NGL-1][NGL-1]
    for k in range(NGL-2, 0, -1):
        U[k] = F[k]
        for j in range(k+1, NGL):
            U[k] = U[k] - SG[k][j] * U[j]
        U[k] = U[k] / SG[k][k]
    
    return U

