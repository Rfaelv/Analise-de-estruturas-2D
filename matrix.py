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

def gauss():
    pass