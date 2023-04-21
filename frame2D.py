import numpy as np


def stiffnessMatrix(E, A, IZ, L):
    S = np.array(np.zeros((6,6)))

    S[0][0] = E*A/L
    S[0][3] =-E*A/L
    S[1][1] = 12*E*IZ/L**3
    S[1][2] =  6*E*IZ/L**2
    S[1][4] =-12*E*IZ/L**3
    S[1][5] =  6*E*IZ/L**2
    S[2][2] =  4*E*IZ/L
    S[2][4] = -6*E*IZ/L**2
    S[2][5] =  2*E*IZ/L
    S[3][3] = E*A/L
    S[4][4] = 12*E*IZ/L**3
    S[4][5] = -6*E*IZ/L**2
    S[5][5] =  4*E*IZ/L
    
    for i in range(1,6):
        for j in range(0,i):
            S[i][j]=S[j][i]

    return S

def eqvLoad(QX, QY, FP, L):
    FEQV = np.array(np.zeros(6))

    FEQV[0] = QX*L/2     + FP
    FEQV[1] = QY*L/2
    FEQV[2] = QY*L**2/12
    FEQV[3] = QX*L/2     - FP
    FEQV[4] = QY*L/2
    FEQV[5] =-QY*L**2/12

    return FEQV
