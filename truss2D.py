import numpy as np


def stiffnessMatrix(E, A, L):
    S = np.array(np.zeros((6,6)))

    S[0][0] = E*A/L
    S[0][3] =-E*A/L
    S[3][0] =-E*A/L
    S[3][3] = E*A/L

    return S

def eqvLoad(FP):
    FEQV = np.array(np.zeros(6))

    FEQV[0] = FP
    FEQV[1] = 0
    FEQV[2] = 0
    FEQV[3] =-FP
    FEQV[4] = 0
    FEQV[5] = 0

    return FEQV
