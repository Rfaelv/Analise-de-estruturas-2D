import numpy as np
from model import Model
import frame2D
import truss2D
import matrix


class Anest:
    def __init__(self, model):
        self.model = model

    def createSystem(self):   

        for elem_number in range(self.model.NEL):
            R, L = matrix.rotation2D(self.model.X, self.model.Y, self.model.NO, elem_number)
            E, A, IZ, FP, QX, QY = self.model.getElementProperties(elem_number)

            if self.model.TE[elem_number]==1:
                S = frame2D.stiffnessMatrix(E, A, IZ, L)
                F = frame2D.eqvLoad(QX, QY, FP, L)
                
            elif self.model.TE[elem_number]==2:
                S = truss2D.stiffnessMatrix(E, A, L)
                F = truss2D.eqvLoad(FP)
            
            self.model.addElementStiffMatrixToGlobalStiffMatrix(R.T @ S @ R, elem_number)
            self.model.addElementEqvForcesToGlobalForceVector(R.T @ F, elem_number)

        self.model.addNodalForcesToForceVector()
        self.model.recordStiffnessMatrix()
        self.model.recordForceVector()
        self.model.applyBonds()

            # for i in range(6):
            #     for j in range(6):
            #         SR[i][j]=0
            #         for k in range(6):
            #             SR[i][j]=SR[i][j]+S[i][k]*R[k][j]

            # Multiplicação RSR=R'.SR
            # for i in range(6):
            #     for j in range(6):
            #         RSR[i][j]=0
            #         for k in range(6):
            #             RSR[i][j]=RSR[i][j]+R[k][i]*SR[k][j]

            # Asicionar matriz de rigidez do elemento na matriz global
            # for i in range(6):
            #     for j in range(6):
            #         Gi = self.model.G[elem_number][i]-1
            #         Gj = self.model.G[elem_number][j]-1
            #         self.model.SG[Gi][Gj] = self.model.SG[Gi][Gj] + RSR[i][j]

            # Forças equivalentes no referencial global
            # for i in range(6):
            #     FEQVG[i]=0
            #     for j in range(6):
            #         FEQVG[i]=FEQVG[i]+R[j][i]*F[j]

            # # Adicionar forças nodais equivalentes do elemento em {F}
            # for i in range(6): 
            #     Gi = self.model.G[elem_number][i]-1
            #     self.model.F[Gi]=self.model.F[Gi]+FEQVG[i]

        # Adicionar forças nodais no vetor de forças {F}
        # for i in range(self.model.NNO):
        #     for j in range(self.model.GLN):
        #         self.model.F[self.model.D[i][j]-1]= self.model.F[self.model.D[i][j]-1] + self.model.FA[i][j]


        def solveSystem(self):
            pass