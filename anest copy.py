import elements.frame2D as frame2D
import elements.truss2D as truss2D
import matrix
import numpy as np


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

    def solveSystem(self):
        self.model.U = np.linalg.solve(self.model.SG, self.model.F)
    
    def calcInternLoads(self):
        for elem_number in range(self.model.NEL):
            R, L = matrix.rotation2D(self.model.X, self.model.Y, self.model.NO, elem_number)
            E, A, IZ, FP, QX, QY = self.model.getElementProperties(elem_number)

            if self.model.TE[elem_number]==1:
                S = frame2D.stiffnessMatrix(E, A, IZ, L)
                F = frame2D.eqvLoad(QX, QY, FP, L)
                
            elif self.model.TE[elem_number]==2:
                S = truss2D.stiffnessMatrix(E, A, L)
                F = truss2D.eqvLoad(FP)

            self.model.ESF[:,elem_number] = S @ R @ self.model.elementDisplacementVector(elem_number) - F
    
    def calcReactions(self):
        self.model.RA = self.model.SG0 @ self.model.U - self.model.F0