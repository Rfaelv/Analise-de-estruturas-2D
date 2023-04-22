import elements.frame2D as frame2D
import elements.truss2D as truss2D
import matrix
import numpy as np


def createSystem(model):   
    for elem_number in range(model.NEL):
        R, L = matrix.rotation2D(model.X, model.Y, model.NO, elem_number)
        E, A, IZ, FP, QX, QY = model.getElementProperties(elem_number)

        if model.TE[elem_number]==1:
            S = frame2D.stiffnessMatrix(E, A, IZ, L)
            F = frame2D.eqvLoad(QX, QY, FP, L)
            
        elif model.TE[elem_number]==2:
            S = truss2D.stiffnessMatrix(E, A, L)
            F = truss2D.eqvLoad(FP)
        
        model.addElementStiffMatrixToGlobalStiffMatrix(R.T @ S @ R, elem_number)
        model.addElementEqvForcesToGlobalForceVector(R.T @ F, elem_number)

    model.addNodalForcesToForceVector()
    model.recordStiffnessMatrix()
    model.recordForceVector()
    model.applyBonds()

def solveSystem(model):
    model.U = np.linalg.solve(model.SG, model.F)

def calcInternLoads(model):
    for elem_number in range(model.NEL):
        R, L = matrix.rotation2D(model.X, model.Y, model.NO, elem_number)
        E, A, IZ, FP, QX, QY = model.getElementProperties(elem_number)

        if model.TE[elem_number]==1:
            S = frame2D.stiffnessMatrix(E, A, IZ, L)
            F = frame2D.eqvLoad(QX, QY, FP, L)
            
        elif model.TE[elem_number]==2:
            S = truss2D.stiffnessMatrix(E, A, L)
            F = truss2D.eqvLoad(FP)

        model.ESF[:,elem_number] = S @ R @ model.elementDisplacementVector(elem_number) - F

def calcReactions(model):
    model.RA = model.SG0 @ model.U - model.F0