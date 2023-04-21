import json
import numpy as np


class Model:
    def __init__(self, path):

        with open(path, 'r') as data:
            data  = json.load(data)

        self.GLN = data['nDegreeOfFreedom']        # NÚMERO DE GRAUS DE LIBERDADE POR NÓ     
        self.NELP = len(data['frameElements'])     # NÚMERO DE ELEMENTOS DE PÓRTICO
        self.NELT = len(data['trussElements'])     # NÚMERO DE ELEMENTOS DE TRELIÇA
        self.NEL = self.NELP + self.NELT           # NÚMERO DE ELEMENTOS
        self.NNO = len(data['nodes'])              # NÚMERO DE NÓS
        self.NGL = self.NNO * self.GLN             # NÚMERO DE GRAUS DE LIBERDADE
        self.NNV = len(data['fixedNodes'])         # NÚMERO DE NÓS VINCULADOS
        self.NLL = None                            # NÚMERO DE LIBERACOES LOCAIS
        self.NTC = len(data['elementProperties'])  # NÚMERO DE TIPOS DE CARACTERÍSTICAS DE ELEMENTO
        self.SG = np.array(np.zeros((self.NGL, self.NGL)))  # MATRIZ DE RIGIDEZ GLOBAL
        self.F = np.array(np.zeros((self.NGL)))             # VETOR DE FORÇAS GLOBAIS
        self.U = np.array(np.zeros(self.NGL))               # VETOR DE DESLOCAMENTOS
        self.RA = np.array(np.zeros((self.NGL)))            # VETOR DE REAÇÕES DE APOIO = [SG].{U} - {F}
        self.ESF = np.array(np.zeros((6, self.NEL)))        # ESFORÇO DIREÇÃO I ELEMENTO J  
        
        self.FA = np.array(np.zeros((self.NNO, self.NGL)))  # FORÇA APLICADA (CARGA NODAL) NO NÓ I NA DIREÇÃO J
        for nodalForce in data['nodalForces']:
            for i in range(self.GLN):
                self.FA[nodalForce[0]-1][i] = nodalForce[i+1]

        self.NO = np.array(np.zeros((2, self.NEL)), dtype=np.int)  # NÓ I DO ELEMENTO J
        self.QX = np.array(np.zeros(self.NEL))  # CARGA DISTRIBUÍDA DIR. X NO ELEMENTO I
        self.QY = np.array(np.zeros(self.NEL))  # CARGA DISTRIBUÍDA DIR. Y NO ELEMENTO I
        self.FP = np.array(np.zeros(self.NEL))  # FORCA DE PROTENSÇÃO NO ELEMENTO I
        self.NC = np.array(np.zeros(self.NEL), dtype=np.int)  # NúMERO DA CARACTERÍSTICA DO ELEMENTO I
        self.TE = np.array(np.zeros(self.NEL), dtype=np.int)  # TIPO DO ELEMENTO I (PÓRTICO=1, TRELIÇA=2)
        for element in data['frameElements']:
            for i in range(2):
                self.NO[i, element[0]-1] = int(element[i+1])
            self.NC[element[0]-1] = element[3]
            self.FP[element[0]-1] = element[4]
            self.QX[element[0]-1] = element[5]
            self.QY[element[0]-1] = element[6]
            self.TE[element[0]-1] = 1
        
        for element in data['trussElements']:
            for i in range(2):
                self.NO[i, element[0]-1] = int(element[i+1])
            self.NC[element[0]-1] = element[3]
            self.FP[element[0]-1] = element[4]
            self.TE[element[0]-1] = 2

        self.X = np.array(np.zeros(self.NNO))  # COORDENADA X DO NÓ I
        self.Y = np.array(np.zeros(self.NNO))  # COORDENADA Y DO NÓ I
        for i, node in enumerate(data['nodes']):
            self.X[i] = node[0]
            self.Y[i] = node[1]
        
        self.E = np.array(np.zeros(self.NTC))   # MÓDULO DE ELASTICIDADE DA CARACTERÍSTICA I
        self.A = np.array(np.zeros(self.NTC))   # ÁREA DA SEÇÃO DA CARACTERÍSTICA I
        self.IZ = np.array(np.zeros(self.NTC))  # INÉRCIA DA CARACTERÍSTICA I
        for i, elemPropertie in enumerate(data['elementProperties']):
            self.E[i] = elemPropertie[0]
            self.A[i] = elemPropertie[1]
            self.IZ[i] = elemPropertie[2]
        
        self.NV = np.array(np.zeros(len(data['fixedNodes'])), dtype=np.int)  # NÚMERO DO I-ÉSIMO NÓ VINCULADO
        self.V = np.array(np.zeros((self.NNO, self.GLN)), dtype=np.int)  # VÍNCULO DO NÓ I NA DIREÇÃO J (0=LIVRE, 1=IMPEDIDO)
        for i, fixedNode in enumerate(data['fixedNodes']):
            self.NV[i] = fixedNode[0]
            for j in range(self.GLN):
                self.V[i][j] = fixedNode[j+1]

        self.D = np.array(np.zeros((self.NNO, self.GLN)), dtype=np.int)  # DIREÇÃO DO GRAU DE LIBERDADE DO NÓ I NA DIREÇÃO J
        for i in range(self.NNO):
            for j in range(self.GLN):
                self.D[i][j]=self.GLN*(i)+j+1 
        
        self.G = np.array(np.zeros((self.NEL, 6)), dtype=np.int)  # G.L. GLOBAL DO ELEMENTO I NA DIREÇÃO J
        for j in range(self.NEL):
            for i in range(self.GLN):
                self.G[j][i]=       self.GLN*(self.NO[0][j]-1)+i+1
                self.G[j][int(i+self.GLN)]=self.GLN*(self.NO[1][j]-1)+i+1

        self.IL = np.array(np.zeros((self.NEL, 6)))
    
    def getElementProperties(self, elem_number):
        return [
            self.E [self.NC[elem_number]-1],
            self.A [self.NC[elem_number]-1],
            self.IZ[self.NC[elem_number]-1],
            self.FP[elem_number],
            self.QX[elem_number],
            self.QY[elem_number],
        ]

    def addElementStiffMatrixToGlobalStiffMatrix(self, RSR, elem_number):
        for i in range(6):
            for j in range(6):
                Gi = self.G[elem_number][i]-1
                Gj = self.G[elem_number][j]-1
                self.SG[Gi][Gj] = self.SG[Gi][Gj] + RSR[i][j]
    
    def addElementEqvForcesToGlobalForceVector(self, RtF, elem_number):
        for i in range(6): 
            Gi = self.G[elem_number][i]-1
            self.F[Gi]=self.F[Gi]+RtF[i]
    
    def addNodalForcesToForceVector(self):
        for i in range(self.NNO):
            for j in range(self.GLN):
                self.F[self.D[i][j]-1]= self.F[self.D[i][j]-1] + self.FA[i][j]

    def recordStiffnessMatrix(self):
        self.SG0 = self.SG
    
    def recordForceVector(self):
        self.F0 = self.F
    
    def applyBonds(self):
        for i in range(self.NNO):
            for j in range(self.GLN):
                if self.V[i][j] == 1:
                    for k in range(self.NGL):
                        self.SG[k][self.D[i][j]-1] = 0
                        self.SG[self.D[i][j]-1][k] = 0
                    self.SG[self.D[i][j]-1][self.D[i][j]-1] = 1
                    self.F[self.D[i][j]-1] = 0

    # def NO(self, i, j):  
        # return int(self.NO[i][j])

    # def X(self, i):   
        pass

    # def Y(self, i):  # COORDENADA Y DO NÓ I 
        pass

    # def A(self, i):  # ÁREA DA SEÇÃO DA CARACTERÍSTICA I
        pass

    # def E(self, i):  # MÓDULO DE ELASTICIDADE DA CARACTERÍSTICA I
        pass

    # def IZ(self, i):  # INÉRCIA DA CARACTERÍSTICA I
        pass

    # def QX(self, i):  # CARGA DISTRIBUÍDA DIR. X NO ELEMENTO I
        pass

    # def QY(self, i):  # CARGA DISTRIBUÍDA DIR. Y NO ELEMENTO I
        pass

    # def FP(self, i):  # FORCA DE PROTENSÇÃO NO ELEMENTO I
        pass

    # def NV(self, i):  # NÚMERO DO I-ÉSIMO NÓ VINCULADO
        pass

    # def TE(self, i):  # TIPO DO ELEMENTO I (PÓRTICO=1, TRELIÇA=2)
        pass

    # def NC(self, i):  # NúMERO DA CARACTERÍSTICA DO ELEMENTO I
        pass

    # def D(self, i, j):  # DIREÇÃO DO GRAU DE LIBERDADE DO NÓ I NA DIREÇÃO J
        pass

    # def G(self, i, j):  # G.L. GLOBAL DO ELEMENTO I NA DIREÇÃO J
        pass

    def IL(self, i, j):  # INDICADOR DE LIBERACAO LOCAL NÓ I, DIR. J (0=LIVRE,1=IMP)
        pass

    # def V(self, i, j):  # VÍNCULO DO NÓ I NA DIREÇÃO J (0=LIVRE, 1=IMPEDIDO)
        pass
