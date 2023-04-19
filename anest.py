import numpy as np


class Anest:
    def __init__(self):
        self.input_file = input('Arquivo de dados (sem extensão):')

        #preliminary intern value attribution
        self.NEL = 3  # Número de elementos
        self.NNO = 3  # Número de nós
        self.NTC = 1  # Número de tipos de características de elemento
        self.GLN = 2  # Graus de liberdade por nó
        self.NO = np.array(np.zeros((2,self.NEL)))  # Nó i do elemento j
        self.NO[0][0] = 1
        self.NO[1][0] = 2
        self.NO[0][1] = 2
        self.NO[1][1] = 3
        self.NO[0][2] = 1
        self.NO[1][2] = 3
        self.X = np.array(np.zeros(self.NNO)) # Coordenada x do nó i
        self.X[0] = 0
        self.X[1] = 2
        self.X[2] = 0
        self.Y = np.array(np.zeros(self.NNO)) # Coordenada y do nó i
        self.Y[0] = 0
        self.Y[1] = 0
        self.Y[2] = 1.5
        self.R = np.array(np.zeros((6,6)))  # Matriz de rotação
        self.TE = np.array(np.zeros(self.NEL))  # Tipo de elemento (pórtico=1, treliça=2)
        self.TE[0] = 2
        self.TE[1] = 2
        self.TE[2] = 2

        self.FP = np.array(np.zeros(self.NEL))  # Força de protensão no elemento I
        self.QX = np.array(np.zeros(self.NEL))  # Carga distribuída dir. x elemento I
        self.QY = np.array(np.zeros(self.NEL))  # Carga distribuída dir. y elemento I
        self.FP[0] = 0
        self.FP[1] = 0
        self.FP[2] = 0

        self.A = np.array(np.zeros(self.NTC))  # Área da seção da característica I
        self.E = np.array(np.zeros(self.NTC))  # Módulo de elasticidade da característica I
        self.IZ = np.array(np.zeros(self.NTC))  # Inécia da característica I
        self.A[0] = 0.0015
        self.E[0] = 200*10e6
        self.IZ[0] = 200
        # self.A[1] = 0.02
        # self.E[1] = 20e10
        # self.IZ[1] = 300

        self.NC = np.array(np.zeros(self.NEL))  # Número da característica do elemento I
        self.NC[0] = 0
        self.NC[1] = 0
        self.NC[2] = 0

        self.NGL = 6  # Número de graus de lioberdade
        self.SG = np.array(np.zeros((self.NGL, self.NGL)))  # Matriz de rigidez global
        self.F = np.array(np.zeros((self.NGL)))  # Vetor de forças global
        self.G = np.array(np.zeros((self.NEL, 6)))  # G.L. global do elemento I na direção J 
        # self.G[0][0] = 3
        # self.G[0][1] = 4
        # self.G[0][2] = 1
        # self.G[0][3] = 2
        # self.G[0][1] = 2
        # self.G[0][1] = 2
        # self.G[0][1] = 2
        # self.G[1][0] = 1
        # self.G[2][0] = 1

        self.D = np.array(np.zeros((self.NNO, self.GLN)))  # Direção do grau de liberdade no nó I na direção J 
        # self.D[0][0] = 1
        # self.D[0][1] = 2
        # self.D[1][0] = 3
        # self.D[1][1] = 4
        # self.D[2][0] = 5
        # self.D[2][1] = 6

        self.FA = np.array(np.zeros((self.NNO, self.GLN)))  # Força aplicada (carga nodal) no nó I na direção J  
        self.FA[0][1] = -20
        self.FA[1][1] = -45

    def preparq(self):
        pass

    def nos(self):
        # FORMAR MATRIZ D COM NUMERA��O DAS DIRECOES DOS G.L. DO N�S
        for i in range(self.NNO):
            for j in range(self.GLN):
                self.D[i][j]=self.GLN*(i)+j+1 
        print(f'D: {self.D}')

    def elementos(self):
        # FORMAR MATRIZ G - NUMERA��O DOS GRAUS DE LIBERDADE DOS ELEMENTOS
        for j in range(self.NEL):
            for i in range(self.GLN):
                self.G[j][i]=       self.GLN*(self.NO[0][j]-1)+i+1
                self.G[j][int(i+self.GLN)]=self.GLN*(self.NO[1][j]-1)+i+1
        print(f'G: {self.G}')

    def sistema(self):
        self.S = np.array(np.zeros((6,6)))  # Matriz de rigides do elemento
        self.FEQV = np.array(np.zeros(6))  # Vetor de forças nodais equivalentes no sistema local
        FEQVG = np.array(np.zeros(6))  # Vetor de forças nodais equivalentes no sistema global
        self.L = None
        SR = np.array(np.zeros((6,6)))  # Produto S.R
        RSR = np.array(np.zeros((6,6)))  # Produto R'.SR

        for elem_number in range(self.NEL):
            R, L = self.rotação(self.X, self.Y, self.NO, elem_number)
            E = self.E[int(self.NC[elem_number])]
            A = self.A[int(self.NC[elem_number])]
            FP = self.FP[elem_number]
            IZ = self.IZ[int(self.NC[elem_number])]
            QX = self.QX[elem_number]
            QY = self.QY[elem_number]

            if self.TE[elem_number]==1:
                S, FEQV = self.port2d(E, A, IZ, L, QX, QY, FP)
                
            elif self.TE[elem_number]==2:
                S, FEQV = self.treli2d(E, A, L, FP)
            print(f'S: {S}')
            # MULTIPLICACAO SR=S.R
            for i in range(6):
                for j in range(6):
                    SR[i][j]=0
                    for k in range(5):
                        SR[i][j]=SR[i][j]+S[i][k]*R[k][j]
            print(f'SR: {SR}')
            # MULTIPLICACAO RSR=R'.SR
            for i in range(6):
                for j in range(6):
                    RSR[i][j]=0
                    for k in range(5):
                        RSR[i][j]=RSR[i][j]+R[k][i]*SR[k][j]
            print(f'RSR: {RSR}')
            # ADICIONAR MATRIZ DE RIGIDEZ DO ELEMENTO NA MATRIZ GLOBAL
            for i in range(6):
                for j in range(6):
                    self.SG[int(self.G[elem_number][i]-1)][int(self.G[elem_number][j]-1)] = self.SG[int(self.G[elem_number][i]-1)][int(self.G[elem_number][j]-1)] + RSR[i][j]

            # FORCAS EQUIVALENTES NO REFERENCIAL GLOBAL
            for i in range(6):
                FEQVG[i]=0
                for j in range(6):
                    FEQVG[i]=FEQVG[i]+R[j][i]*FEQV[j]
            
            # ADICIONAR FORCAS NODAIS EQUIVALENTES DO ELEMENTO EM {F}
            for i in range(6): 
                self.F[int(self.G[elem_number][i]-1)]=self.F[int(self.G[elem_number][i]-1)]+FEQVG[i]

        # ADICIONAR FORCAS NODAIS NO VETOR DE FOR�AS {F}
        for i in range(self.NNO):
            for j in range(self.GLN):
                self.F[int(self.D[i][j]-1)]= self.F[int(self.D[i][j]-1)] + self.FA[i][j]
        
        print(f'SG: {self.SG}')
        print(f'F: {self.F}')

    def port2d(self, E, A, IZ, L, QX, QY, FP):
        S = np.array(np.zeros((6,6)))
        FEQV = np.array(np.zeros(6))

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

        FEQV[0] = QX*L/2     + FP
        FEQV[1] = QY*L/2
        FEQV[2] = QY*L**2/12
        FEQV[3] = QX*L/2     - FP
        FEQV[4] = QY*L/2
        FEQV[5] =-QY*L**2/12

    def treli2d(self, E, A, L, FP):
        S = np.array(np.zeros((6,6)))
        FEQV = np.array(np.zeros(6))

        S[0][0] = E*A/L
        S[0][3] =-E*A/L
        S[3][0] =-E*A/L
        S[3][3] = E*A/L

        FEQV[0] = FP
        FEQV[1] = 0
        FEQV[2] = 0
        FEQV[3] =-FP
        FEQV[4] = 0
        FEQV[5] = 0

        return [S, FEQV]

    def rotação(self, x_coord, y_coord, node_number, elem_number):
        R = np.array(np.zeros((6,6)))
        
        dx = x_coord[int(node_number[1][elem_number]-1)] - x_coord[int(self.NO[0][elem_number]-1)]
        dy = y_coord[int(node_number[1][elem_number]-1)] - y_coord[int(self.NO[0][elem_number]-1)]
        L = (dx**2 + dy**2)**(0.5)
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
