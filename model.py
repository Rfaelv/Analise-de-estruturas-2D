import json
import numpy as np
import matplotlib.pyplot as plt
import math


class Model:
    def __init__(self, path):

        with open(path, 'r') as self.data:
            self.data  = json.load(self.data)

        self.GLN = self.data['nDegreeOfFreedom']        # NÚMERO DE GRAUS DE LIBERDADE POR NÓ     
        self.NELP = len(self.data['frameElements'])     # NÚMERO DE ELEMENTOS DE PÓRTICO
        self.NELT = len(self.data['trussElements'])     # NÚMERO DE ELEMENTOS DE TRELIÇA
        self.NEL = self.NELP + self.NELT           # NÚMERO DE ELEMENTOS
        self.NNO = len(self.data['nodes'])              # NÚMERO DE NÓS
        self.NGL = self.NNO * self.GLN             # NÚMERO DE GRAUS DE LIBERDADE
        self.NNV = len(self.data['fixedNodes'])         # NÚMERO DE NÓS VINCULADOS
        self.NLL = None                            # NÚMERO DE LIBERACOES LOCAIS
        self.NTC = len(self.data['elementProperties'])  # NÚMERO DE TIPOS DE CARACTERÍSTICAS DE ELEMENTO
        self.SG = np.array(np.zeros((self.NGL, self.NGL)))  # MATRIZ DE RIGIDEZ GLOBAL
        self.F = np.array(np.zeros((self.NGL)))             # VETOR DE FORÇAS GLOBAIS
        self.U = np.array(np.zeros(self.NGL))               # VETOR DE DESLOCAMENTOS
        self.RA = np.array(np.zeros((self.NGL)))            # VETOR DE REAÇÕES DE APOIO = [SG].{U} - {F}
        self.ESF = np.array(np.zeros((6, self.NEL)))        # ESFORÇO DIREÇÃO I ELEMENTO J  
        
        self.FA = np.array(np.zeros((self.NNO, self.NGL)))  # FORÇA APLICADA (CARGA NODAL) NO NÓ I NA DIREÇÃO J
        for nodalForce in self.data['nodalForces']:
            for i in range(self.GLN):
                self.FA[nodalForce[0]-1][i] = nodalForce[i+1]

        self.NO = np.array(np.zeros((2, self.NEL)), dtype=np.int)  # NÓ I DO ELEMENTO J
        self.QX = np.array(np.zeros(self.NEL))  # CARGA DISTRIBUÍDA DIR. X NO ELEMENTO I
        self.QY = np.array(np.zeros(self.NEL))  # CARGA DISTRIBUÍDA DIR. Y NO ELEMENTO I
        self.FP = np.array(np.zeros(self.NEL))  # FORCA DE PROTENSÇÃO NO ELEMENTO I
        self.NC = np.array(np.zeros(self.NEL), dtype=np.int)  # NúMERO DA CARACTERÍSTICA DO ELEMENTO I
        self.TE = np.array(np.zeros(self.NEL), dtype=np.int)  # TIPO DO ELEMENTO I (PÓRTICO=1, TRELIÇA=2)
        for element in self.data['frameElements']:
            for i in range(2):
                self.NO[i, element[0]-1] = int(element[i+1])
            self.NC[element[0]-1] = element[3]
            self.FP[element[0]-1] = element[4]
            self.QX[element[0]-1] = element[5]
            self.QY[element[0]-1] = element[6]
            self.TE[element[0]-1] = 1
        
        for element in self.data['trussElements']:
            for i in range(2):
                self.NO[i, element[0]-1] = int(element[i+1])
            self.NC[element[0]-1] = element[3]
            self.FP[element[0]-1] = element[4]
            self.TE[element[0]-1] = 2

        self.X = np.array(np.zeros(self.NNO))  # COORDENADA X DO NÓ I
        self.Y = np.array(np.zeros(self.NNO))  # COORDENADA Y DO NÓ I
        for i, node in enumerate(self.data['nodes']):
            self.X[i] = node[0]
            self.Y[i] = node[1]
        
        self.E = np.array(np.zeros(self.NTC))   # MÓDULO DE ELASTICIDADE DA CARACTERÍSTICA I
        self.A = np.array(np.zeros(self.NTC))   # ÁREA DA SEÇÃO DA CARACTERÍSTICA I
        self.IZ = np.array(np.zeros(self.NTC))  # INÉRCIA DA CARACTERÍSTICA I
        for i, elemPropertie in enumerate(self.data['elementProperties']):
            self.E[i] = elemPropertie[0]
            self.A[i] = elemPropertie[1]
            self.IZ[i] = elemPropertie[2]
        
        self.NV = np.array(np.zeros(len(self.data['fixedNodes'])), dtype=np.int)  # NÚMERO DO I-ÉSIMO NÓ VINCULADO
        self.V = np.array(np.zeros((self.NNO, self.GLN)), dtype=np.int)  # VÍNCULO DO NÓ I NA DIREÇÃO J (0=LIVRE, 1=IMPEDIDO)
        for i, fixedNode in enumerate(self.data['fixedNodes']):
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
        self.SG0 = np.copy(self.SG)
    
    def recordForceVector(self):
        self.F0 = np.copy(self.F)
    
    def applyBonds(self):
        for i in range(self.NNO):
            for j in range(self.GLN):
                if self.V[i][j] == 1:
                    for k in range(self.NGL):
                        self.SG[k][self.D[i][j]-1] = 0
                        self.SG[self.D[i][j]-1][k] = 0
                    self.SG[self.D[i][j]-1][self.D[i][j]-1] = 1
                    self.F[self.D[i][j]-1] = 0
    
    def elementDisplacementVector(self, elem_number):
        RtU = np.array(np.zeros(6))
        for i in range(6): 
            Gi = self.G[elem_number][i]-1
            RtU[i] = self.U[Gi]
        return RtU

    def plot(self):
        plt.axes()

        for frameElement in self.data['frameElements']:
            xi = self.data['nodes'][frameElement[1]-1][0]
            xf = self.data['nodes'][frameElement[2]-1][0]
            yi = self.data['nodes'][frameElement[1]-1][1]
            yf = self.data['nodes'][frameElement[2]-1][1]
            plt.gca().add_line(plt.Line2D((xi, xf), (yi, yf)))
            plt.gca().annotate(frameElement[0], xy=((xi+xf)/2, (yi+yf)/2),
                 bbox=dict(boxstyle="round", fc="w"), ha="center", va='center', fontsize=8,
                 annotation_clip=True)
            
            tag = ''

            if frameElement[4] != 0:
                tag += f'FP:{frameElement[4]}'

            if frameElement[5] != 0:
                if tag != '':
                    tag += '; '
                tag += f'QX:{frameElement[5]}'

            if frameElement[6] != 0:
                if tag != '':
                    tag += '; '
                tag += f'QY:{frameElement[6]}'

            if tag != '':
                plt.gca().annotate(tag, xy=((xi+3*xf)/4, (yi+3*yf)/4),
                    bbox=dict(boxstyle="round", fc="w"), ha="center", va='center', fontsize=8,
                    annotation_clip=True, rotation=math.asin((yf-yi)/((xf-xi)**2+(yf-yi)**2)**(1/2))*180/math.pi)
            
        for trussElement in self.data['trussElements']:
            xi = self.data['nodes'][trussElement[1]-1][0]
            xf = self.data['nodes'][trussElement[2]-1][0]
            yi = self.data['nodes'][trussElement[1]-1][1]
            yf = self.data['nodes'][trussElement[2]-1][1]
            plt.gca().add_line(plt.Line2D((xi, xf), (yi, yf), marker="", markersize=10, markerfacecolor='r', 
                         markeredgecolor='r', color='b'))
            plt.gca().annotate(trussElement[0], xy=((xi+xf)/2, (yi+yf)/2),
                 bbox=dict(boxstyle="round", fc="w"), ha="center", va='center', fontsize=8,
                 annotation_clip=True)

            if trussElement[4] != '':
                plt.gca().annotate(f'FP: {trussElement[4]}', xy=((xi+3*xf)/4, (yi+3*yf)/4),
                    bbox=dict(boxstyle="round", fc="w"), ha="center", va='center', fontsize=8,
                    annotation_clip=True, rotation=math.asin((yf-yi)/((xf-xi)**2+(yf-yi)**2)**(1/2))*180/math.pi)
            
        for fixedNode in self.data['fixedNodes']:
            nodeCoord = self.data['nodes'][fixedNode[0]-1]
            if fixedNode[1] == 1 and fixedNode[2] == 1 and fixedNode[3] == 1:
                scale = 0.2
                line1  = plt.Line2D((nodeCoord[0]-1.2*scale, nodeCoord[0]+1.2*scale), (nodeCoord[1], nodeCoord[1]), lw=2)
                line2 = plt.Line2D((nodeCoord[0]-scale, nodeCoord[0]-1.5*scale), (nodeCoord[1], nodeCoord[1]-0.5*scale), lw=1)
                line3 = plt.Line2D((nodeCoord[0]-0.5*scale, nodeCoord[0]-scale), (nodeCoord[1], nodeCoord[1]-0.5*scale), lw=1)
                line4 = plt.Line2D((nodeCoord[0], nodeCoord[0]-0.5*scale), (nodeCoord[1], nodeCoord[1]-0.5*scale), lw=1)
                line5 = plt.Line2D((nodeCoord[0]+0.5*scale, nodeCoord[0]), (nodeCoord[1], nodeCoord[1]-0.5*scale), lw=1)
                line6 = plt.Line2D((nodeCoord[0]+scale, nodeCoord[0]+0.5*scale), (nodeCoord[1], nodeCoord[1]-0.5*scale), lw=1)
                plt.gca().add_line(line1)
                plt.gca().add_line(line2)
                plt.gca().add_line(line3)
                plt.gca().add_line(line4)
                plt.gca().add_line(line5)
                plt.gca().add_line(line6)
            
            if fixedNode[1] == 0 and fixedNode[2] == 1 and fixedNode[3] == 0:
                scale = 0.2
                line1  = plt.Line2D((nodeCoord[0]-1.2*scale, nodeCoord[0]+1.2*scale), (nodeCoord[1]-1.2*scale, nodeCoord[1]-1.2*scale), lw=2)
                line2  = plt.Line2D((nodeCoord[0]-scale, nodeCoord[0]+scale), (nodeCoord[1]-scale, nodeCoord[1]-scale), lw=1)
                line3  = plt.Line2D((nodeCoord[0], nodeCoord[0]+scale), (nodeCoord[1], nodeCoord[1]-scale), lw=1)
                line4  = plt.Line2D((nodeCoord[0], nodeCoord[0]-scale), (nodeCoord[1], nodeCoord[1]-scale), lw=1)
                plt.gca().add_line(line1)
                plt.gca().add_line(line2)
                plt.gca().add_line(line3)
                plt.gca().add_line(line4)
            
            if fixedNode[1] == 1 and fixedNode[2] == 1 and fixedNode[3] == 0:
                scale = 0.2
                line2  = plt.Line2D((nodeCoord[0]-1.2*scale, nodeCoord[0]+1.2*scale), (nodeCoord[1]-scale, nodeCoord[1]-scale), lw=2)
                line3  = plt.Line2D((nodeCoord[0], nodeCoord[0]+scale), (nodeCoord[1], nodeCoord[1]-scale), lw=1)
                line4  = plt.Line2D((nodeCoord[0], nodeCoord[0]-scale), (nodeCoord[1], nodeCoord[1]-scale), lw=1)
                line4  = plt.Line2D((nodeCoord[0], nodeCoord[0]-scale), (nodeCoord[1], nodeCoord[1]-scale), lw=1)
                line5 = plt.Line2D((nodeCoord[0]-scale, nodeCoord[0]-1.5*scale), (nodeCoord[1]-scale, nodeCoord[1]-1.5*scale), lw=1)
                line6 = plt.Line2D((nodeCoord[0]-0.5*scale, nodeCoord[0]-scale), (nodeCoord[1]-scale, nodeCoord[1]-1.5*scale), lw=1)
                line7 = plt.Line2D((nodeCoord[0], nodeCoord[0]-0.5*scale), (nodeCoord[1]-scale, nodeCoord[1]-1.5*scale), lw=1)
                line8 = plt.Line2D((nodeCoord[0]+0.5*scale, nodeCoord[0]), (nodeCoord[1]-scale, nodeCoord[1]-1.5*scale), lw=1)
                line9 = plt.Line2D((nodeCoord[0]+scale, nodeCoord[0]+0.5*scale), (nodeCoord[1]-scale, nodeCoord[1]-1.5*scale), lw=1)
                plt.gca().add_line(line2)
                plt.gca().add_line(line3)
                plt.gca().add_line(line4)
                plt.gca().add_line(line5)
                plt.gca().add_line(line6)
                plt.gca().add_line(line7)
                plt.gca().add_line(line8)
                plt.gca().add_line(line9)
            
            if fixedNode[1] == 1 and fixedNode[2] == 0 and fixedNode[3] == 0:
                scale = 0.2
                line2  = plt.Line2D((nodeCoord[0]-scale, nodeCoord[0]+scale), (nodeCoord[1]-scale, nodeCoord[1]-scale), lw=1)
                line3  = plt.Line2D((nodeCoord[0], nodeCoord[0]+scale), (nodeCoord[1], nodeCoord[1]-scale), lw=1)
                line4  = plt.Line2D((nodeCoord[0], nodeCoord[0]-scale), (nodeCoord[1], nodeCoord[1]-scale), lw=1)
                plt.gca().add_line(line1)
                plt.gca().add_line(line2)
                plt.gca().add_line(line3)
                plt.gca().add_line(line4)
        
        for i, node in enumerate(self.data['nodes']):
            scale = 0.15
            plt.gca().annotate(i+1, xy=(node[0]-scale, node[1]+scale),
                 bbox=dict(boxstyle="circle", fc="w"), ha="center", va='center', fontsize=8,
                 annotation_clip=True)
            
            fixedNodes = np.array(self.data['fixedNodes'])
            nFixedNodes = fixedNodes[:, 0]

            if i+1 not in nFixedNodes:
                plt.gca().add_patch(plt.Circle((node[0], node[1]), 0.05))
        
        for nodalForce in self.data['nodalForces']:
            scale = 0.5
            x = self.data['nodes'][nodalForce[0]-1][0]
            y = self.data['nodes'][nodalForce[0]-1][1]

            fx = nodalForce[1]
            fy = nodalForce[2]
            mz = nodalForce[3]
            if fx != 0:
                plt.annotate(fx,
                    horizontalalignment = 'center',
                    verticalalignment = 'center',
                    xytext = (x-fx/abs(fx)*scale, y),
                    xy = (x-0.1*fx/abs(fx)*scale, y),
                    arrowprops = dict(facecolor = 'black', shrink = 0.01, headwidth=6, headlength=6, width=1)
                )

            if fy != 0:
                plt.annotate(fy,
                    horizontalalignment = 'center',
                    verticalalignment = 'center',
                    xytext = (x, y-fy/abs(fy)*scale),
                    xy = (x, y-0.1*fy/abs(fy)*scale),
                    arrowprops = dict(facecolor = 'black', shrink = 0.01, headwidth=6, headlength=6, width=1)
                )
            
            if mz != 0:
                plt.annotate(mz,
                    horizontalalignment = 'center',
                    verticalalignment = 'center',
                    xytext = (x+mz/abs(mz)*scale, y+mz/abs(mz)*scale),
                    xy = (x+0.1*mz/abs(mz)*scale, y+0.1*mz/abs(mz)*scale),
                    arrowprops = dict(facecolor = 'black', shrink = 0.01, headwidth=6, headlength=6, width=1)
                )
                plt.annotate(mz,
                    horizontalalignment = 'center',
                    verticalalignment = 'center',
                    xytext = (x+mz/abs(mz)*scale, y+mz/abs(mz)*scale),
                    xy = (x+0.2*mz/abs(mz)*scale, y+0.2*mz/abs(mz)*scale),
                    arrowprops = dict(facecolor = 'black', shrink = 0.01, headwidth=6, headlength=6, width=1)
                )

        plt.axis('scaled')
        plt.axis('off')
        plt.show()
    
    def plotNFD(self, scale):
        plt.axes()

        for i, node in enumerate(self.data['nodes']):          
            fixedNodes = np.array(self.data['fixedNodes'])
            nFixedNodes = fixedNodes[:, 0]

            if i+1 not in nFixedNodes:
                plt.gca().add_patch(plt.Circle((node[0], node[1]), 0.05))

        for frameElement in self.data['frameElements']:
            xi = self.data['nodes'][frameElement[1]-1][0]
            xf = self.data['nodes'][frameElement[2]-1][0]
            yi = self.data['nodes'][frameElement[1]-1][1]
            yf = self.data['nodes'][frameElement[2]-1][1]
            plt.gca().add_line(plt.Line2D((xi, xf), (yi, yf)))

            Ni = self.ESF[0][frameElement[0]-1]
            Nf =-self.ESF[3][frameElement[0]-1]

            Nix = Ni * (yf-yi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)
            Niy =-Ni * (xf-xi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)

            Nfx = Nf * (yf-yi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)
            Nfy =-Nf * (xf-xi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)

            plt.gca().add_line(plt.Line2D((xi-Nix*scale, xf-Nfx*scale), (yi-Niy*scale, yf-Nfy*scale), lw=0.5, color='r'))
            plt.gca().add_line(plt.Line2D((xi, xi-Nix*scale), (yi, yi-Niy*scale), lw=0.5, color='r'))
            plt.gca().add_line(plt.Line2D((xf, xf-Nfx*scale), (yf, yf-Nfy*scale), lw=0.5, color='r'))

            plt.gca().annotate(f'{Ni:.2f}', xy=(xi-Nix*scale, yi-Niy*scale),
                ha="right", va='center', fontsize=8, bbox=dict(boxstyle="square", fc="b", edgecolor='b', alpha=0.7),
                annotation_clip=True, color='w', fontweight='bold')
            plt.gca().annotate(f'{Nf:.2f}', xy=(xf-Nfx*scale, yf-Nfy*scale),
                ha="right", va='bottom', fontsize=8, bbox=dict(boxstyle="square", fc="b", edgecolor='b', alpha=0.7),
                annotation_clip=True, color='w', fontweight='bold')
            
        for trussElement in self.data['trussElements']:
            xi = self.data['nodes'][trussElement[1]-1][0]
            xf = self.data['nodes'][trussElement[2]-1][0]
            yi = self.data['nodes'][trussElement[1]-1][1]
            yf = self.data['nodes'][trussElement[2]-1][1]
            plt.gca().add_line(plt.Line2D((xi, xf), (yi, yf), marker="", markersize=10, markerfacecolor='r', 
                         markeredgecolor='r', color='b'))

            Ni = self.ESF[0][trussElement[0]-1]
            Nf =-self.ESF[3][trussElement[0]-1]

            Nix = Ni * (yf-yi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)
            Niy =-Ni * (xf-xi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)

            Nfx = Nf * (yf-yi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)
            Nfy =-Nf * (xf-xi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)

            plt.gca().add_line(plt.Line2D((xi-Nix*scale, xf-Nfx*scale), (yi-Niy*scale, yf-Nfy*scale), lw=0.5, color='r'))
            plt.gca().add_line(plt.Line2D((xi, xi-Nix*scale), (yi, yi-Niy*scale), lw=0.5, color='r'))
            plt.gca().add_line(plt.Line2D((xf, xf-Nfx*scale), (yf, yf-Nfy*scale), lw=0.5, color='r'))

            plt.gca().annotate(f'{Ni:.2f}', xy=(xi-Nix*scale, yi-Niy*scale),
                ha="right", va='center', fontsize=8, bbox=dict(boxstyle="square", fc="b", edgecolor='b', alpha=0.7),
                annotation_clip=True, color='w', fontweight='bold')
            plt.gca().annotate(f'{Nf:.2f}', xy=(xf-Nfx*scale, yf-Nfy*scale),
                ha="right", va='bottom', fontsize=8, bbox=dict(boxstyle="square", fc="b", edgecolor='b', alpha=0.7),
                annotation_clip=True, color='w', fontweight='bold')

        self.plotReactions()
        plt.axis('scaled')
        plt.axis('off')
        plt.show()

    def plotSFD(self, scale):
        plt.axes()

        for i, node in enumerate(self.data['nodes']):          
            fixedNodes = np.array(self.data['fixedNodes'])
            nFixedNodes = fixedNodes[:, 0]

            if i+1 not in nFixedNodes:
                plt.gca().add_patch(plt.Circle((node[0], node[1]), 0.05))

        for frameElement in self.data['frameElements']:
            xi = self.data['nodes'][frameElement[1]-1][0]
            xf = self.data['nodes'][frameElement[2]-1][0]
            yi = self.data['nodes'][frameElement[1]-1][1]
            yf = self.data['nodes'][frameElement[2]-1][1]
            plt.gca().add_line(plt.Line2D((xi, xf), (yi, yf)))

            Qi = self.ESF[1][frameElement[0]-1]
            Qf =-self.ESF[4][frameElement[0]-1]

            Qix = Qi * (yf-yi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)
            Qiy =-Qi * (xf-xi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)

            Qfx = Qf * (yf-yi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)
            Qfy =-Qf * (xf-xi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)

            plt.gca().add_line(plt.Line2D((xi-Qix*scale, xf-Qfx*scale), (yi-Qiy*scale, yf-Qfy*scale), lw=0.5, color='r'))
            plt.gca().add_line(plt.Line2D((xi, xi-Qix*scale), (yi, yi-Qiy*scale), lw=0.5, color='r'))
            plt.gca().add_line(plt.Line2D((xf, xf-Qfx*scale), (yf, yf-Qfy*scale), lw=0.5, color='r'))

            plt.gca().annotate(f'{-Qi:.2f}', xy=(xi-Qix*scale, yi-Qiy*scale),
                ha="right", va='center', fontsize=8, bbox=dict(boxstyle="square", fc="b", edgecolor='b', alpha=0.7),
                annotation_clip=True, color='w', fontweight='bold')
            plt.gca().annotate(f'{-Qf:.2f}', xy=(xf-Qfx*scale, yf-Qfy*scale),
                ha="right", va='bottom', fontsize=8, bbox=dict(boxstyle="square", fc="b", edgecolor='b', alpha=0.7),
                annotation_clip=True, color='w', fontweight='bold')
            
        self.plotReactions()
        plt.axis('scaled')
        plt.axis('off')
        plt.show()
    
    def plotBMD(self, scale, n):
        plt.axes()

        for i, node in enumerate(self.data['nodes']):          
            fixedNodes = np.array(self.data['fixedNodes'])
            nFixedNodes = fixedNodes[:, 0]

            if i+1 not in nFixedNodes:
                plt.gca().add_patch(plt.Circle((node[0], node[1]), 0.05))

        for frameElement in self.data['frameElements']:
            xi = self.data['nodes'][frameElement[1]-1][0]
            xf = self.data['nodes'][frameElement[2]-1][0]
            yi = self.data['nodes'][frameElement[1]-1][1]
            yf = self.data['nodes'][frameElement[2]-1][1]
            plt.gca().add_line(plt.Line2D((xi, xf), (yi, yf)))

            Mi = self.ESF[2][frameElement[0]-1]
            Mf =-self.ESF[5][frameElement[0]-1]

            Mix = Mi * (yf-yi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)
            Miy =-Mi * (xf-xi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)

            Mfx = Mf * (yf-yi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)
            Mfy =-Mf * (xf-xi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)

            xix = xi * (xf-xi) / ((xf-xi)**2+(yf-yi)**2)**(1/2) + yi * (yf-yi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)
            yiy =-xi * (yf-yi) / ((xf-xi)**2+(yf-yi)**2)**(1/2) + yi * (xf-xi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)

            xfx = xf * (xf-xi) / ((xf-xi)**2+(yf-yi)**2)**(1/2) + yf * (yf-yi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)
            yfy =-xf * (yf-yi) / ((xf-xi)**2+(yf-yi)**2)**(1/2) + yf * (xf-xi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)

            Qi = -self.ESF[1][frameElement[0]-1]

            a = (Mf-Mi-Qi*(xfx-xix))/(xfx**2-xix**2-2*(xix*xfx-xix**2))
            b = Qi-2*a*xix
            c = Mi-a*xix**2-b*xix

            for i in range(n):
                x1 = xix + (xfx-xix)/n*i
                x2 = xix + (xfx-xix)/n*(i+1)
                y1 = yiy + (yfy-yiy)/n*i
                y2 = yiy + (yfy-yiy)/n*(i+1)

                M1 = a*x1**2 + b*x1 + c
                M2 = a*x2**2 + b*x2 + c

                x1 = xi + (xf-xi)/n*i
                x2 = xi + (xf-xi)/n*(i+1)
                y1 = yi + (yf-yi)/n*i
                y2 = yi + (yf-yi)/n*(i+1)

                M1x = M1 * (y2-y1) / ((x2-x1)**2+(y2-y1)**2)**(1/2)
                M1y =-M1 * (x2-x1) / ((x2-x1)**2+(y2-y1)**2)**(1/2)

                M2x = M2 * (y2-y1) / ((x2-x1)**2+(y2-y1)**2)**(1/2)
                M2y =-M2 * (x2-x1) / ((x2-x1)**2+(y2-y1)**2)**(1/2)
                plt.gca().add_line(plt.Line2D((x1-M1x*scale, x2-M2x*scale), (y1-M1y*scale, y2-M2y*scale), lw=0.5, color='r'))
            plt.gca().add_line(plt.Line2D((xi, xi-Mix*scale), (yi, yi-Miy*scale), lw=0.5, color='r'))
            plt.gca().add_line(plt.Line2D((xf, xf-Mfx*scale), (yf, yf-Mfy*scale), lw=0.5, color='r'))

            plt.gca().annotate(f'{Mi:.2f}', xy=(xi-Mix*scale, yi-Miy*scale),
                ha="right", va='center', fontsize=8, bbox=dict(boxstyle="square", fc="b", edgecolor='b', alpha=0.7),
                annotation_clip=True, color='w', fontweight='bold')
            plt.gca().annotate(f'{Mf:.2f}', xy=(xf-Mfx*scale, yf-Mfy*scale),
                ha="right", va='bottom', fontsize=8, bbox=dict(boxstyle="square", fc="b", edgecolor='b', alpha=0.7),
                annotation_clip=True, color='w', fontweight='bold')
        self.plotReactions()
        plt.axis('scaled')
        plt.axis('off')
        plt.show()
    
    def plotReactions(self):
        for i, node in enumerate(self.data['nodes']):
            scale = 0.5
            x = node[0]
            y = node[1]

            fx = self.RA[(i+1)*3-3]
            fy = self.RA[(i+1)*3-2]
            mz = self.RA[(i+1)*3-1]
            if abs(fx) > 10e-5:
                plt.annotate(f'{fx:.2f}',
                    horizontalalignment = 'center',
                    verticalalignment = 'center',
                    xytext = (x-fx/abs(fx)*scale, y),
                    xy = (x-0.1*fx/abs(fx)*scale, y),
                    arrowprops = dict(facecolor = 'black', shrink = 0.01, headwidth=6, headlength=6, width=1)
                )

            if abs(fy) > 10e-5:
                plt.annotate(f'{fy:.2f}',
                    horizontalalignment = 'center',
                    verticalalignment = 'center',
                    xytext = (x, y-fy/abs(fy)*scale),
                    xy = (x, y-0.1*fy/abs(fy)*scale),
                    arrowprops = dict(facecolor = 'black', shrink = 0.01, headwidth=6, headlength=6, width=1)
                )
            
            if abs(mz) > 10e-5:
                plt.annotate(f'{mz:.2f}',
                    horizontalalignment = 'center',
                    verticalalignment = 'center',
                    xytext = (x+mz/abs(mz)*scale, y+mz/abs(mz)*scale),
                    xy = (x+0.1*mz/abs(mz)*scale, y+0.1*mz/abs(mz)*scale),
                    arrowprops = dict(facecolor = 'black', shrink = 0.01, headwidth=6, headlength=6, width=1)
                )
                plt.annotate(f'{mz:.2f}',
                    horizontalalignment = 'center',
                    verticalalignment = 'center',
                    xytext = (x+mz/abs(mz)*scale, y+mz/abs(mz)*scale),
                    xy = (x+0.2*mz/abs(mz)*scale, y+0.2*mz/abs(mz)*scale),
                    arrowprops = dict(facecolor = 'black', shrink = 0.01, headwidth=6, headlength=6, width=1)
                )
        
    def plotDisplacement(self, fact, n):
        deformedNodes = np.copy(self.data['nodes'])
        for i, node in enumerate(self.data['nodes']):
            scale = 0.5
            x = node[0]
            y = node[1]

            dx = self.U[(i+1)*3-3]
            dy = self.U[(i+1)*3-2]
            
            deformedNodes[i][0] = x+dx*fact
            deformedNodes[i][1] = y+dy*fact

        for frameElement in self.data['frameElements']:
            xi = self.data['nodes'][frameElement[1]-1][0]
            xf = self.data['nodes'][frameElement[2]-1][0]
            yi = self.data['nodes'][frameElement[1]-1][1]
            yf = self.data['nodes'][frameElement[2]-1][1]
            plt.gca().add_line(plt.Line2D((xi, xf), (yi, yf)))

            xi = deformedNodes[frameElement[1]-1][0]
            xf = deformedNodes[frameElement[2]-1][0]
            yi = deformedNodes[frameElement[1]-1][1]
            yf = deformedNodes[frameElement[2]-1][1]
            plt.gca().add_line(plt.Line2D((xi, xf), (yi, yf), marker='', color='r', ls='dashed', lw=1))

            # Mi = self.ESF[2][frameElement[0]-1]
            # Mf =-self.ESF[5][frameElement[0]-1]

            # Mix = Mi * (yf-yi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)
            # Miy =-Mi * (xf-xi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)

            # Mfx = Mf * (yf-yi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)
            # Mfy =-Mf * (xf-xi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)

            # xix = xi * (xf-xi) / ((xf-xi)**2+(yf-yi)**2)**(1/2) + yi * (yf-yi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)
            # yiy =-xi * (yf-yi) / ((xf-xi)**2+(yf-yi)**2)**(1/2) + yi * (xf-xi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)

            # xfx = xf * (xf-xi) / ((xf-xi)**2+(yf-yi)**2)**(1/2) + yf * (yf-yi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)
            # yfy =-xf * (yf-yi) / ((xf-xi)**2+(yf-yi)**2)**(1/2) + yf * (xf-xi) / ((xf-xi)**2+(yf-yi)**2)**(1/2)

            # Qi = -self.ESF[1][frameElement[0]-1]

            # a = (Mf-Mi-Qi*(xfx-xix))/(xfx**2-xix**2-2*(xix*xfx-xix**2))
            # b = Qi-2*a*xix
            # c = Mi-a*xix**2-b*xix
            # EI = self.E[frameElement[3]-1] * self.IZ[frameElement[3]-1]
            # d = self.U[frameElement[1]*3-1]*EI

            # for i in range(n):
            #     x1 = xix + (xfx-xix)/n*i
            #     x2 = xix + (xfx-xix)/n*(i+1)
            #     y1 = yiy + (yfy-yiy)/n*i
            #     y2 = yiy + (yfy-yiy)/n*(i+1)

            #     M1 = (a*x1**4/12 + b*x1**3/6 + c*x1**2/2+d*x1)/EI*fact
            #     M2 = (a*x2**4/12 + b*x2**3/6 + c*x2**2/2+d*x2)/EI*fact

            #     x1 = xi + (xf-xi)/n*i
            #     x2 = xi + (xf-xi)/n*(i+1)
            #     y1 = yi + (yf-yi)/n*i
            #     y2 = yi + (yf-yi)/n*(i+1)

            #     M1x = M1 * (y2-y1) / ((x2-x1)**2+(y2-y1)**2)**(1/2)
            #     M1y =-M1 * (x2-x1) / ((x2-x1)**2+(y2-y1)**2)**(1/2)

            #     M2x = M2 * (y2-y1) / ((x2-x1)**2+(y2-y1)**2)**(1/2)
            #     M2y =-M2 * (x2-x1) / ((x2-x1)**2+(y2-y1)**2)**(1/2)
            #     plt.gca().add_line(plt.Line2D((x1-M1x, x2-M2x), (y1-M1y, y2-M2y), lw=0.5, color='r'))

            
            
        for trussElement in self.data['trussElements']:
            xi = self.data['nodes'][trussElement[1]-1][0]
            xf = self.data['nodes'][trussElement[2]-1][0]
            yi = self.data['nodes'][trussElement[1]-1][1]
            yf = self.data['nodes'][trussElement[2]-1][1]
            plt.gca().add_line(plt.Line2D((xi, xf), (yi, yf), marker="", markersize=10, markerfacecolor='r', 
                         markeredgecolor='r', color='b'))
            
            xi = deformedNodes[trussElement[1]-1][0]
            xf = deformedNodes[trussElement[2]-1][0]
            yi = deformedNodes[trussElement[1]-1][1]
            yf = deformedNodes[trussElement[2]-1][1]
            plt.gca().add_line(plt.Line2D((xi, xf), (yi, yf), marker="", markersize=10, markerfacecolor='r', 
                         markeredgecolor='r', color='r', ls='dashed', lw=1))
        plt.axis('scaled')
        plt.axis('off')
        plt.show()

    # def IL(self, i, j):  # INDICADOR DE LIBERACAO LOCAL NÓ I, DIR. J (0=LIVRE,1=IMP)
        # pass


