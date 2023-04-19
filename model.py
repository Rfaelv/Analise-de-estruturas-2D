import json


class Model:
    def __init__(self):
        self.GLN = 3     # NÚMERO DE GRAUS DE LIBERDADE POR NÓ     

        with open('input.json', 'r') as data:
            data  = json.load(data)

            self.NELP = len(data['frameElements'])     # NÚMERO DE ELEMENTOS DE PÓRTICO
            self.NELT = len(data['trussElements'])     # NÚMERO DE ELEMENTOS DE TRELIÇA
            self.NEL = self.NELP + self.NELT           # NÚMERO DE ELEMENTOS
            self.NNO = len(data['nodes'])              # NÚMERO DE NÓS
            self.NGL = self.NNO * self.NGL             # NÚMERO DE GRAUS DE LIBERDADE
            self.NNV = len(data['fixedNodes'])         # NÚMERO DE NÓS VINCULADOS
            self.NLL = None                            # NÚMERO DE LIBERACOES LOCAIS
            self.NTC = len(data['elementProperties'])  # NÚMERO DE TIPOS DE CARACTERÍSTICAS DE ELEMENTO

# *                                                                      *
# * REAIS - VETORES                                                      *
# *  F....... VETOR DE FOR�AS GLOBAIS                                    *
# *  RA...... VETOR DE REA��ES DE APOIO = [SG].{U} - {F}                 *
# *  U....... VETOR DE DESLOCAMENTOS                                     *
# *                                                                      *
# * REAIS - MATRIZES                                                     *
# *  ESF(I,J) ESFOR�O DIRE��O I ELEMENTO J                               *
# *  FA(I,J). FOR�A APLICADA (CARGA NODAL) NO N� I NA DIRE��O J          *
# *  SG...... MATRIZ DE RIGIDEZ GLOBAL 

    def NO(self, i, j):  # NÓ I DO ELEMENTO J
        pass

    def X(self, i):  # COORDENADA X DO NÓ I 
        pass

    def Y(self, i):  # COORDENADA Y DO NÓ I 
        pass

    def A(self, i):  # ÁREA DA SEÇÃO DA CARACTERÍSTICA I
        pass

    def E(self, i):  # MÓDULO DE ELASTICIDADE DA CARACTERÍSTICA I
        pass

    def IZ(self, i):  # INÉRCIA DA CARACTERÍSTICA I
        pass

    def QX(self, i):  # CARGA DISTRIBUÍDA DIR. X NO ELEMENTO I
        pass

    def QY(self, i):  # CARGA DISTRIBUÍDA DIR. Y NO ELEMENTO I
        pass

    def FP(self, i):  # FORCA DE PROTENSÇÃO NO ELEMENTO I
        pass

    def NV(self, i):  # NÚMERO DO I-ÉSIMO NÓ VINCULADO
        pass

    def TE(self, i):  # TIPO DO ELEMENTO I (PÓRTICO=1, TRELIÇA=2)
        pass

    def NC(self, i):  # NúMERO DA CARACTERÍSTICA DO ELEMENTO I
        pass

    def D(self, i, j):  # DIREÇÃO DO GRAU DE LIBERDADE DO NÓ I NA DIREÇÃO J
        pass

    def G(self, i, j):  # G.L. GLOBAL DO ELEMENTO I NA DIREÇÃO J
        pass

    def IL(self, i, j):  # INDICADOR DE LIBERACAO LOCAL NÓ I, DIR. J (0=LIVRE,1=IMP)
        pass

    def V(self, i, j):  # VÍNCULO DO NÓ I NA DIREÇÃO J (0=LIVRE, 1=IMPEDIDO)
        pass
