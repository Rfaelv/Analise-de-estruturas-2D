from model import Model
from anest import Anest


# backbone
model = Model(path='input.json')
# model.plot()
anest = Anest(model)
anest.createSystem()
print(model.SG0)
# anest.solveSystem()
# anest.calcLoads()
# anest.plotResults()
