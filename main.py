from model import Model
import anest


model = Model('input.json')
model.plot()
anest.createSystem(model)
anest.solveSystem(model)
anest.calcInternLoads(model)
anest.calcReactions(model)
model.plotNFD(0.01)
model.plotSFD(0.01)
model.plotBMD(0.01, 10)
model.plotDisplacement(100, 10)
