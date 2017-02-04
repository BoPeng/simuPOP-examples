import simuPOP as sim
import math
from simuPOP.utils import *
# example 1
traj = simulateBackwardTrajectory(N=50000, fitness=[1, 1, 1],
     endGen=10000, endFreq=0.1)
# example 2
traj = simulateBackwardTrajectory(N=50000, fitness=[1, 1.001, 1.002],
     endGen=10000, endFreq=0.1)
# example 3
traj = simulateBackwardTrajectory(N=50000, fitness=[1, 0.999, 0.998],
     endGen=10000, endFreq=0.1)
# example 4
def fitnessFunc(gen, subPop):
    if gen > 8000:
        return (1, 0.999, 0.998)
    else:
        return (1, 1.001, 1.002)

traj = simulateBackwardTrajectory(N=50000, fitness=fitnessFunc,
     endGen=10000, endFreq=0.1)
# example 5
def demoFunc(gen):
    if gen < 5000:
        return 10000
    else:
        return int(10000*math.exp(0.000921*(gen-5000)))

traj = simulateBackwardTrajectory(N=demoFunc, fitness=[1, 1, 1],
     endGen=10000, endFreq=0.1)
# example 6
traj = simulateBackwardTrajectory(N=demoFunc, fitness=[1, 0.999, 0.998],
     endGen=10000, endFreq=0.1)
