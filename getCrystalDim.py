import numpy as np
from pylab import *
from find_intersections import seg_intersect

def plotPt(pointArray, vectOut = [0, 0, 0], marker = 'o'): # plot in x, y
    plot(pointArray[0], pointArray[1], marker = marker, markersize=4, color='r')
    plot([vectOut[0] + pointArray[0], pointArray[0]],
         [vectOut[1] + pointArray[1], pointArray[1]], 'b-')

def getUnitVector(vector):
    size = np.sqrt(np.sum([x**2 for x in vector]))
    return vector/size

# numerical parameters
m = 21. # number of points along crystal to use

# physical setup
# lengths in mm, angles in radians
L = 1160. # for now, assume vertical detector plane
D = 100.
x0 = 170. # 870.
detWidth = 30.
detLength = 30.
detTilt = np.deg2rad(10.) # flat detector for now

# points in space [x, y, z]
source = [0, D, 0]
detCenter = [x0, 0, 0]

nodalUnit = [np.cos(detTilt), np.sin(detTilt), 0]
delta = np.arctan(D/x0) # angle from source to detector w.r.t. x-axis
bragg = delta + detTilt
sourceToNodal = np.sqrt(D**2 + x0**2) * np.sin(bragg)
sourceImage = [2*sourceToNodal*np.sin(detTilt), 
               D - 2*sourceToNodal*np.cos(detTilt),
               0]

# get points on detector for multicone
mInds = [float(x) for x in np.arange(0, 2*m)]
crystalPoints = [[detWidth*cos(detTilt)*(i/m-1), 
                  detWidth*sin(detTilt)*(i/m-1), 0] for i in mInds]
crystalPoints = np.add(crystalPoints, detCenter)

raysOut, detPts, coneAxes = {}, {}, {}
for i, pt in enumerate(crystalPoints):
    rayOutUnit = getUnitVector(np.subtract(pt, sourceImage))
    raysOut[i] = rayOutUnit
    
    detHeight = (L - pt[0])*rayOutUnit[1]/rayOutUnit[0]
    detPt = [L, detHeight, 0]
    detPts[i] = detPt
    
    coneAxis = np.subtract(detPt, source)
    axisNormal = getUnitVector([-1*coneAxis[1], coneAxis[0], coneAxis[2]])
    
    arcCenter = seg_intersect(detPt, source, \
                           pt, np.add(axisNormal*1000.,pt))
    
show()