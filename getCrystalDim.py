import numpy as np
from pylab import *
from find_intersections import seg_intersect
from mpl_toolkits.mplot3d import Axes3D


def plotPt(pointArray, vectOut = [0, 0, 0], marker = 'o'): # plot in x, y, z
    xs, ys, zs = [[x] for x in pointArray]
    print(xs)
    plot(xs, ys, zs, marker = marker, markersize=4, color='r')
    # plot([vectOut[0] + pointArray[0], pointArray[0]],
    #      [vectOut[1] + pointArray[1], pointArray[1]], 'b-')

def getUnitVector(vector):
    size = np.sqrt(np.sum([x**2 for x in vector]))
    return vector/size
    
def getDistance(point1, point2):
    return np.sqrt(np.sum([(x-y)**2 for x in point1 for y in point2]))

def getCrystalDim(D, x0, colorMap = 'copper'):
    # numerical parameters
    m = 51. # number of points along crystal to use
    cmap = get_cmap(colorMap)
    offset = [x0, -1*D, 0] # for plotting substrate on the origin
    
    
    # physical setup
    # lengths in mm, angles in radians
    L = 1280. # for now, assume vertical detector plane
    # D = 240.
    # x0 = 290.
    detWidth = 60.
    detLength = 60.
    detTilt = np.deg2rad(0.) # flat detector for now

    def getArc(alpha, R, centerPoint = [0, 0, 0], detLength = detLength):
        [x0, y0, z0] = centerPoint
        maxTheta = np.arcsin(detLength/(2*R))
        def x(theta):
            return x0 - (1-np.cos(theta))*R*np.sin(alpha)
        def y(theta):
            return y0 + (1-np.cos(theta))*R*(np.cos(alpha))
        def z(theta):
            return z0 + R*np.sin(theta)
        theta = np.arange(-1*maxTheta,maxTheta,0.0001)
        curveCoords = [[x(t), y(t), z(t)] for t in theta]
        return curveCoords

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

    params = []
    raysOut, detPts, coneAxes = {}, {}, {}
    for i, pt in enumerate(crystalPoints):
        rayOutUnit = getUnitVector(np.subtract(pt, sourceImage))
        raysOut[i] = rayOutUnit
        
        detHeight = (L - pt[0])*rayOutUnit[1]/rayOutUnit[0]
        detPt = [L, detHeight, 0]
        detPts[i] = detPt
        
        coneAxis = np.subtract(detPt, source)
        coneAngle = np.arctan(coneAxis[1]/coneAxis[0])
        axisNormal = getUnitVector([-1*coneAxis[1], coneAxis[0], coneAxis[2]])

        arcCenter = seg_intersect(detPt, source,
                               pt, np.add(axisNormal*1000.,pt))
                               
        R = getDistance(arcCenter, pt)
        detCurve = getArc(coneAngle, R, np.subtract(pt,np.add(offset,source)))
        [detCurveX, detCurveY, detCurveZ] = [[x[0] for x in detCurve],
                                            [x[1] for x in detCurve],
                                            [x[2] for x in detCurve]]

        plot(detCurveX, detCurveZ, detCurveY, color=cmap(float(i)/len(crystalPoints)))
        params.append(R)
    
# plotPt(source, '*')
    
fig = figure()
ax = fig.add_subplot(111, projection='3d')
getCrystalDim(240, 580, 'spectral')
getCrystalDim(120, 290, 'copper')
# getCrystalDim(10,780,   'RdBu')
legend(['L = 24, x0 = 58', 'L = 12, x0 = 29'], loc='best')
show()
    