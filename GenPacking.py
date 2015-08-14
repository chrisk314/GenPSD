#!/usr/bin/python

import sys
import numpy as np
import numpy.random as npr
import math

# Takes the position of a particle. Returns indices of containing subdomains.
def getCubeIndex(pos):
    cubeIndex = []
    subCubes = []
    for i in range(lCube):
        if pos[2] > cubeData[i*lCube2].bounds[4] and pos[2] < cubeData[i*lCube2].bounds[5]:
            subCubes = [i*lCube2 + x for x in range(lCube2)]
            if not i == lCube - 1:
                if pos[2] > cubeData[(i+1)*lCube2].bounds[4]:
                    subCubes.extend([(i+1)*lCube2 + x for x in range(lCube2)])
            break
    for i in subCubes:
        if (pos[0] > cubeData[i].bounds[0] and pos[0] < cubeData[i].bounds[1])\
                and (pos[1] > cubeData[i].bounds[2] and pos[1] < cubeData[i].bounds[3]):
            cubeIndex.append(i)
    return cubeIndex

# Generate random position
def genPos():
    return [npr.rand() * Lx, npr.rand() * Ly, npr.rand() * Lz] 

# Check if particle crosses periodic boundaries
def checkPbc(pos, pbc):
    check = False
    if pos[0] < partRad[0]:
        pbc[0] = 1;
        check = True
    elif pos[0] > Lx - partRad[0]:
        pbc[0] = -1
        check = True
    if pos[1] < partRad[0]:
        pbc[1] = 1;
        check = True
    elif pos[1] > Ly - partRad[0]:
        pbc[1] = -1
        check = True
    if pos[2] < partRad[0]:
        pbc[2] = 1;
        check = True
    elif pos[2] > Lz - partRad[0]:
        pbc[2] = -1
        check = True
    return check
        
def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

# default values
gradFile = "Grading_ToyouraSand_q0.txt"
poros = 0.4
Lx = 2.e-3
Ly = 2.0e-3
Lz = 2.0e-3
minDia = 1.2e-4

if len(sys.argv) != 6:
    print "Usage: GenPacking.py <grading file> <porosity> <Lx> <Ly> <Lz>."
    print "using default options."
else:
    gradFile = sys.argv[1]
    poros = float(sys.argv[2])
    Lx = float(sys.argv[3])
    Ly = float(sys.argv[4])
    Lz = float(sys.argv[5])
    minDia = float(sys.argv[6])

domainVol = Lx * Ly * Lz
domainPartVol = domainVol * (1 - poros)

# Load grading data
gradData = np.genfromtxt(gradFile, usecols=(0,1,2))
gradData[:,:2] = np.divide(gradData[:,:2], 1000) # Convert mm to m
gradData[:,2] = np.divide(gradData[:,2], 100) # Convert % to prob

# Eliminate small particle classes
for gradClass in gradData:
    if gradClass[1] < minDia:
        gradClass[2] = 0.0
    elif gradClass[1] > minDia and gradClass[0] < minDia:
        gradClass[2] *= (gradClass[1] - minDia) / (gradClass[1] - gradClass[0])
        gradClass[0] = minDia

meanDia = 0.5 * (gradData[:,0] + gradData[:,1])             # Class mean particle diameters
classPartVol = (math.pi/6) * meanDia[:]**3 * gradData[:,2]  # Class mean particle volumes
meanPartVol = sum(classPartVol)                             # Global mean particle volume
estNumParts = domainPartVol / meanPartVol                   # Estimated number of particles
classVol = estNumParts * classPartVol                       # Class volumes

# Generate particle diameters
partDia = []
overflowVol = 0.0
for i in reversed(range(len(gradData))):
    accumVol = 0.0
    if classVol[i] > 0.0:
        while accumVol < classVol[i] + overflowVol:
            partDia.append(gradData[i,0] + npr.rand() * (gradData[i,1] - gradData[i,0]))
            addVol = (math.pi/6) * partDia[-1]**3
            accumVol += addVol
        if accumVol - classVol[i] > classVol[i] - (accumVol - addVol):
            partDia.pop()
            overflowVol += classVol[i] - (accumVol - addVol)

# Sort diameters from largest to smallest
partDia.sort()
partDia.reverse()
partDia = np.array(partDia)
partRad = 0.5*partDia
partPos = np.zeros((len(partDia),3))

volActual = (math.pi/6) * sum(partDia**3)

print("%d particles generated with volume %+1.4e within %+1.4e of target volume\n" \
        %(len(partDia), volActual, abs(domainPartVol - volActual)/domainPartVol))

# Split domain into subdomains to reduce neighbour search overhead
lCube = 4
lCube2 = lCube**2
numCubes = lCube**3

# Subdomain lengths
Sx = Lx/lCube
Sy = Ly/lCube
Sz = Lz/lCube

tol = 1.0e-8*partDia[-1]

estNumPartsCube = math.floor((1.4 * estNumParts) / numCubes)

class Cube:
    def __init__(self, index):
        self.partData = np.empty((estNumPartsCube,4))
        self.bounds = np.empty(6)
        self.numParts = 0
        self.bounds[0] = (index % lCube) * Sx - (partRad[0] + tol);
        self.bounds[2] = (math.floor(index/lCube) % lCube) * Sy - (partRad[0] + tol);
        self.bounds[4] = (math.floor(index/lCube2)) * Sz - (partRad[0] + tol);
        self.bounds[1] = self.bounds[0] + Sx + partDia[0] + tol;
        self.bounds[3] = self.bounds[2] + Sy + partDia[0] + tol;
        self.bounds[5] = self.bounds[4] + Sz + partDia[0] + tol;

    def addPart(self, data):
        if self.numParts == self.partData.shape[0]:
            self.partData = np.vstack(self.partData, data)
        else:
            self.partData[self.numParts] = data
        self.numParts += 1

    def checkOverlaps(self):
        if self.numParts < 2:
            return False
        else:
            for i in range(self.numParts-1):
                dist = sum((self.partData[i,1:] - self.partData[-1,1:])**2)
                if dist - (self.partData[i,0] + self.partData[-1,0] + tol) < 0.0:
                    self.numParts -= 1
                    return True
            return False

cubeData = [Cube(i) for i in range(numCubes)]

# All possible ghost particle translation unit vectors
pbcPerms = []
for i in range(2):
    for j in range(2):
        for k in range(2):
            pbcPerms.append([i, j, k])

# Fill domain with particles
placed = 0
for i in range(len(partDia)):
    retry = True
    
    # Attempt to place particle. Loop until success.
    while retry:
        pos = genPos()                  # Generate a random position.
        cubeIndex = getCubeIndex(pos)   # Determine subdomain membership.
        
        # Check for overlaps with particles in containing subdomains.
        for j in cubeIndex:
            cubeData[j].addPart([partRad[i], pos[0], pos[1], pos[2]])
            retry = cubeData[j].checkOverlaps()
            if retry:
                break
            
        # Determine if particle crosses PBCs
        if not retry:
            pbc = np.zeros(3)
            if checkPbc(pos, pbc):
                # Determine ghost translations.
                posPerms = np.zeros((8,3))
                posPerms = pos + np.multiply(np.multiply(pbcPerms, pbc), [Lx, Ly, Lz])
                posPerms = unique_rows(posPerms)
                posPerms = np.delete(posPerms, np.where((posPerms == pos).all(axis=1)), axis=0) 
                
                # Check for overlaps in ghost positions
                for ghost in posPerms:
                    cubeIndex = getCubeIndex(ghost)
                    for j in cubeIndex: 
                        cubeData[j].addPart([partRad[i], ghost[0], ghost[1], ghost[2]])
                        retry = cubeData[j].checkOverlaps()
                        if retry:
                            break
                    if retry:
                        break
        
        # Succesful placement
        if not retry:
            partPos[i] = pos
            print ("placed particle %d" %(i))
            print pos
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
