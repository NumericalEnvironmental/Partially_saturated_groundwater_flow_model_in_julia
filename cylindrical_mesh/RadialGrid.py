################################################################################
#
# RadialGrid.py - generate cyclindrical grid (nodes and connections) for jFlow
#
# by Walt McNab (2021)
#
#######################################################################


from numpy import *
import pandas as pd


### model classes ###


class Connection:           # cell-to-cell connections
    
    def __init__(self, radial, layerModel):
        
        self.node0 = []     # connection node indices
        self.node1 = []
        self.d0 = []        # distances from connecting nodes to interface
        self.d1 = []
        self.area = []      # connection interfacial area
        self.baseArea = []  # reference base area (to be passed to Node object arrays)
        self.zConnect = []  # vertical connection (boolean)
        self.deltaB = 0.0               # distance into boundary node (for calculating gradient)
        
        # horizontal grid spacing
        rFace = self.Gridder(radial)                                # array of outside grid cell interface radii
        self.rNode = 0.5*rFace[1:] + 0.5*rFace[:-1]           # radius of node point associated with each cell
        self.rNode = insert(self.rNode, 0, 0.)                          # add cell representing well
       
        # connection indices arrays and associated parameters
        for i in range(radial.numLayers):          # horizontal connections
            for j in range(radial.numCols-1):
                node0 = j + i*radial.numCols
                self.node0.append(node0)
                self.node1.append(node0 + 1)
                self.d0.append(rFace[j] - self.rNode[j])
                self.d1.append(self.rNode[j+1] - rFace[j])
                self.area.append(2.*pi*rFace[j]*(layerModel.top[i]-layerModel.bottom[i]))
                self.zConnect.append(False)
              
            # horizontal exterior boundary connections
            self.node0.append(radial.numCols-1 + i*radial.numCols)
            self.node1.append(radial.numNodes - 1)
            self.d0.append(rFace.max() - self.rNode.max())
            self.d1.append(self.deltaB)
            self.area.append(2.*pi*radial.boundaryRadius*(layerModel.top[i]-layerModel.bottom[i]))
            self.zConnect.append(False)         

        # vertical connections
        rB = rFace
        rB = insert(rB, 0, 0.)
        rB = append(rB, radial.boundaryRadius)
        for i in range(radial.numCols): self.baseArea.append(pi * (rB[i+1]**2 - rB[i]**2))
        for j in range(radial.numLayers-1): self.area.extend(self.baseArea)     
        for i in range(radial.numCols):
            for j in range(radial.numLayers-1):
                node0 = j*radial.numCols + i
                self.node0.append(node0)
                self.node1.append(node0 + radial.numCols)
                d0 = 0.5*(layerModel.top[j] - layerModel.bottom[j])
                d1 = 0.5*(layerModel.top[j+1] - layerModel.bottom[j+1])
                self.d0.append(d0)                 
                self.d1.append(d1)
                self.zConnect.append(True)
                
        print('Set up numerical model grid cell connections.')
        
    def Gridder(self, radial):
        # generate radial grid
        index = arange(0, radial.numCols, 1)
        f = 10.**(log10(radial.boundaryRadius/radial.wellRadius)/(radial.numCols-1))   # sequential scaling factor
        r = radial.wellRadius * f**index
        return r        
        
    def WriteConnects(self):
        # convert to data frame and write to file
        self.node0 = array(self.node0)
        self.node1 = array(self.node1)
        output = pd.DataFrame(data={'node0':self.node0+1, 'node1':self.node1+1, 'del1':self.d0,
            'del2':self.d1, 'A':self.area, 'zConnect':self.zConnect})
        output['zConnect'] = output['zConnect'].astype(str)        
        output['zConnect'] = output['zConnect'].str.lower()
        output.to_csv('connections.csv', index=False)
        print('Wrote connection information to file.')


class Node:         # volume element properties
    
    def __init__(self, radial, layerModel, connection):
        self.index = []
        self.r = []
        self.z = []
        self.b = []
        self.h0 = []
        self.material = []
        self.vol = []
        for i in range(radial.numLayers):
            self.r.extend(connection.rNode)
            for j in range(radial.numCols):
                self.index.append(j + i*radial.numCols)
                self.z.append(0.5*layerModel.bottom[i] + 0.5*layerModel.top[i])
                self.b.append(layerModel.b[i])
                self.h0.append(radial.h0)
                if j: self.material.append(layerModel.mat[i])
                else: self.material.append(radial.wellMaterial)
                self.vol.append(connection.baseArea[j] * (layerModel.top[i]-layerModel.bottom[i]))
                
        # boundary node
        dVert = max(layerModel.top) - min(layerModel.bottom)
        self.index.append(radial.numNodes-1)
        self.r.append(radial.boundaryRadius + connection.deltaB)
        self.z.append(min(layerModel.bottom) + 0.5*dVert)
        self.b.append(dVert)
        self.h0.append(radial.h0)
        self.material.append(radial.boundaryMaterial)
        self.vol.append(1e+30)

        print('Set up numerical model cells.')


    def WriteNodes(self):
        # convert to data frame and write to file
        self.index = array(self.index)
        y = zeros(len(self.index), float)
        output = pd.DataFrame(data={'cellIndex':self.index+1, 'x':self.r, 'y':y, 'z':self.z,
            'vol':self.vol, 'b':self.b, 'h0':self.h0, 'material':self.material})
        output.to_csv('nodes.csv', index=False)
        print('Wrote grid cell information to file.')


class Radial:        # radial flow model grid characteristics
    
    def __init__(self):
        # read radial grid constraints from file
        lineInput = []
        inputFile = open('radial_params.txt','r')
        for line in inputFile:
            lineInput.append(line.split()[1])
        inputFile.close()
        self.numLayers = int(lineInput[0])
        self.numCols = int(lineInput[1])
        self.wellRadius = float(lineInput[2])
        self.wellMaterial = int(lineInput[3])
        self.boundaryRadius = float(lineInput[4])
        self.boundaryMaterial = int(lineInput[5])
        self.h0 = float(lineInput[6])
        self.numNodes = self.numLayers * self.numCols + 1
        print('Read radial grid parameters.')


class LayerModel:      # subsurface layer structure, numbered from bottom to top
    
    def __init__(self, model):
        layers = pd.read_csv('layers.csv')
        self.bottom = array(layers['bottom'])               # layer structure
        self.top = array(layers['top'])      
        self.b = array(layers['b'])                        # layer thickness (set =0 for Van Genuchten formulation)
        self.mat = array(layers['material'])                     # layer material index number
        print('Read aquifer structure attributes.')


### main script ###


def RadialGrid():


    ### note convert indices to base 1


    # radial system attributes
    radial = Radial()

    # aquifer structure
    layerModel = LayerModel(radial)
    
    # set up radial grid
    connection = Connection(radial, layerModel)
    node = Node(radial, layerModel, connection)

    # write grid data to external files for reference
    node.WriteNodes()
    connection.WriteConnects()    
    
    print('Done.')
    

### run script ###
RadialGrid()
