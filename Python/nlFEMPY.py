import timeit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
from matplotlib.cm import ScalarMappable


def map_DOF(x):
    return x*2, x*2+1

def gather_index(a, b):
    """return index pairs between index locations a and b, including a, b."""
    a1, a2 = a
    b1, b2 = b
    c = np.max(abs(np.array(a)-np.array(b))) + 1 # +1 to include last pair
    output = [None]*c
    if a1 == b1 and a2 < b2:
        for i in enumerate(output):
            output[i[0]] = [a1, a2 + i[0]]
    elif a1 == b1 and a2 > b2:
        for i in enumerate(output):
            output[i[0]] = [a1, a2 - i[0]]
    elif a2 == b2 and a1 < b1:
        for i in enumerate(output):
            output[i[0]] = [a1 + i[0], a2]
    elif a2 == b2 and a1 > b1:
        for i in enumerate(output):
            output[i[0]] = [a1 - i[0], a2]
    else:
        raise Exception('Input arrays must share 1 common values.')
    return output

def xy_shape(xieta, xy):
    """Takes in the local coordinate [xi,eta] and the node points [x,y]_i
      of a 8 noded isoparametric element.
      Returns the position in [x,y]"""
    xi = xieta[0]
    eta = xieta[1]
    N = np.zeros(8)
    x = xy[0]
    y = xy[1]
    N[0] = -(1 - xi)*(1 - eta)*(1 + xi + eta)/4
    N[1] = -(1 + xi)*(1 - eta)*(1 - xi + eta)/4
    N[2] = -(1 + xi)*(1 + eta)*(1 - xi - eta)/4
    N[3] = -(1 - xi)*(1 + eta)*(1 + xi - eta)/4
    N[4] = (1 - xi**2)*(1 - eta)/2
    N[5] = (1 + xi)*(1 - eta**2)/2
    N[6] = (1 - xi**2)*(1 + eta)/2
    N[7] = (1 - xi)*(1 - eta**2)/2
    return np.array([np.dot(N,x), np.dot(N,y)])

class Mesher:
    def __init__(self):
        pass
    def create(self, blocks_nums, div_nums, void = None, merge = None):
        NS = blocks_nums[0] #number of blocks in S direction
        NW = blocks_nums[1] #number of blocks in W direction
        NSW = NS*NW #total number of blocks

        NSD = div_nums[0] #array containing number of divisions in each block in S
        NWD = div_nums[1] #array containing number of divisions in each block in W

        NNS = 1 + np.sum(NSD) #total nodes along S
        NNW = 1 + np.sum(NWD) #toal nodes along W

        NNT = NNS*NNW #total number of nodes

        NNAR = np.zeros((NNS,NNW)) #initialize node list with zeros

        # Loop through and set -1 to all blocks not set to void
        BLOCKID = 0
        WPLACE = 0
        for KW in range(NW):
            SPLACE = 0
            for KS in range(NS):
                if BLOCKID in void:
                    pass
                else:
                    NNAR[SPLACE : SPLACE + NSD[KS] + 1, WPLACE : WPLACE + NWD[KW] + 1 ] = -1
                SPLACE += NSD[KS]
                BLOCKID += 1
            WPLACE += NWD[KW]
        print(NNAR)

        # Create a map of the corner points
        self.MAP = np.reshape(list(range((NS + 1)*(NW + 1))),(NS + 1, NW + 1))

        # Create an array indicating the indices of corner points of each block
        indexW = [0]
        WPLACE = 0
        indexS = [0]
        SPLACE = 0
        for k in NWD:
            WPLACE += k
            indexW.append(WPLACE)
        for p in NSD:
            SPLACE += p
            indexS.append(SPLACE)
        X, Y = np.meshgrid(indexS, indexW)
        self.POS = [None]*len(np.ravel(X))
        for i in range(len(np.ravel(X))):
            self.POS[i] = [np.ravel(X)[i], np.ravel(Y)[i]]
        # print(self.POS)

        # Set merge nodes
        if merge != None:
            for pair in merge:
                # get node numbers of nodes on merge lists
                line1 = pair[0:2]
                line2 = pair[2:4]
                # print(line1,line2)
                # print(self.POS[line1[0]], self.POS[line1[1]])
                # Convert nodes to lists of index positions to be merged
                line1 = np.array(gather_index(self.POS[line1[0]], self.POS[line1[1]]))
                # print(line1)
                # print(self.POS[line2[0]], self.POS[line2[1]])
                line2 = np.array(gather_index(self.POS[line2[0]], self.POS[line2[1]]))
                # print(line2)
                # print(np.max(line1[:,0])>np.max(line2[:,0]))
                # Find the line with the highest node number and conver that to the lower node numbers
                if np.max(line1[:,0]) < np.max(line2[:,0]):
                    for point in enumerate(line1):
                        # print(point)
                        # print(line2[point[0]])
                        node_num = np.ravel_multi_index(line2[point[0]], NNAR.shape, order="F")
                        NNAR[point[1][0], point[1][1]] = node_num
                else:
                    for point in enumerate(line2):
                        node_num = np.ravel_multi_index(line1[point[0]], NNAR.shape, order="F")
                        NNAR[point[1][0], point[1][1]] = node_num
                print(NNAR)
        # Assign final node numbers
        DUMMY = np.zeros(NNT)
        NCOUNT = 0
        for node in enumerate(np.ravel(NNAR, order="F")):
            if node[1] < 0:
                DUMMY[node[0]] = NCOUNT # Assign the node number according to sequence
                NCOUNT += 1
            elif node[1] > 0:
                DUMMY[node[0]] = node[1] # Assign the merged node number
                if node[1] == NCOUNT:
                    NCOUNT += 1
            else:
                pass
        NNAR = np.reshape(DUMMY, NNAR.shape, order="F")
        print(NNAR)



class Mesh:
    """Mesh object of a meshed region with arrays containing the nodal positions and element numbering."""
    def __init__(self):
        self.nodes = None
        self.elements = None
        self.gauss_points = None
        self.material = None

    def assign_material(self, mat_model):
        self.material = mat_model

    def make_rect_mesh(self, size, num_elements):
        """Creates a rectangular meshed region of defined by a width, height, and number of elements per side."""
        x_length = size[0]
        y_length = size[1]
        num_elements_x = num_elements[0]
        num_elements_y = num_elements[1]

        # Pre-size matrices
        num_nodes = (num_elements_x + 1)*(num_elements_y + 1)
        self.nodes = np.zeros((2, num_nodes*num_nodes), dtype='float')
        self.elements = np.zeros(((num_elements_x)*(num_elements_y), 4), dtype='int')

        #make x and y coordinates
        x_coordinates = np.linspace(0.0, x_length, num_elements_x + 1)
        y_coordinates = np.linspace(0.0, y_length, num_elements_y + 1)
        x_locations, y_locations = np.meshgrid(x_coordinates, y_coordinates, indexing='xy')

        self.nodes = np.stack((np.ndarray.flatten(x_locations), np.ndarray.flatten(y_locations)))

        k = 0 #Create list of elements with their node numbers.
        for j in range(num_elements_y):
            for i in range(num_elements_x):
                self.elements[k, :] = np.array([i+j*(num_elements_x+1), (i+1)+(j*(num_elements_x+1)), (i+1)+(j+1)*(num_elements_x+1), (i)+(j+1)*(num_elements_x+1)], dtype='int')
                k += 1
        # Create the gauss points for each element and store in self.gauss_points.  The integration rule used is 2X2.    
        GP = np.array([[-1/np.sqrt(3), -1/np.sqrt(3), 1],
                       [1/np.sqrt(3), -1/np.sqrt(3), 1], 
                       [1/np.sqrt(3), 1/np.sqrt(3), 1],
                       [-1/np.sqrt(3), 1/np.sqrt(3), 1]])
        
        N = np.zeros((8,8))
        for i in range(0,4):
            xi = GP[i,0]
            eta = GP[i,1]
            N[i*2:i*2+2,:] = [[(1/4)*(1-xi)*(1-eta), 0,  (1/4)*(1+xi)*(1-eta), 0,  (1/4)*(1+xi)*(1+eta), 0,  (1/4)*(1-xi)*(1+eta), 0], 
                               [ 0,  (1/4)*(1-xi)*(1-eta), 0,  (1/4)*(1+xi)*(1-eta), 0,  (1/4)*(1+xi)*(1+eta), 0,  (1/4)*(1-xi)*(1+eta)]]

        num_elements = self.elements.shape[0]
        Q = self.nodes[:, self.elements.ravel(order='F')]
        Q = np.reshape(Q, (8, num_elements), order='F')
        self.gps = np.reshape(N@Q, (2, num_elements*4), order='F')
        self.gauss_points = np.vstack((self.gps, np.tile(GP[:,2].T, [1, num_elements])))

        # Create the B matrix and determinate of the Jacobian for each gauss point.
        def BmatdetJ(self, x, loc):
            """Computes the B matrices and value of the Jacobian determinates for each element and guass point."""
            xi, eta = loc[0:2]
            x1, x2, x3, x4, y1, y2, y3, y4 = np.ravel(x)
            Jac = 1/4*np.array([[-(1-eta)*x1 + (1-eta)*x2 + (1+eta)*x3 - (1+eta)*x4, -(1-eta)*y1 + (1-eta)*y2 + (1+eta)*y3 - (1+eta)*y4],
                                  [-(1-xi)*x1 - (1+xi)*x2 + (1+xi)*x3 + (1-xi)*x4, -(1-xi)*y1 - (1+xi)*y2 + (1+xi)*y3 + (1-xi)*y4]])
            J11 = Jac[0,0]
            J12 = Jac[0,1]
            J21 = Jac[1,0]
            J22 = Jac[1,1]
            detJ = J11*J22-J12*J21
            A = (1.0/detJ)*np.array([[J22, -J12, 0., 0.], [0., 0., -J21, J11], [-J21, J11, J22, -J12]])
            G = (1/4)*np.array([[-(1-eta), 0, (1-eta), 0, (1+eta), 0, -(1+eta), 0],
                    [-(1-xi), 0, -(1+xi), 0, (1+xi), 0, (1-xi), 0],
                    [0, -(1-eta), 0, (1-eta), 0, (1+eta), 0, -(1+eta)],
                    [0, -(1-xi), 0, -(1+xi), 0, (1+xi), 0, (1-xi)]])
            Bmat = A@G
            return Bmat, detJ
        self.B = np.zeros((3,8,num_elements*4))
        self.detJ = np.zeros(num_elements*4)
        for elem in range(num_elements):
            for i in range(4):
                # output1, output2 = BmatdetJ(self, np.reshape(self.nodes[:, self.elements[:, elem]], (8, 1), order='F'), GP[i,:])
                output1, output2 = BmatdetJ(self, self.nodes[:, self.elements[elem]], GP[i,:])
                self.B[:,:, elem*4-(4-i)] = output1
                self.detJ[elem*4-(4-i)] = output2

    def plot(self):
        node_x_locations = self.nodes[0]
        node_y_locations = self.nodes[1]
        max_x = np.max(node_x_locations)
        max_y = np.max(node_y_locations)

        fig, ax = plt.subplots()
        for element in self.elements:
            x = self.nodes[0, element]
            y = self.nodes[1, element]
            ax.fill(x,y, facecolor='None', edgecolor='black')
        ax.scatter(node_x_locations, node_y_locations, color="black")
        for i in range(self.nodes.shape[1]):
            ax.annotate('node: '+str(i), (node_x_locations[i]+0.02*max_x, node_y_locations[i]+0.12*max_y))

            plt.arrow(node_x_locations[i], node_y_locations[i], 0.04*np.mean([max_x, max_y]), 0, width=0.01, color='green')
            ax.annotate(str(i*2), (node_x_locations[i]+0.04*max_x, node_y_locations[i] + 0.04*max_y), color='green')

            plt.arrow(node_x_locations[i], node_y_locations[i], 0, 0.04*np.mean([max_x, max_y]), width=0.01, color='orange')
            ax.annotate(str(i*2+1), (node_x_locations[i]-0.04*max_x, node_y_locations[i] + 0.04*max_y), color='orange')
        for i in range(self.elements.shape[0]):
            ax.annotate(f'Element: {i}', (np.mean(self.nodes[0, self.elements[i]]), np.mean(self.nodes[1, self.elements[i]])), ha='center', va='center', color="red")
        # plt.subplots_adjust(top=2, right=2, bottom=1.5)
        ax.set_aspect('equal', 'box')
        plt.show()

class Material_model:
    """Material model class object which contains the stiffness matrix, D, and other constitutive matrices.\n
       arguments:\n
       model_inputs = a value or array of values which populate a particular model\n
       model_type = a string which specified the type of material model.\n
       Available model input and types:\n
       [E, nu] "linear elastic"\n
       """
    def __init__(self, model_inputs: float, model_type: str):
        self.type = model_type
        if self.type == "linear elastic":
            E = model_inputs[0]
            nu = model_inputs[1]
            self.D = E/(1-nu**2)*np.array([[1, nu, 0.0], [nu, 1.0, 0.0], [0.0, 0.0, 0.5*(1.0-nu)]])

class Global_K_matrix:
    """The global stiffness matrix object."""
    def __init__(self, input_mesh: Mesh):
        self.mesh = input_mesh
        nodes = self.mesh.nodes
        elements = self.mesh.elements
        self.K_global = np.zeros((nodes.shape[1]*2, nodes.shape[1]*2))
        self.DOF_mapping = np.zeros((elements.shape[0], 8), dtype='int')
        i = 0
        for element in elements:
            self.DOF_mapping[i,:] = np.ndarray.flatten(np.array(list(map(map_DOF, element))).T, order='F')
            i += 1
    def build(self, mat_model: Material_model):
        """Constructs the global stiffness matrix."""
        D = mat_model.D
        for p in enumerate(self.mesh.elements):
            gauss_index = np.arange(4*p[0], 4*p[0]+4)
            index = self.DOF_mapping[p[0],:]
            for k in gauss_index:
                B = self.mesh.B[:,:,k]
                weight = self.mesh.gauss_points[2,k]
                detJ = self.mesh.detJ[k]
                k_element = B.T@D@B*weight*detJ
                for n in enumerate(index):
                    for m in enumerate(index):
                        self.K_global[n[1], m[1]] += k_element[n[0],m[0]]
        self.mesh.assign_material(mat_model)


class Global_T_matrix:
    """The global internal nodal force vector."""
    def __init__(self, input_mesh: Mesh):
            self.mesh = input_mesh
            nodes = self.mesh.nodes
            elements = self.mesh.elements
            self.T_global = np.zeros(nodes.shape[1]*2)
            self.DOF_mapping = np.zeros((elements.shape[0], 8), dtype='int')
            i = 0
            for element in elements:
                self.DOF_mapping[i,:] = np.ndarray.flatten(np.array(list(map(map_DOF, element))).T, order='F')
                i += 1
    def build(self, S):
        """Constructs the global internal force vector."""
        i = 0
        for p in enumerate(self.mesh.elements):
            gauss_index = np.arange(4*p[0],4*p[0]+4)
            for k in gauss_index:
                B = self.mesh.B[:,:,k]
                weight = self.mesh.gauss_points[2,k]
                detJ = self.mesh.detJ[k]
                t = B.T@S.return_all()[:,k]*weight*detJ # S is the stress field
                index = self.DOF_mapping[i,:]
                n = 0
                for position1 in index:
                    self.T_global[position1] += t[n]
                    n += 1
            i += 1

class Global_F_matrix:
    """The global applied nodal force vector."""
    def __init__(self, input_mesh: Mesh):
            self.mesh = input_mesh
            nodes = self.mesh.nodes
            elements = self.mesh.elements
            self.F_global = np.zeros(nodes.shape[1]*2)
            self.DOF_mapping = np.zeros((elements.shape[0], 8), dtype='int')
            i = 0
            for element in elements:
                self.DOF_mapping[i,:] = np.ndarray.flatten(np.array(list(map(map_DOF, element))).T, order='F')
                i += 1
    def build(self, trac_nodes, trac_value, trac_dir):
        """Constructs the global applied force vector."""
        for k in range(trac_nodes.shape[0]):
            n1 = trac_nodes[k,0]
            x1 = self.mesh.nodes[0,n1]
            y1 = self.mesh.nodes[1,n1]
            n2 = trac_nodes[k,1]
            x2 = self.mesh.nodes[0,n2]
            y2 = self.mesh.nodes[1,n2]
            len23 = np.sqrt((x2-x1)**2 + (y2-y1)**2)
            ty = 0
            tx = 0
            if trac_dir == 'x':
                tx = trac_value
            elif trac_dir == 'y':
                ty = trac_value
            self.F_global[n1*2] += tx*len23/2
            self.F_global[n1*2 + 1] += ty*len23/2
            self.F_global[n2*2] += tx*len23/2
            self.F_global[n2*2 + 1] += ty*len23/2

##### Tensor and scalar quantities #####

class Nodal_quantity:
    """An object that stores the values of nodes."""
    def __init__(self, mesh: Mesh, number_of_components: int):
        self.mesh = mesh
        self.length = mesh.nodes.shape[1]
        self.num = number_of_components
        self.values = np.zeros((self.num, self.length))
    def update(self, new_values):
        self.values = new_values

class Displacement(Nodal_quantity):
    def __init__(self, mesh: Mesh, number_of_components = 2):
        super().__init__(mesh, number_of_components)
        self.length = self.length*2
        self.values = {'U1': np.zeros(self.length), 'U2': np.zeros(self.length)}
    def update(self, new_values):
        self.values['U1'], self.values['U2'] = np.array([new_values[range(0,self.length,2)], new_values[range(1,self.length,2)]])
    def return_all(self):
        return np.ravel(np.array([self.values['U1'], self.values['U2']]), order='F')

class Elemental_quantity:
    """An object that stores the values of elements at the gauss points."""
    def __init__(self, mesh: Mesh, number_of_components: int):
        self.mesh = mesh
        self.length = mesh.gauss_points.shape[1]
        self.num = number_of_components
        self.values = np.zeros((self.num, self.length))
    def update(self, new_values):
        self.values = new_values
    def return_all(self):
        return self.values
        
class Stress(Elemental_quantity):
    """Stress tensor object."""
    def __init__(self, mesh: Mesh, number_of_components = 3):
        super().__init__(mesh, number_of_components)
        self.values = {'S11': np.zeros(self.length), 'S22': np.zeros(self.length), 'S12': np.zeros(self.length)}
        self.tensor = np.array([self.values['S11'], self.values['S22'], self.values['S12']])
    def update(self, new_values):
        self.values['S11'], self.values['S22'], self.values['S12'] = new_values
    def compute(self, U):
        newS = np.zeros((3, self.length))
        D = self.mesh.material.D
        for k in enumerate(self.mesh.elements):
            i = k[0]
            element = k[1]
            n1, n2, n3, n4 = element
            index = [n1*2, n1*2+1, n2*2, n2*2+1, n3*2, n3*2+1, n4*2, n4*2+1]
            U_element = U[index]
            gauss_index = np.arange(4*i, 4*i+4)
            for k in gauss_index:
                B = self.mesh.B[:,:,k]
                newS[:, k] = D@B@U_element
        self.values['S11'], self.values['S22'], self.values['S12'] = newS
    def return_all(self):
        return np.array([self.values['S11'], self.values['S22'], self.values['S12']])

class Strain(Elemental_quantity):
    """Strain tensor object."""
    def __init__(self, mesh: Mesh, number_of_components = 3):
        super().__init__(mesh, number_of_components)
        self.values = {'E11': np.zeros(self.length), 'E22': np.zeros(self.length), 'E12': np.zeros(self.length)}
        self.tensor = np.array([self.values['E11'], self.values['E22'], self.values['E12']])
    def update(self, new_values):
        self.values['E11'], self.values['E22'], self.values['E12'] = new_values
    def return_all(self):
        return np.array([self.values['E11'], self.values['E22'], self.values['E12']])
    
class delta_Strain(Elemental_quantity):
    """Delta strain tensor object."""
    def __init__(self, mesh: Mesh, number_of_components = 3):
        super().__init__(mesh, number_of_components)
        self.values = {'dE11': np.zeros(self.length), 'dE22': np.zeros(self.length), 'dE12': np.zeros(self.length)}
        self.tensor = np.array([self.values['dE11'], self.values['dE22'], self.values['dE12']])
    def update(self, new_values):
        self.values['dE11'], self.values['dE22'], self.values['dE12'] = new_values
    def compute(self, delU):
        delE = np.zeros((3, self.length))
        for k in enumerate(self.mesh.elements):
            i = k[0]
            element = k[1]
            n1, n2, n3, n4 = element
            index = [n1*2, n1*2+1, n2*2, n2*2+1, n3*2, n3*2+1, n4*2, n4*2+1]
            delU_element = delU[index]
            gauss_index = np.arange(4*i, 4*i+4)
            for k in gauss_index:
                B = self.mesh.B[:,:,k]
                delE[:, k] = B@delU_element
        self.values['dE11'], self.values['dE22'], self.values['dE12'] = delE     
    def return_all(self):
        return np.array([self.values['dE11'], self.values['dE22'], self.values['dE12']])

class Boundary_condition:
    """Defines a boundary condition for the model."""
    def __init__(self, DOFs: np.ndarray, values: np.ndarray, Kg: Global_K_matrix):
        self.DOFs = DOFs
        self.values = values
        K = Kg.K_global
        self.num_DOFs = len(DOFs)
        self.dim_Kglobal = K.shape[0]
        self.type = None
    # def create_LagrangeMultipliers(self):
        """Creates the C and Q matrices necessary for the Langrange multiplier approach."""
        self.type = "Lagrange multipliers"
        self.C = np.zeros((self.num_DOFs, self.dim_Kglobal))
        self.Q = np.zeros(self.num_DOFs)
        for k in range(self.num_DOFs):
            self.C[k, self.DOFs[k]] = 1.0
            self.Q[k] = self.values[k]

class Solver:
    def __init__(self, solver_type):
        self.type = solver_type

class Standard(Solver):
    def __init__(self, Kg: Global_K_matrix, Tg: Global_T_matrix, Fg: Global_F_matrix, BC: Boundary_condition, 
                 S: Stress, E: Strain, U: Displacement):
        self.K = Kg.K_global
        self.F = Fg.F_global
        self.T = Tg.T_global
        self.BC = BC
        self.S = S
        self.E = E
        self.U = U
        # self.delS = delS
        # self.delE = delE
        # self.delq = delq
        # self.dq = dq
    def start(self):
        # Construct the problem [K C';C zeros(size(C,1))]\[R;Q]
        R = self.F - self.T
        Q = self.BC.Q

        a = np.block([[self.K, self.BC.C.T],
                 [self.BC.C, np.zeros((self.BC.C.shape[0], self.BC.C.shape[0]))]])
        b = np.append(R, Q.T)

        c = np.linalg.solve(a, b)
        self.U.update(c[0:self.K.shape[0]])
        print('Solve complete.')

def plot_result(mesh: Mesh, result, component: str, U, deformed=True):
    fig, ax = plt.subplots()
    zmin = result.values[component].min()
    zmax = result.values[component].max()
    levels = MaxNLocator(nbins=10).tick_values(zmin, zmax)
    cmap = plt.colormaps['gist_rainbow_r']
    norm = BoundaryNorm(levels, ncolors=cmap.N, extend='both')
    if deformed == True:
        node_positions_x = mesh.nodes[0,:] + U.values['U1']
        node_positions_y = mesh.nodes[1,:] + U.values['U2']
    else:
        node_positions_x = mesh.nodes[0,:]
        node_positions_y = mesh.nodes[1,:]
    for element in enumerate(mesh.elements):
        x1, x2, x3, x4 = node_positions_x[element[1]]
        y1, y2, y3, y4 = node_positions_y[element[1]]
        X = [[x4, x3], [x1, x2]]
        Y = [[y4, y3], [y1, y2]]
        Z = np.array([[ np.mean(result.values[component][element[0]*4:element[0]*4+4]) ]])
        im = ax.pcolormesh(X, Y, Z, edgecolor='black', cmap=cmap, norm=norm)
    fig.colorbar(im, ticks=levels, drawedges=True, extend='both', extendrect=True)
    ax.set_title(component)
    plt.show()



# Test code
start = timeit.default_timer()

mesh1 = Mesh()
mesh1.make_rect_mesh([1,1], [2,2])
print('nodes\n', mesh1.nodes)
print('elements:\n', mesh1.elements)
steel = Material_model([30e6, 0.30], "linear elastic")
K = Global_K_matrix(mesh1)
K.build(steel)
# mesh1.plot()

E = Strain(mesh1)
dE = delta_Strain(mesh1)
S = Stress(mesh1)
U = Displacement(mesh1)
T = Global_T_matrix(mesh1)
F = Global_F_matrix(mesh1)
F.build(np.array([[7, 8]]), 150e3, 'y')
BC1 = Boundary_condition([0,1], [0, 0, 0, 0], K)
solution = Standard(K,T,F,BC1,S,E,U)
solution.start()

stop = timeit.default_timer()
print('Elapsed time: ', f'{stop-start} seconds.')

dE.compute(U.return_all())
S.compute(U.return_all())
# plot_result(mesh1, S, 'S22', U=U)
# plot_result(mesh1, S, 'S11', U=U)
# plot_result(mesh1, S, 'S12', U=U)

mesher1 = Mesher()
mesher1.create([2,2], [[2,2],[3,2]], [3], [[4,7,4,5]])
print(xy_shape([0,0], [[0,2,2,0,1,2,1,0], [0,0,4,4,0,2,4,2]]))