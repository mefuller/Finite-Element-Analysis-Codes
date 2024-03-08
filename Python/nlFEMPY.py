import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as ptch


def map_DOF(x):
    return x*2, x*2+1

class Mesh:
    """Mesh object of a meshed region with arrays containing the nodal positions and element numbering."""
    def __init__(self):
        self.nodes = None
        self.elements = None
        self.gauss_points = None

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
            J = 1/4*np.array([[-(1-eta)*x1 + (1-eta)*x2 + (1+eta)*x3 - (1+eta)*x4, -(1-eta)*y1 + (1-eta)*y2 + (1+eta)*y3 - (1+eta)*y4],
                                  [-(1-xi)*x1 - (1+xi)*x2 + (1+xi)*x3 + (1-xi)*x4, -(1-xi)*y1 - (1+xi)*y2 + (1+xi)*y3 + (1-xi)*y4]])
            J11 = J[0,0]
            J12 = J[0,1]
            J21 = J[1,0]
            J22 = J[1,1]
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
    def make_plot(self):
        node_x_locations = self.nodes[0]
        node_y_locations = self.nodes[1]
        max_x = np.max(node_x_locations)
        max_y = np.max(node_y_locations)

        fig, ax = plt.subplots()
        for element in self.elements:
            x = self.nodes[0, element]
            y = self.nodes[1, element]
            ax.fill(x,y, facecolor='None', edgecolor='black')
        ax.scatter(node_x_locations, node_y_locations, color="blue")
        for i in range(self.nodes.shape[1]):
            ax.annotate(str(i), (node_x_locations[i]+0.02*max_x, node_y_locations[i]+0.02*max_y))
        for i in range(self.elements.shape[0]):
            ax.annotate(f'Element: {i}', (np.mean(self.nodes[0,self.elements[i]]), np.mean(self.nodes[1,self.elements[i]])), ha='center', va='center', color="red")
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
    def __init__(self, mesh: Mesh):
        nodes = mesh.nodes
        elements = mesh.elements
        self.K_global = np.zeros((nodes.shape[1]*2, nodes.shape[1]*2))
        self.DOF_mapping = np.zeros((elements.shape[0], 8), dtype='int')
        i = 0
        for element in elements:
            self.DOF_mapping[i,:] = np.ndarray.flatten(np.array(list(map(map_DOF, element))).T, order='F')
            i += 1
    def build(self, mesh: Mesh, material_model: Material_model):
        """Constructs the global stiffness matrix."""
        D = material_model.D
        i = 0
        for element in mesh.elements:
            for k in range(4*i,4*i+4):
                B = mesh.B[:,:,k]
                weight = mesh.gauss_points[2,k]
                detJ = mesh.detJ[k]
                k = B.T@D@B*weight*detJ
            # Scatter k to Kg
            index = self.DOF_mapping[i,:]
            n = 0
            for position1 in index:
                m = 0
                for position2 in index:
                    self.K_global[position1, position2] = k[n,m]
                    m += 1
                n += 1
            i += 1

class Global_T_matrix:
    """The global internal nodal force vector."""
    def __init__(self, mesh: Mesh):
            nodes = mesh.nodes
            elements = mesh.elements
            self.T_global = np.zeros((nodes.shape[1]*2, 1))
            self.DOF_mapping = np.zeros((elements.shape[0], 8), dtype='int')
            i = 0
            for element in elements:
                self.DOF_mapping[i,:] = np.ndarray.flatten(np.array(list(map(map_DOF, element))).T, order='F')
                i += 1
    def build(self, mesh: Mesh, S):
        """Constructs the global internal force vector."""
        i = 0
        for element in mesh.elements.T:
            for k in range(4*i,4*i+4):
                B = mesh.B[:,:,k]
                weight = mesh.gauss_points[2,k]
                detJ = mesh.detJ[k]
                t = B.T@S[:,k]*weight*detJ # S is the stress field
            # Scatter t to Tg !!!! Check this !!!!!
            index = self.DOF_mapping[i,:]
            n = 0
            for position1 in index:
                self.T_global[position1, 0] = t[n,0]
                n += 1
            i += 1

class Global_F_matrix:
    """The global applied nodal force vector."""
    def __init__(self, mesh: Mesh):
            nodes = mesh.nodes
            elements = mesh.elements
            self.F_global = np.zeros((nodes.shape[1]*2, 1))
            self.DOF_mapping = np.zeros((elements.shape[0], 8), dtype='int')
            i = 0
            for element in elements:
                self.DOF_mapping[i,:] = np.ndarray.flatten(np.array(list(map(map_DOF, element))).T, order='F')
                i += 1
    def build(self, mesh: Mesh):
        """Constructs the global applied force vector."""

class Nodal_quantity:
    """An object that stores the values of nodes."""
    def __init__(self, mesh: Mesh, number_of_components: int):
        length = mesh.nodes.shape[1]
        num = number_of_components
        self.values = np.zeros((num, length))
    def update(self, values):
        self.values = values

class Elemental_quantity:
    """An object that stores the values of elements at the gauss points."""
    def __init__(self, mesh: Mesh, number_of_components: int):
        length = mesh.gauss_points.shape[1]
        num = number_of_components
        self.values = np.zeros((num, length))
    def update(self, values):
        self.values = values

class Boundary_condition:
    """Defines a boundary condition for the model."""
    def __init__(self, DOFs: np.ndarray, values: np.ndarray, Kg: Global_K_matrix):
        self.DOFs = DOFs
        self.values = values
        K = Kg.K_global
        self.num_DOFs = DOFs.shape[0]
        self.dim_Kglobal = K.shape[0]
        self.type = None
    def create_LagrangeMultipliers(self):
        """Creates the C and Q matrices necessary for the Langrange multiplier approach."""
        self.type = "Lagrange multipliers"
        self.C = np.zeros((self.num_DOFs, self.dim_Kglobal))
        self.Q = np.zeros((self.num_DOFs, 1))
        for k in range(self.num_DOFs):
            self.C[k, self.DOFs[k]] = 1.0
            self.Q[k, 0] = self.values[k]
        # do dadlamb = [K C';C zeros(size(C,1))]\[R;Q] in the solver to apply the BC.


class Solver:
    def __init__(self):
        self.output = None

## Plotting Functions ##
        
def plot(mesh: Mesh):
    node_x_locations = mesh.nodes[0]
    node_y_locations = mesh.nodes[1]
    max_x = np.max(node_x_locations)
    max_y = np.max(node_y_locations)

    fig, ax = plt.subplots()
    for element in mesh.elements:
        x = mesh.nodes[0, element]
        y = mesh.nodes[1, element]
        ax.fill(x,y, facecolor='None', edgecolor='black')
    ax.scatter(node_x_locations, node_y_locations, color="blue")
    for i in range(mesh.nodes.shape[1]):
        ax.annotate(str(i), (node_x_locations[i]+0.02*max_x, node_y_locations[i]+0.02*max_y))
    for i in range(mesh.elements.shape[0]):
        ax.annotate(f'Element: {i}', (np.mean(mesh.nodes[0,mesh.elements[i]]), np.mean(mesh.nodes[1,mesh.elements[i]])), ha='center', va='center', color="red")
    plt.show()


# Test code
mesh1 = Mesh()
mesh1.make_rect_mesh([2,1], [1,1])
print('nodes\n', mesh1.nodes)
print('elements:\n', mesh1.elements)
steel = Material_model([30e6, 0.3], "linear elastic")
print('D:\n', steel.D)
Kg = Global_K_matrix(mesh1)
print('Kg.DOF_mapping:\n',Kg.DOF_mapping)
Kg.build(mesh1, steel)
print('B:', mesh1.B)