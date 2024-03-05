import numpy as np

class Mesh:
    """Creates a Mesh object of a rectangular meshed region with arrays containing the nodal positions and element numbering."""
    def __init__(self):
        self.nodes = []
        self.elements = []
        self.gauss_points = []

    def make_rect_mesh(self, size, num_elements):
        x_length = size[0]
        y_length = size[1]
        num_elements_x = num_elements[0]
        num_elements_y = num_elements[1]

        # Pre-size matrices
        num_nodes = (num_elements_x + 1)*(num_elements_y + 1)
        self.nodes = np.zeros((2, num_nodes*num_nodes), dtype='float')
        self.elements = np.zeros((4, (num_elements_x)*(num_elements_y)), dtype='int')

        #make x and y coordinates
        x_coordinates = np.linspace(0.0, x_length, num_elements_x + 1)
        y_coordinates = np.linspace(0.0, y_length, num_elements_y + 1)
        x_locations, y_locations = np.meshgrid(x_coordinates, y_coordinates, indexing='xy')

        self.nodes = np.stack((np.ndarray.flatten(x_locations), np.ndarray.flatten(y_locations)))

        k = 0 #counter for element matrix loop
        for j in range(num_elements_y):
            for i in range(num_elements_x):
                self.elements[:, k] = np.array([i+j*(num_elements_x+1), (i+1)+(j*(num_elements_x+1)), (i+1)+(j+1)*(num_elements_x+1), (i)+(j+1)*(num_elements_x+1) ]).T
                k += 1

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
            
        num_elements = self.elements.shape[1]
        Q = self.nodes[:, self.elements.ravel(order='F')]
        Q = np.reshape(Q, (8, num_elements), order='F')
        self.gps = np.reshape(N@Q, (2, num_elements*4), order='F')
        print(GP[:,2].T)
        print('output', np.tile(GP[:,2].T, [1, num_elements]))
        self.gauss_points = np.vstack((self.gps, np.tile(GP[:,2].T, [1, num_elements])))

class Material_model:
    """Creates the model's constitutive materices."""
    def __init__(self, model_inputs, model_type):
        self.type = model_type

        if self.type == "linear elastic":
            E = model_inputs[0]
            nu = model_inputs[1]
            self.D_matrix = E/(1-nu**2)*np.array([[1, nu, 0.0], [nu, 1.0, 0.0], [0.0, 0.0, 0.5*(1.0-nu)]])

class Global_K_matrix:
    def __init__(self, nodes, elements, material_model):
        self.K_global = np.zeros(nodes.shape[1]*2, nodes.shape[1]*2)



mesh1 = Mesh()
mesh1.make_rect_mesh([10,2], [10,2])
print('nodes\n', mesh1.nodes)
print('elements:\n', mesh1.elements)