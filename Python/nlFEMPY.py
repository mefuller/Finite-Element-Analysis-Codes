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

        k = 0 #Create list of elements with their node numbers.
        for j in range(num_elements_y):
            for i in range(num_elements_x):
                self.elements[:, k] = np.array([i+j*(num_elements_x+1), (i+1)+(j*(num_elements_x+1)), (i+1)+(j+1)*(num_elements_x+1), (i)+(j+1)*(num_elements_x+1) ]).T
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

        num_elements = self.elements.shape[1]
        Q = self.nodes[:, self.elements.ravel(order='F')]
        Q = np.reshape(Q, (8, num_elements), order='F')
        self.gps = np.reshape(N@Q, (2, num_elements*4), order='F')
        print(GP[:,2].T)
        print('output', np.tile(GP[:,2].T, [1, num_elements]))
        self.gauss_points = np.vstack((self.gps, np.tile(GP[:,2].T, [1, num_elements])))

        # Create the B matrix and determinate of the Jacobian for each gauss point.
        # def BmatdetJ(self, x, loc):
        #     xi = loc[0]; eta = loc[1]
        #     x1 = x[0]
        #     y1 = x[1]
        #     x2 = x[2]
        #     y2 = x[3]
        #     x3 = x[4]
        #     y3 = x[5]
        #     x4 = x[6]
        #     y4 = x[7]
        #     Jac = (1/4)*np.array([[-(1-eta)*x1 + (1-eta)*x2 + (1+eta)*x3 - (1+eta)*x4, -(1-eta)*y1 + (1-eta)*y2 + (1+eta)*y3 - (1+eta)*y4],
        #                   -(1-xi)*x1 - (1+xi)*x2 + (1+xi)*x3 + (1-xi)*x4, -(1-xi)*y1 - (1+xi)*y2 + (1+xi)*y3 + (1-xi)*y4]])
        #     J11 = Jac[0,0]; J12 = Jac[0,1]; J21 = Jac[1,0]; J22 = Jac[1,1]
        #     detJ = J11*J22-J12*J21
        #     A = (1/detJ)*np.array([[J22, -J12, 0, 0], [0, 0, -J21, J11], [-J21, J11, J22, -J12]])
        #     G = (1/4)*np.array([[-(1-eta), 0, (1-eta), 0, (1+eta), 0, -(1+eta), 0],
        #             [-(1-xi), 0, -(1+xi), 0, (1+xi), 0, (1-xi), 0],
        #             [0, -(1-eta), 0, (1-eta), 0, (1+eta), 0, -(1+eta)],
        #             [0, -(1-xi), 0, -(1+xi), 0, (1+xi), 0, (1-xi)]])
        #     Bmat = A@G
        # for elem in range(num_elements):
        
        #     for i in range(4):
        #         [B(:,:,elem*Iorder^2-(Iorder^2-i)), detJ(elem*Iorder^2-(Iorder^2-i))]= BmatdetJ(reshape(nodes(:,elements(:,elem)),[8 1]),GP{Iorder}(i,:));
   
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