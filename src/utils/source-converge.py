import os
import h5py
import numpy
import matplotlib.pyplot as plt

directory = '../../examples/pincell/'

mesh_cells = raw_input('Enter number of mesh cells in each direction as a tuple in the form (x,y,z)')
mesh_size = raw_input('Enter the dimensions of one mesh cell as a tuple in the form (x,y,z)')
mesh_center = mesh_size = raw_input('Enter the centerpoint of the mesh as a tuple in the form (x,y,z)')

for filename in os.listdir(directory):
  if '.h5' in filename:
    f = h5py.File(directory + filename,'r')
    positions = f['source_bank']['xyz']

    #plotting the particle positions
    xvals = [e[0] for e in positions]
    yvals = [e[1] for e in positions]
    fig = plt.figure()
    plt.plot(xvals,yvals,'o')
    plt.title('pretty.png')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid()
    plt.show()
    fig.savefig()

    #calculating mesh probabilities
    
      


    f.close()
