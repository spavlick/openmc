import os
import h5py
import numpy
import math
import matplotlib.pyplot as plt

#note - find better way to strip file name
#also order file names in order of increasing batch size
directory = '../../examples/pincell/'
filenames = [f for f in os.listdir(directory) if '.h5' in f]
h5_rm = [f.strip('.h5') for f in filenames]
entropies = numpy.zeros(max([int(f.strip('statepoint.')) for f in h5_rm]))
kl_divs = numpy.zeros(max([int(f.strip('statepoint.')) for f in h5_rm])-1)


mesh_cells = raw_input('Enter number of mesh cells in each direction as a tuple in the form (x,y,z)')
mesh_size = raw_input('Enter the dimensions of one mesh cell as a tuple in the form (x,y,z)')
mesh_center = raw_input('Enter the centerpoint of the mesh as a tuple in the form (x,y,z)')

prev_probs = None
cur_probs = None


for filename in filenames:
  f = h5py.File(directory + filename,'r')
  positions = f['source_bank']['xyz']
  num_neutrons = f['n_particles']

  #plotting the particle positions
  xvals = [e[0] for e in positions]
  yvals = [e[1] for e in positions]
  fig = plt.figure()
  plt.plot(xvals,yvals,'o')
  plt.title('Particle Positions in the XY Plane')
  plt.xlabel('X')
  plt.ylabel('Y')
  plt.grid()
  plt.show()
  fig.savefig('batch'+batch_num+'.png')

  #calculating mesh probabilities
  (num_x,num_y,num_z) = mesh_cells
  (x_width,y_width,z_width) = mesh_size
  cur_probs = numpy.zeros((mesh_cells))
  lower_left = () #find way of calculating lower left
  (last_x,last_y,last_z)=lower_left
  (cur_x,cur_y,cur_z)= lower_left
  for i in range(num_x):
    cur_x+=x_width
    for j in range(num_y):
      cur_y+=y_width
      for k in range(num_z):
        cur_z+=z_width
        counter=0
        for pos in positions:
          (x_pos,y_pos,z_pos)=pos
          if x_pos>=last_x and x_pos<cur_x and y_pos>=last_y and y_pos<cur_y and z_pos>=last_z and z_pos<cur_z:
            counter+=1
          cur_probs[i][j][k]=float(counter)/float(num_neutrons)
          last_z=cur_z
        last_y=cur_y
      last_x=cur_x

  #calculating the shannon entropies
  entropy = 0
  for i in range(num_x):
    for j in range(num_y):
      for k in range(num_z):
        entropy+=cur_probs[i][j][k]*math.log(cur_probs[i][j][k])
  entropies[index]=entropy #must initialize index

  #calculating KL Divergence
  if prev_probs != None:
    kl_div = 0
    for i in range(num_x):
      for j in range(num_y):
        for k in range(num_z):
          kl_div+=cur_probs[i][j][k]*math.log(cur_probs[i][j][k]/prev_probs[i][j][k])
    kl_divs[index]=kl_div

  prev_probs = cur_probs
  f.close()


#plotting Shannon entropies
fig = plt.figure()
plt.plot(range(len(entropies)),entropies,'o-')
plt.title('Shannon Entropy v Batch Number')
plt.xlabel('Batch Number')
plt.ylabel('Shannon Entropy')
plt.grid()
plt.show()
fig.savefig('shannon-entropies.png')

#plotting KL Divergence
fig = plt.figure()
plt.plot([e+1 for e in range(len(kl_divs))],kl_divs,'o-')
plt.title('KL Divergence v Batch Number')
plt.xlabel('Batch Number')
plt.ylabel('KL Divergence')
plt.grid()
plt.show()
fig.savefig('kl-divergences.png')
