import os
import h5py
import numpy
import math
import matplotlib.pyplot as plt
from ast import literal_eval

#get batch numbers from file names
directory = '../../examples/pincell/'
filenames = [f for f in os.listdir(directory) if '.h5' in f]
h5_rm = [f.strip('.h5') for f in filenames]
batch_nums = [int(f.strip('statepoint.')) for f in h5_rm]
batch_nums = sorted(batch_nums,key=int) #sort in ascending order
ord_filenames = [] #numpy.empty(len(filenames),dtype=str)
for i in range(len(batch_nums)):
  ord_filenames.append('statepoint.' + str(batch_nums[i]) + '.h5') #[i] = 'statepoint.' + str(batch_nums[i]) + '.h5'

#make empty lists to store entropies and kl divergences
entropies = numpy.zeros(max(batch_nums))
kl_divs = numpy.zeros(max(batch_nums)-1)


mesh_cells = literal_eval(raw_input('Enter number of mesh cells in each direction as a tuple in the form (x,y,z)'))
mesh_size = literal_eval(raw_input('Enter the dimensions of the mesh as a tuple in the form (x,y,z)'))
mesh_center = literal_eval(raw_input('Enter the centerpoint of the mesh as a tuple in the form (x,y,z)'))


#calculate dimensions for a single cell
(num_x,num_y,num_z) = mesh_cells
(whole_x_width,whole_y_width,whole_z_width) = mesh_size
x_width = whole_x_width/float(num_x)
y_width = whole_y_width/float(num_y)
z_width = whole_z_width/float(num_z)

#find lower left of mesh
(x_center,y_center,z_center)=mesh_center
x_left = x_center-((num_x/2.0)*x_width)
y_left = y_center-((num_y/2.0)*y_width)
z_left = z_center-((num_z/2.0)*z_width)
lower_left = (x_left,y_left,z_left)
print lower_left

#create variables for mesh probabilities
prev_probs = None
cur_probs = None


for index, filename in enumerate(ord_filenames):
  f = h5py.File(directory + filename,'r')
  positions = f['source_bank']['xyz']
  num_neutrons = int(f['n_particles'][0])

  #plotting the particle positions
  xvals = [e[0] for e in positions]
  yvals = [e[1] for e in positions]
  fig = plt.figure()
  plt.plot(xvals,yvals,'o')
  plt.title('Particle Positions in the XY Plane')
  plt.xlabel('X')
  plt.ylabel('Y')
  plt.grid()
  fig.savefig('batch'+str(batch_nums[index])+'.png')

  #calculating mesh probabilities
  cur_probs = numpy.zeros(mesh_cells)
  (last_x,last_y,last_z)=lower_left
  (cur_x,cur_y,cur_z)= lower_left
  for i in range(num_x):
    cur_x+=x_width
    cur_y=y_left
    cur_z=z_left
    for j in range(num_y):
      cur_y+=y_width
      cur_z=z_left
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
        if not cur_probs[i][j][k] == 0.0:
          entropy+=cur_probs[i][j][k]*math.log(cur_probs[i][j][k])
  entropies[index]=entropy 

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
plt.plot(batch_nums,entropies,'o-')
plt.title('Shannon Entropy v Batch Number')
plt.xlabel('Batch Number')
plt.ylabel('Shannon Entropy')
plt.grid()
plt.show()
fig.savefig('shannon-entropies.png')

#plotting KL Divergence
fig = plt.figure()
kl_batches = copy(batch_nums)
kl_batches.pop(0)
plt.plot(kl_batches,kl_divs,'o-')
plt.title('KL Divergence v Batch Number')
plt.xlabel('Batch Number')
plt.ylabel('KL Divergence')
plt.grid()
plt.show()
fig.savefig('kl-divergences.png')
