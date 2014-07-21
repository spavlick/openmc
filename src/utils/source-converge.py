import os
import h5py
import numpy
import math
import matplotlib.pyplot as plt
from ast import literal_eval

#get batch numbers from file names
directory = '../../examples/pincell/'
filenames = [f for f in os.listdir(directory) if '.h5' in f]
h5_rm = [f[:-3] for f in filenames]
batch_nums = [int(f[11:]) for f in h5_rm]
batch_nums = sorted(batch_nums,key=int) #sort in ascending order
ord_filenames = []
for i in range(len(batch_nums)):
  ord_filenames.append('statepoint.' + str(batch_nums[i]) + '.h5')

#make empty lists to store entropies and kl divergences
entropies = numpy.zeros(len(batch_nums))
kl_divs = numpy.zeros(len(batch_nums)-1)


mesh_cells = (4,4,1)
mesh_size = (1.25984,1.25984,100000000000.0)
mesh_center = (0,0,0)


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
    cur_y = y_left
    for j in range(num_y):
      cur_y+=y_width
      cur_z = z_left
      for k in range(num_z):
        cur_z+=z_width
        for pos in positions:
          (x_pos,y_pos,z_pos)=pos
          if x_pos>=last_x and x_pos<cur_x and y_pos>=last_y and y_pos<cur_y and z_pos>=last_z and z_pos<cur_z:
            cur_probs[i][j][k]+=1
        last_z=cur_z
      last_y=cur_y
      last_z=z_left
    last_x=cur_x
    last_y=y_left
  cur_probs[:,:,:] /= float(num_neutrons)

  #calculating the shannon entropies
  entropy=cur_probs[:,:,:]*numpy.log(cur_probs[:,:,:])
  entropy=numpy.nan_to_num(entropy)
  entropies[index] = -entropy.sum()

  #calculating KL Divergence
  if prev_probs != None:
    kl_div=cur_probs[:,:,:]*numpy.log(cur_probs[:,:,:]/prev_probs[:,:,:])
    kl_div=numpy.nan_to_num(kl_div)
    kl_divs[index-1]=kl_div.sum()

  prev_probs = cur_probs
  f.close()

#plotting Shannon entropies
fig = plt.figure()
plt.plot(batch_nums,entropies,'o-')
plt.title('Shannon Entropy v Batch Number')
plt.xlabel('Batch Number')
plt.ylabel('Shannon Entropy')
plt.grid()
fig.savefig('shannon-entropies.png')

#plotting KL Divergence
fig = plt.figure()
kl_batches = numpy.copy(batch_nums)
kl_batches = numpy.delete(kl_batches,numpy.array([0]))
plt.plot(kl_batches,kl_divs,'o-')
plt.title('KL Divergence v Batch Number')
plt.xlabel('Batch Number')
plt.ylabel('KL Divergence')
plt.grid()
fig.savefig('kl-divergences.png')
