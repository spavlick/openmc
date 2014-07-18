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
batch_nums = [int(f.strip('statepoint.')) for f in h5_rm]
batch_nums = sorted(batch_nums,key=int) #sort in ascending order
entropies = numpy.zeros(max(batch_nums))
kl_divs = numpy.zeros(max(batch_nums)-1)
ord_filenames = numpy.zeros(max(batch_nums))
for i in range(len(ord_filenames)):
  ord_filenames[i] = 'statepoint.' + batch_nums[i] + '.h5'



mesh_cells = raw_input('Enter number of mesh cells in each direction as a tuple in the form (x,y,z)')
mesh_size = raw_input('Enter the dimensions of the mesh as a tuple in the form (x,y,z)')
mesh_center = raw_input('Enter the centerpoint of the mesh as a tuple in the form (x,y,z)')

#calculate dimensions for a single cell
(num_x,num_y,num_z) = mesh_cells
(whole_x_width,whole_y_width,whole_z_width) = mesh_size
x_width = whole_x_width/float(num_x)
y_width = whole_y_width/float(num_y)
z_width = whole_z_width/float(num_z)

#find lower left of mesh
(x_center,y_center,z_center)=mesh_center
x_left = x_center-((num_x/2)*x_width)
y_left = y_center-((num_y/2)*y_width)
z_left = z_center-((num_z/2)*z_width)
lower_left = (x_left,y_left,z_left)

prev_probs = None
cur_probs = None


for index, filename in enumerate(ord_filenames):
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
  cur_probs = numpy.zeros((mesh_cells))
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
