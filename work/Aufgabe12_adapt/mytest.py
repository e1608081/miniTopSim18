import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt

x1=np.array([1.,3.,6.,7.,10.,16.])
y1=np.array([0.,0.,0.,0.,0.,0.])

MAX_SEGLEN = 2.

x1=np.linspace(0,1,6)
y1=7*np.sin(2*np.pi*x1)

xdist = np.abs(x1[1:]-x1[:-1])
ydist = np.abs(y1[1:]-y1[:-1])
#prepend zero because insert adds before given index
dist = np.insert(np.sqrt(xdist**2 + ydist**2), 0, 0)
#dist_mask = dist>MAX_SEGLEN
dist_index_mask, = np.where(dist > MAX_SEGLEN)


num_insert_nodes=np.int32(np.ceil(dist[dist_index_mask]/MAX_SEGLEN))
#num_insert_nodes=np.floor_divide(dist[dist_index_mask],MAX_SEGLEN)

newpoints = np.split(x1, dist_index_mask)
#newpoints = np.zeros(np.sum(num_insert_nodes)-len(num_insert_nodes))

for i in range(len(dist_index_mask)):
    line = np.linspace(x1[dist_index_mask[i]-1],
                                 x1[dist_index_mask[i]],
                                 num_insert_nodes[i],
                                 endpoint = False)
    line=np.delete(line, 0)
    newpoints[i] = np.append(newpoints[i],line)
newarray = np.concatenate(newpoints)
        #np.linspace(x1[dist_index_mask-1],x1[dist_index_mask],num_insert_nodes)

x2=newarray

f2=scipy.interpolate.interp1d(x1,y1,'quadratic')
f3=scipy.interpolate.interp1d(x1,y1,'cubic')
y2=f2(x2)
y3=f3(x2)
plt.plot(x1, y1, '*g', x2, y2, '1r')
plt.show()