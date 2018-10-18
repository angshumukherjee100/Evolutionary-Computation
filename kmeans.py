import random
import numpy as np
from numpy import vstack,array
from numpy.random import rand
from scipy.cluster.vq import kmeans,vq

arr = np.array([[1.0,2.0,5.3,2.8],[3.0,4.0,1.9,3.6],[6.0,8.0,9.3,2.9],[7.0,3.0,4.8,7.3],[5.0,1.0,8.2,6.4],[2.0,6.0,1.2,4.5],[9.0,1.0,7.0,2.3],[4.0,7.0,8.2,4.6]])
centroids,distortion = kmeans(arr,2)

idx, distortion = vq(arr,centroids)

bestpos = 4

 



print(arr[idx == 0])
print("===============")
print(arr[idx == 1])

'''
plot(data[idx==0,0],data[idx==0,1],'ob',
     data[idx==1,0],data[idx==1,1],'or')
plot(centroids[:,0],centroids[:,1],'sg',markersize=8)
show()
plot(x[idx == 0,0],'ob', x[idx == 1, 0], 'or')
plot(centroids[:,0],'sg,markersize = 8')

'''	


