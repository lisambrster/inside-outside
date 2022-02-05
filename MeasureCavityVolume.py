
# for certain embryo
# given point in blastocoel
# measure volume of cavity
# from point - find alll closest points in all directions
# uses these points to constuct convex hull
# measure its volume

import numpy as np
import os
import math
from scipy.ndimage.morphology import distance_transform_edt
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import matplotlib.pyplot as plt

def  FindClosestNucleus(ctr, v, label_img):
    bFound = False
    pt = []
    d,h,w = label_img.shape
    # for next point along vector v - check if it is a nucleus
    for iTest in range(2000):
        vx = int(v[0]*iTest + ctr[0])
        vy = int(v[1]*iTest + ctr[1])
        vz = int(v[2]*iTest + ctr[2])
        if (vx < 0) or (vy < 0) or (vz < 0):
            return bFound, pt
        if (vx >= d) or (vy >= h) or (vz >= w):
            return bFound, pt
        if (label_img[vx,vy,vz] > 0):
            bFound = True
            pt = [vx,vy,vz]
            return bFound, pt
    return bFound, pt

embryo_path = '/mnt/ceph/users/lbrown/Labels3DMouse/GTSets/2022_Full/'
embryo_file = '/mnt/ceph/users/lbrown/Labels3DMouse/GTSets/2022_Full//F55_185/F55_185/masks/F55_185_masks_0001.npy'
#embryo_file = '/mnt/ceph/users/lbrown/Labels3DMouse/GTSets/2022_Full/F46_110/F46_110/masks/F46_110_masks_0001.npy'
label_img = np.load(embryo_file)
print(label_img.shape)
labels = np.unique(label_img)
nlabels = labels.shape[0]
print('Number of Labels ',nlabels-1)  # 102 (why isn't it 104 - is metaphase counted in interphase or not included??)

# point inside blastocoel
blastocoel_center = [15,693.471]  # z,y,x - best estimate over all slices
blastocoel_center = [11,699,501] # good slice for 2D
bx = blastocoel_center[2]
by = blastocoel_center[1]

# do in 2D first.. (slice 11 only)
all_pts = []
for iangle in range(0,360,5):
    # create vector in direction of current angle
    # vector starts at blastocoel center
    # and goes in direction cos(iangle),sin(iangle)
    iangle_rad = (iangle*np.pi)/180
    v = [0, math.cos(iangle_rad), math.sin(iangle_rad)]
    # find closest point in a nucleus in this direction
    bFound, pt = FindClosestNucleus(blastocoel_center, v, label_img)
    if (bFound):
        print('for this angle, found point ',iangle, pt)
        all_pts.append(pt)


points = np.zeros([len(all_pts)+1,2])
for i in range(len(all_pts)):
    points[i,0] = all_pts[i][1]
    points[i,1] = all_pts[i][2]
points[len(all_pts),0] = all_pts[i][1]
points[len(all_pts),1] = all_pts[i][2]
print(points.shape)

# get the convex hull
#rng = np.random.default_rng()
#points = rng.random((30, 2))   # 30 random points in 2-D
hull = ConvexHull(points)
d,h,w = label_img.shape

plt.plot(points[:,1], h - points[:,0], 'o-')
for simplex in hull.simplices:
   plt.plot(points[simplex, 1], h - points[simplex, 0], 'k-')
#plt.plot(points[hull.vertices[0],0], points[hull.vertices[0],1], 'ro')
plt.savefig('blastocoel_slice_11.jpg')
