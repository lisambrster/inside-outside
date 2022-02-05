
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

def dist(x1,y1,x2,y2):
    d = math.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))
    return d

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
npts = 52
points = np.zeros([52,2])
# make new set of point 1/r, theta
inv_points = np.zeros([52,2])
i = 0
for iangle in range(0,360,5):
    # create vector in direction of current angle
    # vector starts at blastocoel center
    # and goes in direction cos(iangle),sin(iangle)
    iangle_rad = (iangle*np.pi)/180
    v = [0, math.cos(iangle_rad), math.sin(iangle_rad)]
    # find closest point in a nucleus in this direction
    bFound, pt = FindClosestNucleus(blastocoel_center, v, label_img)
    if (bFound):
        all_pts.append(pt)
        points[i,0] = all_pts[i][1]
        points[i,1] = all_pts[i][2]
        # get r for this point (have theta already)
        r = dist(points[i,0],points[i,1], by, bx) # r is length (origin at center)
        # create new point at 1/r, theta
        invr = 1500 - r
        theta = iangle_rad
        inv_points[i,0] = invr * math.cos(theta) + by
        inv_points[i,1] = invr  * math.sin(theta) + bx
        print('angle (deg), point, theta (rad), r ',iangle, pt[1:], theta, r)
        print('inv point ',inv_points[i,:])


        i = i + 1

# repeat the first point in the last position so it is a closed polygon
points[len(all_pts),0] = all_pts[0][1]
points[len(all_pts),1] = all_pts[0][2]
inv_points[len(all_pts),0] = all_pts[0][1]
inv_points[len(all_pts),1] = all_pts[0][2]
print(points.shape)

# plot the inverted points ...

# get the convex hull - of the inverted points
hull = ConvexHull(inv_points)

# invert points in convex hull
for ivertex in hull.vertices:
   d = dist(inv_points[ivertex, 0], inv_points[ivertex, 1], 0, 0)
   norm_x = (inv_points[ivertex,0] - by)/d
   norm_y = (inv_points[ivertex,1] - bx)/d
   theta = math.atan2( norm_y , norm_x )
   print('ivertex, r, theta, inv_x, inv_y ', ivertex, d, theta, inv_points[ivertex, 0], inv_points[ivertex, 1])
   inv_points[ivertex,0] = (1500 - d) * math.cos(theta)
   inv_points[ivertex,1] = (1500 - d) * math.sin(theta)
   print('inside hull ', ivertex, inv_points[ivertex, 0], inv_points[ivertex, 1])


# now plot inverted simplices
d,h,w = label_img.shape

for vertex in hull.vertices:
    plt.plot(inv_points[hull.vertices, 1] + bx , h - inv_points[hull.vertices, 0] - by, 'r--', lw = 2)

plt.plot(points[:,1], h - points[:,0], 'o-')
#plt.plot(points[hull.vertices[0],0], points[hull.vertices[0],1], 'ro')
plt.savefig('blastocoel_slice_11_inv_hull.jpg')
