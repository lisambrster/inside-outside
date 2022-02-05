

the current code reads a 3D label image that was too big for github
so you need to change it to read the 2D image I put here

there's some bug w.r.t the order of operations and origin ..

MeasureCavityVolume.py - just finds the points and the convex hull

InvoluteCavityVolume.py - was the attempt to invert the points, find the convex hull and then invert them back again..
