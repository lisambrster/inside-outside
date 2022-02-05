

the current code reads a 3D label image that was too big for github
so you need to change it to read the 2D image I put here


MeasureCavityVolume.py - just finds the points and the convex hull

InvoluteCavityVolume.py - I invert the points, find the convex hull and then re-invert them back again..seems to work now!

the bottom right should be included - I tried adding more points using more angles but it didn't help - because that part is tangent to the vectors and so few points land on the sides of the part that sticks out
