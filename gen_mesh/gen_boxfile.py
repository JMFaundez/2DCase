import numpy as np
import matplotlib.pyplot as plt
import subprocess

#Number of elements x and y direction
Nx = 100
Ny = 20
# Upper and lower point where the airfoil is cut
xup = 0.5
xlow = 0.3
# Height to free stream
yh = 0.3
#ratio for Boundary layer
ratio = 1.8

pi = np.pi
x = np.zeros(Nx+1)

beta = np.linspace(-xlow*pi,xup*pi, Nx+1)
for i in range(len(beta)):
    if beta[i]<0:
        x[i] = -0.5*(1-np.cos(beta[i]))
    else:
        x[i] = 0.5*(1-np.cos(beta[i]))

# Write box file
f = open("airfoil_box.box","w+")
f.write("%i\n" % -2)
f.write("%i\n"% 1)
f.write("Box \n")
f.write("%3i  %3i\n"% (Nx, -Ny))
count = 0
for i in range(len(x)):
    f.write("%.16f" %  x[i])
    count = count + 1
    if count>5:
        f.write("\n")
        count = 0
    else:
        f.write(" ")
f.write("\n")
f.write("%.16f %.16f %.3f\n" % (0, yh, ratio))
f.write("P  ,P  ,W  ,W  ")
f.close()

p = subprocess.Popen(['genbox'],  stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
p_stdout = p.communicate(input='airfoil_box.box'.encode())
#p = subprocess.Popen(['genmap'],  stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
#p_stdout = p.communicate(input='box'.encode())
