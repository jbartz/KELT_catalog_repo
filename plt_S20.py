#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

colors = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854']


#e_data = np.genfromtxt('cat_S20_east.dat')
#w_data = np.genfromtxt('cat_S20_west.dat')

S20_data = np.genfromtxt('./canon/fields_/S20.dat')

S20_ra = S20_data[:,2]
S20_dec = S20_data[:,3]
S20_x = S20_data[:,0]
S20_y = S20_data[:,1]


S25_data = np.genfromtxt('./canon/fields_/S25.dat')

S25_ra = S25_data[:,2]
S25_dec = S25_data[:,3]
S25_x = S25_data[:,0]
S25_y = S25_data[:,1]

S26_data = np.genfromtxt('./canon/fields_/S26.dat')

S26_ra = S26_data[:,3]
S26_dec = S26_data[:,4]
S26_x = S26_data[:,0]
S26_y = S26_data[:,1]

S24_data = np.genfromtxt('./canon/fields_/S24.dat')

S24_ra = S24_data[:,2]
S24_dec = S24_data[:,3]
S24_x = S24_data[:,0]
S24_y = S24_data[:,1]

S27_data = np.genfromtxt('./canon/fields_/S27.dat')

S27_ra = S27_data[:,2]
S27_dec = S27_data[:,3]
S27_x = S27_data[:,0]
S27_y = S27_data[:,1]

S28_data = np.genfromtxt('./canon/fields_/S28.dat')

S28_ra = S28_data[:,2]
S28_dec = S28_data[:,3]
S28_x = S28_data[:,0]
S28_y = S28_data[:,1]

S19_data = np.genfromtxt('./canon/fields_/S19.dat')

S19_ra = S19_data[:,2]
S19_dec = S19_data[:,3]
S19_x = S19_data[:,0]
S19_y = S19_data[:,1]
#w_ra = w_data[:,5]
#w_dec = w_data[:,6]
#w_x = w_data[:,0]
#w_y = w_data[:,1]


fig = plt.figure()
fig.set_size_inches(10.0,10.0)
plt.scatter(S20_ra, S20_dec,c='r',marker='o',s=0.55, label='S20')
plt.scatter(S19_ra, S19_dec,c='b',marker='o',s=0.55, alpha = 0.5, label='S19')

#plt.scatter(S25_ra, S25_dec,c='k',marker='o',s=0.55, alpha = 0.5)
#plt.scatter(S26_ra, S26_dec,c=colors[1],marker='o',s=0.55, alpha = 0.5, label='S26')
#plt.scatter(S24_ra, S24_dec,c=colors[0],marker='o',s=0.55, alpha = 0.5, label='S24')
#plt.scatter(S27_ra, S27_dec,c=colors[2],marker='o',s=0.55, alpha = 0.5, label='S27')
#plt.scatter(S28_ra, S28_dec,c=colors[3],marker='o',s=0.55, alpha = 0.5, label='S28')
#plt.scatter(w_ra, w_dec,c='b', marker=',',s=2.39, alpha=0.5)
#plt.show()
plt.gca().invert_xaxis()
plt.xlabel('RA')
plt.ylabel('Dec')
#plt.show()

plt.title('S19 and S20 canon')

plt.gca().legend()
for label in plt.gca().legend().get_lines():
    label.set_markersize(10.0)
plt.savefig('/home/jbartz/Desktop/S19_S20_canon.png')
plt.close()