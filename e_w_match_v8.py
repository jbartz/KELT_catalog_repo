#!/usr/bin/python
import os,sys
import numpy as np
import time
import matplotlib.pyplot as plt
from math import *
from os import listdir
from os.path import isfile, join


#for output file, if no matches, just fill in the rest of columns with dummy index, use 'NULL'
start_time = time.clock() #used to time the code
#Kelt columns are: x, y, ID, mag, mag_err, RA, Dec

###-----This section from Josh Pepper's code-----###
def rAngSep(ra1r, dec1r, ra2r, dec2r):
   """
   Compute angular separation(s) with a dot product.  All input/output is
   in radians.  Inputs are converted to Cartesian coordinates and their
   dot product is computed.  The arccosine of the dot product is the
   angular separation (since A dot B = |A||B| * cos(angular_separation).
   """
   x1 = np.cos(dec1r) * np.cos(ra1r)
   y1 = np.cos(dec1r) * np.sin(ra1r)
   z1 = np.sin(dec1r)
   x2 = np.cos(dec2r) * np.cos(ra2r)
   y2 = np.cos(dec2r) * np.sin(ra2r)
   z2 = np.sin(dec2r)
   dot = x1*x2 + y1*y2 + z1*z2
   #print dot
   return np.arccos(dot)

## Angular separation in degrees (a wrapper for the above):
def dAngSep(ra1d, dec1d, ra2d, dec2d):
   """
   Compute angular separation(s) using a dot product.  This is a wrapper
   for the rAngSep() function.  See its docstring for more info.
   """
   ra1r, dec1r = np.radians([ra1d, dec1d])
   ra2r, dec2r = np.radians([ra2d, dec2d])
   return np.degrees(rAngSep(ra1r, dec1r, ra2r, dec2r))
###-----End copied section for computing angular seperation-----###


kelt_cat_files = [ f for f in listdir('/media/sf_Astro/k_cat/kelt_cat_raw/') if isfile(join('/media/sf_Astro/k_cat/kelt_cat_raw/',f))]
for file_no in range(0,len(kelt_cat_files),2):
#for file_no in range(len(kelt_cat_files[36:]),len(kelt_cat_files),2):
#for file_no in range(0,2,2):
  east_file = kelt_cat_files[file_no]
  west_file = kelt_cat_files[file_no+1]
  field = east_file[4:7]
  #f = open('/media/sf_Astro/k_cat/kelt_cat_raw/{0}'.format(west_file),mode='r')
  #file_w_line = f.readlines() #is a list of lines from west_file
  #ar_file_w_line = np.array(file_w_line) #is an array of lines from west_file
  ##print 'read all lines in file 2...'
  #f.close()  
  east_data = np.genfromtxt('/media/sf_Astro/k_cat/kelt_cat_raw/{0}'.format(east_file), dtype='|S44')
  west_data = np.genfromtxt('/media/sf_Astro/k_cat/kelt_cat_raw/{0}'.format(west_file), dtype='|S44')
  #tacks on east or west, and field number: lc_xxxxx.dat -> e_02_lc_xxxxx.dat. Takes very little time to run this part.
  for i in range(east_data.shape[0]):
      east_data[i][2] = "e_{0}_".format(field) + east_data[i][2]
  for i in range(west_data.shape[0]):
      west_data[i][2] = "w_{0}_".format(field) + west_data[i][2]
  east_sort = east_data[east_data[:,6].argsort()] #sorts by Dec
  west_sort = west_data[west_data[:,6].argsort()] #sorts by Dec

#creates 1d arrays (dtype=float) from columns in the data files
#These must be arrays (not lists) for the matching conditions below
  w_id = np.array([ row[2] for row in west_sort ])
  e_id = np.array([ row[2] for row in east_sort ])
  w_dec = np.array([ row[6] for row in west_sort ]).astype(np.float)
  e_dec = np.array([ row[6] for row in east_sort ]).astype(np.float)
  w_ra = np.array([ row[5] for row in west_sort ]).astype(np.float)
  e_ra = np.array([ row[5] for row in east_sort ]).astype(np.float)
  w_mag = np.array([ row[3] for row in west_sort ]).astype(np.float)
  e_mag = np.array([ row[3] for row in east_sort ]).astype(np.float)
  e_x = np.array([ row[0] for row in east_sort ]).astype(np.float)
  e_y = np.array([ row[1] for row in east_sort ]).astype(np.float)
  w_x = np.array([ row[0] for row in west_sort ]).astype(np.float)
  w_y = np.array([ row[1] for row in west_sort ]).astype(np.float)
  
  match_0_list = []
  match_1_list = []
  match_2_list = []
  match_3_list = []
  match_4_list = []
  single_e_w = []
  rad = 0.0147 #matching radius in degrees. 0.0147 deg = 52.92" = 2.3 kelt pixels
  for (i,cand) in enumerate(e_id):
    lowdec = e_dec[i] - rad
    highdec = e_dec[i] + rad
    cond = (w_dec >= lowdec) & (w_dec <= highdec)
    if (np.sum(cond)==0): #if no west objects match to the east object's dec
      match_0_list.append("%19s %d %19s %9.5f %9.5f %6.3f %6s %9.3f %9.3f %9s %9s" %(e_id[i],0, 'NULL', e_ra[i], e_dec[i],e_mag[i], 'NULL',e_x[i],e_y[i],'NULL', 'NULL') )
    if (np.sum(cond)>=1): #if one or more west objects matches to the east object's dec
      w_id_ = w_id[cond]
      w_x_ = w_x[cond]
      w_y_ = w_y[cond]
      w_ra_ = w_ra[cond]
      w_dec_ = w_dec[cond]
      w_file_line_ = west_sort[cond]
      w_mag_ = w_mag[cond]
      angsep = dAngSep(e_ra[i], e_dec[i],w_ra_, w_dec_)
      rcond = (angsep<rad)
      if (np.sum(rcond)==0): #if none of the western objects in the dec range are within the matching radius of the east object
        match_0_list.append("%19s %d %19s %9.5f %9.5f %6.3f %6s %9.3f %9.3f %9s %9s" %(e_id[i],0, 'NULL', e_ra[i], e_dec[i],e_mag[i], 'NULL',e_x[i],e_y[i],'NULL', 'NULL') )
      if (np.sum(rcond)==1): #if exactly one of the western objects in the dec range are within the matching radius of the east object
        asep = angsep[rcond]
        w_id_m = w_id_[rcond]
        w_x_m = w_x_[rcond]
        w_y_m = w_y_[rcond]
        w_ra_m = w_ra_[rcond]
        w_dec_m = w_dec_[rcond]
        w_file_line_m = w_file_line_[rcond]
        order = np.argsort(asep)
        asep_sort = asep[order]
        w_id_m_sort = w_id_m[order]
        w_x_m_sort = w_x_m[order]
        w_y_m_sort = w_y_m[order]
        w_ra_m_sort = w_ra_m[order]
        w_dec_m_sort = w_dec_m[order]
        w_mag_m = w_mag_[rcond]
        w_mag_m_sort = w_mag_m[order]
        if (abs(e_mag[i] - w_mag_m_sort[0]) <= 1.0): #counts as a match if there is a magnitude difference of less than one between the eastern object and the western object that matches by coordinates
          mean_ra = 0.5*(e_ra[i] + w_ra_m_sort[0]) #mean coords for east object and the west object that matches to it
          mean_dec = 0.5*(e_dec[i] + w_dec_m_sort[0])
          single_e_w.append(e_id[i] + " " + w_id_m_sort[0]) #used to keep track of all the 1:1 matches
          #match_1_list.append( "%s %d %9.5f %8.5f %2.3f %s" %(e_id[i],np.sum(rcond), e_ra[i], e_dec[i],e_mag[i], ' '.join(map(str,w_file_line_m[0]))) )
          match_1_list.append( "%19s %d %19s %9.5f %9.5f %6.3f %6.3f %9.3f %9.3f %9.3f %9.3f" %(e_id[i],0, w_id_m_sort[0], mean_ra, mean_dec,e_mag[i],w_mag_m_sort[0],e_x[i],e_y[i],w_x_m_sort[0],w_y_m_sort[0]) )
        if (abs(e_mag[i] - w_mag_m_sort[0]) > 1.0): #if the magnitude matching criteria is not met, there is no match
          match_0_list.append("%19s %d %19s %9.5f %9.5f %6.3f %6s %9.3f %9.3f %9s %9s" %(e_id[i],0, 'NULL', e_ra[i], e_dec[i],e_mag[i], 'NULL',e_x[i],e_y[i],'NULL','NULL') )
            
      if (np.sum(rcond)==2): #if 2 west objects match to one east object
        asep = angsep[rcond]
        w_id_m = w_id_[rcond]
        w_x_m = w_x_[rcond]
        w_y_m = w_y_[rcond]
        w_ra_m = w_ra_[rcond]
        w_dec_m = w_dec_[rcond]
        w_file_line_m = w_file_line_[rcond]
        w_mag_m = w_mag_[rcond]
        mag_diff = abs(e_mag[i] - w_mag_m)
        order = np.argsort(mag_diff)
        asep_sort = asep[order]
        w_id_m_sort = w_id_m[order]
        w_x_m_sort = w_x_m[order]
        w_y_m_sort = w_y_m[order]
        w_ra_m_sort = w_ra_m[order]
        w_dec_m_sort = w_dec_m[order]
        w_mag_m_sort = w_mag_m[order]
        mag_diff_sort = mag_diff[order]
        if (mag_diff_sort[1] <= 1.0): #if the greatest magnitude diff. between the 2 matching western objects is less than 1.0 mag, add the eastern object, plus both western objects to match_2_list
          #single_e_w.append(e_id[i] + " " + w_id_m_sort[0]) #not sure if this belongs here- let's see... hmm, probably not. check anyway. Maybe I need to have this sorted by angsep after all, and not by delta_mag
          match_2_list.append( "%19s %d %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f" \
          %(e_id[i],1, e_x[i],e_y[i],e_ra[i], e_dec[i], e_mag[i], w_id_m_sort[0],w_x_m_sort[0],w_y_m_sort[0], w_ra_m_sort[0], w_dec_m_sort[0], w_mag_m_sort[0], w_id_m_sort[1], w_x_m_sort[1],w_y_m_sort[1],w_ra_m_sort[1], w_dec_m_sort[1], w_mag_m_sort[1] ))
        if (mag_diff_sort[1] > 1.0) & (mag_diff_sort[0] <=0): #if the obj. w/ largest mag. diff. doesn't make the mag. cut, but the other does, append match_1_list w/ info for eastern object, and nearest (by mag.) western object
          match_1_list.append( "%19s %d %19s %9.5f %9.5f %6.3f %6.3f %9.3f %9.3f %9.3f %9.3f" %(e_id[i],1, w_id_m_sort[0], mean_ra, mean_dec,e_mag[i],w_mag_m_sort[0],e_x[i],e_y[i],w_x_m_sort[0],w_y_m_sort[0]) )
          #single_e_w.append(e_id[i] + " " + w_id_m_sort[0]) #used to keep track of 1:1 matches for comparison between e to w and w to e matching
        if (mag_diff_sort[0] > 1.0):
          match_0_list.append("%19s %d %19s %9.5f %9.5f %6.3f %6s %9.3f %9.3f %9s %9s" %(e_id[i],0, 'NULL', e_ra[i], e_dec[i],e_mag[i], 'NULL',e_x[i],e_y[i],'NULL', 'NULL') )

      if (np.sum(rcond)==3): #if 3 west objects match to one east object
        asep = angsep[rcond]
        w_id_m = w_id_[rcond]
        w_x_m = w_x_[rcond]
        w_y_m = w_y_[rcond]
        w_ra_m = w_ra_[rcond]
        w_dec_m = w_dec_[rcond]
        w_mag_m = w_mag_[rcond]
        w_file_line_m = w_file_line_[rcond]
        mag_diff = abs(e_mag[i] - w_mag_m)
        order = np.argsort(mag_diff)
        asep_sort = asep[order]
        w_id_m_sort = w_id_m[order]
        w_x_m_sort = w_x_m[order]
        w_y_m_sort = w_y_m[order]
        w_ra_m_sort = w_ra_m[order]
        w_dec_m_sort = w_dec_m[order]
        w_mag_m = w_mag_[rcond]
        w_mag_m_sort = w_mag_m[order]
        mag_diff_sort = mag_diff[order]
        if (mag_diff_sort[2] <= 1.0):
          #single_e_w.append(e_id[i] + " " + w_id_m_sort[0])
          match_3_list.append( "%19s %d %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f" \
          %(e_id[i],1, e_x[i],e_y[i],e_ra[i], e_dec[i], e_mag[i], w_id_m_sort[0],w_x_m_sort[0],w_y_m_sort[0], w_ra_m_sort[0], w_dec_m_sort[0], w_mag_m_sort[0], w_id_m_sort[1], w_x_m_sort[1],w_y_m_sort[1],w_ra_m_sort[1], w_dec_m_sort[1], w_mag_m_sort[1],w_id_m_sort[2], w_x_m_sort[2],w_y_m_sort[2],w_ra_m_sort[2], w_dec_m_sort[2], w_mag_m_sort[2] ))
        if (mag_diff_sort[2] > 1.0) & (mag_diff_sort[1] <=1.0):
          #single_e_w.append(e_id[i] + " " + w_id_m_sort[0])
          match_2_list.append( "%19s %d %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f" \
          %(e_id[i],1, e_x[i],e_y[i],e_ra[i], e_dec[i], e_mag[i], w_id_m_sort[0],w_x_m_sort[0],w_y_m_sort[0], w_ra_m_sort[0], w_dec_m_sort[0], w_mag_m_sort[0], w_id_m_sort[1], w_x_m_sort[1],w_y_m_sort[1],w_ra_m_sort[1], w_dec_m_sort[1], w_mag_m_sort[1] ))
        if (mag_diff_sort[0] <= 1.0):
          #single_e_w.append(e_id[i] + " " + w_id_m_sort[0])
          match_1_list.append( "%19s %d %19s %9.5f %9.5f %6.3f %6.3f %9.3f %9.3f %9.3f %9.3f" %(e_id[i],0, w_id_m_sort[0], mean_ra, mean_dec,e_mag[i],w_mag_m_sort[0],e_x[i],e_y[i],w_x_m_sort[0],w_y_m_sort[0]) )
        if (mag_diff_sort[0] > 1.0):
          match_0_list.append("%19s %d %19s %9.5f %9.5f %6.3f %6s %9.3f %9.3f %9s %9s" %(e_id[i],0, 'NULL', e_ra[i], e_dec[i],e_mag[i], 'NULL',e_x[i],e_y[i],'NULL', 'NULL') )


      if (np.sum(rcond)==4): #if 3 west objects match to one east object
        asep = angsep[rcond]
        w_id_m = w_id_[rcond]
        w_x_m = w_x_[rcond]
        w_y_m = w_y_[rcond]
        w_ra_m = w_ra_[rcond]
        w_dec_m = w_dec_[rcond]
        w_mag_m = w_mag_[rcond]
        w_file_line_m = w_file_line_[rcond]
        mag_diff = abs(e_mag[i] - w_mag_m)
        order = np.argsort(mag_diff)
        asep_sort = asep[order]
        w_id_m_sort = w_id_m[order]
        w_x_m_sort = w_x_m[order]
        w_y_m_sort = w_y_m[order]
        w_ra_m_sort = w_ra_m[order]
        w_dec_m_sort = w_dec_m[order]
        mag_diff_sort = mag_diff[order]
        w_mag_m_sort = w_mag_m[order]
        if (mag_diff_sort[3] <= 1.0):
          #single_e_w.append(e_id[i] + " " + w_id_m_sort[0])
          match_4_list.append(  "%19s %d %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f" \
          %(e_id[i],1, e_x[i],e_y[i],e_ra[i], e_dec[i], e_mag[i], w_id_m_sort[0],w_x_m_sort[0],w_y_m_sort[0], w_ra_m_sort[0], w_dec_m_sort[0], w_mag_m_sort[0], w_id_m_sort[1], w_x_m_sort[1],w_y_m_sort[1],w_ra_m_sort[1], w_dec_m_sort[1], w_mag_m_sort[1],w_id_m_sort[2], w_x_m_sort[2],w_y_m_sort[2],w_ra_m_sort[2], w_dec_m_sort[2], w_mag_m_sort[2],w_id_m_sort[3], w_x_m_sort[3],w_y_m_sort[3],w_ra_m_sort[3], w_dec_m_sort[3], w_mag_m_sort[3] ))
        if (mag_diff_sort[3] > 1.0) & (mag_diff_sort[2] <=1.0):
          #single_e_w.append(e_id[i] + " " + w_id_m_sort[0])
          match_3_list.append( "%19s %d %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f" \
          %(e_id[i],1, e_x[i],e_y[i],e_ra[i], e_dec[i], e_mag[i], w_id_m_sort[0],w_x_m_sort[0],w_y_m_sort[0], w_ra_m_sort[0], w_dec_m_sort[0], w_mag_m_sort[0], w_id_m_sort[1], w_x_m_sort[1],w_y_m_sort[1],w_ra_m_sort[1], w_dec_m_sort[1], w_mag_m_sort[1],w_id_m_sort[2], w_x_m_sort[2],w_y_m_sort[2],w_ra_m_sort[2], w_dec_m_sort[2], w_mag_m_sort[2] ))
        if (mag_diff_sort[2] > 1.0) & (mag_diff_sort[1] <=1.0):
          #single_e_w.append(e_id[i] + " " + w_id_m_sort[0])
          match_2_list.append( "%19s %d %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f" \
          %(e_id[i],1, e_x[i],e_y[i],e_ra[i], e_dec[i], e_mag[i], w_id_m_sort[0],w_x_m_sort[0],w_y_m_sort[0], w_ra_m_sort[0], w_dec_m_sort[0], w_mag_m_sort[0], w_id_m_sort[1], w_x_m_sort[1],w_y_m_sort[1],w_ra_m_sort[1], w_dec_m_sort[1], w_mag_m_sort[1] ))
        if (mag_diff_sort[1] > 1.0) & (mag_diff_sort[0] <= 1.0):
          match_1_list.append( "%19s %d %19s %9.5f %9.5f %6.3f %6.3f %9.3f %9.3f %9.3f %9.3f" %(e_id[i],1, w_id_m_sort[0], mean_ra, mean_dec,e_mag[i],w_mag_m_sort[0],e_x[i],e_y[i],w_x_m_sort[0],w_y_m_sort[0]) )
          #single_e_w.append(e_id[i] + " " + w_id_m_sort[0])
        if (mag_diff_sort[0] > 1.0):
          match_0_list.append("%19s %d %19s %9.5f %9.5f %6.3f %6s %9.3f %9.3f %9s %9s" %(e_id[i],0, 'NULL', e_ra[i], e_dec[i],e_mag[i], 'NULL',e_x[i],e_y[i],'NULL', 'NULL') )


#output columns in ./kelt_e_w_match/ are: e_id, # of coord matches, w_id, ra (avg.), dec (avg.), e_mag, w_mag, e_x, e_y, w_x, w_y
  np.savetxt('/media/sf_Astro/k_cat/single_e_w/{0}_ew_cat.dat'.format(field),(match_0_list + match_1_list),fmt='%s')
  np.savetxt('/media/sf_Astro/k_cat/multi_e_w/{0}_ew_cat.dat'.format(field),(match_2_list + match_3_list + match_4_list),fmt='%s')
  np.savetxt('/media/sf_Astro/k_cat/no_match_e_w/{0}_e_cat.dat'.format(field),match_0_list,fmt='%s') #do this so that later, with another script, I can make a single file with all east objects, and either 0 or 1 west matches, then do the same for west first... 
  print '%.1f' %(time.clock() - start_time) + " seconds"


## Changes from v7: changed all the paths to files to be correct for use in Unix box
