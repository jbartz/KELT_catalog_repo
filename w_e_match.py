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
#for file_no in range(2,4,2):
  east_file = kelt_cat_files[file_no]
  west_file = kelt_cat_files[file_no+1]
  field = east_file[4:7]
  #f = open('/media/sf_Astro/k_cat/kelt_cat_raw/{0}'.format(west_file),mode='r')
  #file_w_line = f.readlines() #is a list of lines from west_file
  #ar_file_w_line = np.array(file_w_line) #is an array of lines from west_file
  #print 'read all lines in file 2...'
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
  for (i,cand) in enumerate(w_id):
    lowdec = w_dec[i] - rad
    highdec = w_dec[i] + rad
    cond = (e_dec >= lowdec) & (e_dec <= highdec)
    if (np.sum(cond)==0): #if no west objects match to the east object's dec
#####THIS LINE BELOW IS GOOD, I BELIEVE. BE SURE TO FIX THE OTHER OUTPUTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      match_0_list.append("%19s %d %19s %9.5f %9.5f %6s %6.3f %9s %9s %9.3f %9.3f" %('NULL',0, w_id[i], w_ra[i], w_dec[i],'NULL',w_mag[i],'NULL','NULL',w_x[i],w_y[i]) )
#############################################################################################################################
    if (np.sum(cond)>=1): #if one or more west objects matches to the east object's dec
      e_id_ = e_id[cond]
      e_x_ = e_x[cond]
      e_y_ = e_y[cond]
      e_ra_ = e_ra[cond]
      e_dec_ = e_dec[cond]
      e_file_line_ = east_sort[cond]
      e_mag_ = e_mag[cond]
      angsep = dAngSep(w_ra[i], w_dec[i],e_ra_, e_dec_)
      rcond = (angsep<rad)
      if (np.sum(rcond)==0): #if none of the western objects in the dec range are within the matching radius of the east object
        match_0_list.append("%19s %d %19s %9.5f %9.5f %6s %6.3f %9s %9s %9.3f %9.3f" %('NULL',0, w_id[i], w_ra[i], w_dec[i],'NULL',w_mag[i],'NULL','NULL',w_x[i],w_y[i]) )
      if (np.sum(rcond)==1): #if exactly one of the western objects in the dec range are within the matching radius of the east object
        asep = angsep[rcond]
        e_id_m = e_id_[rcond]
        e_x_m = e_x_[rcond]
        e_y_m = e_y_[rcond]
        e_ra_m = e_ra_[rcond]
        e_dec_m = e_dec_[rcond]
        e_file_line_m = e_file_line_[rcond]
        order = np.argsort(asep)
        asep_sort = asep[order]
        e_id_m_sort = e_id_m[order]
        e_x_m_sort = e_x_m[order]
        e_y_m_sort = e_y_m[order]
        e_ra_m_sort = e_ra_m[order]
        e_dec_m_sort = e_dec_m[order]
        e_mag_m = e_mag_[rcond]
        e_mag_m_sort = e_mag_m[order]
        if (abs(w_mag[i] - e_mag_m_sort[0]) <= 1.0): #counts as a match if there is a magnitude difference of less than one between the eastern object and the western object that matches by coordinates
          mean_ra = 0.5*(w_ra[i] + e_ra_m_sort[0]) #mean coords for east object and the west object that matches to it
          mean_dec = 0.5*(w_dec[i] + e_dec_m_sort[0])
          single_e_w.append(e_id_m_sort[0]+ " " +w_id[i]) #used to keep track of all the 1:1 matches
          #match_1_list.append( "%s %d %9.5f %8.5f %2.3f %s" %(e_id[i],np.sum(rcond), e_ra[i], e_dec[i],e_mag[i], ' '.join(map(str,w_file_line_m[0]))) )
          match_1_list.append( "%19s %d %19s %9.5f %9.5f %6.3f %6.3f %9.3f %9.3f %9.3f %9.3f" %(e_id_m_sort[0],0, w_id[i], mean_ra, mean_dec,e_mag_m_sort[0],w_mag[i],e_x_m_sort[0],e_y_m_sort[0],w_x[i],w_y[i]) )
        if (abs(w_mag[i] - e_mag_m_sort[0]) > 1.0): #if the magnitude matching criteria is not met, there is no match
          match_0_list.append("%19s %d %19s %9.5f %9.5f %6s %6.3f %9s %9s %9.3f %9.3f" %('NULL',0, w_id[i], w_ra[i], w_dec[i],'NULL',w_mag[i],'NULL','NULL',w_x[i],w_y[i]) )
            
      if (np.sum(rcond)==2): #if 2 west objects match to one east object
        asep = angsep[rcond]
        e_id_m = e_id_[rcond]
        e_x_m = e_x_[rcond]
        e_y_m = e_y_[rcond]
        e_ra_m = e_ra_[rcond]
        e_dec_m = e_dec_[rcond]
        e_file_line_m = e_file_line_[rcond]
        order = np.argsort(asep)

        asep_sort = asep[order]
        e_id_m_sort = e_id_m[order]
        e_x_m_sort = e_x_m[order]
        e_y_m_sort = e_y_m[order]
        e_ra_m_sort = e_ra_m[order]
        e_dec_m_sort = e_dec_m[order]
        e_mag_m = e_mag_[rcond]
        e_mag_m_sort = e_mag_m[order]
        mag_diff = abs(w_mag[i] - e_mag_m)
        mag_diff_sort = mag_diff[order]
        if (mag_diff_sort[1] <= 1.0): #if the greatest magnitude diff. between the 2 matching western objects is less than 1.0 mag, add the eastern object, plus both western objects to match_2_list
          #single_e_w.append(e_id[i] + " " + w_id_m_sort[0]) #not sure if this belongs here- let's see... hmm, probably not. check anyway. Maybe I need to have this sorted by angsep after all, and not by delta_mag
          match_2_list.append( "%19s %d %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f" \
          %(w_id[i],1, w_x[i],w_y[i],w_ra[i], w_dec[i], w_mag[i], e_id_m_sort[0],e_x_m_sort[0],e_y_m_sort[0], e_ra_m_sort[0], e_dec_m_sort[0], e_mag_m_sort[0], e_id_m_sort[1], e_x_m_sort[1],e_y_m_sort[1],e_ra_m_sort[1], e_dec_m_sort[1], e_mag_m_sort[1] ))
        if (mag_diff_sort[1] > 1.0) & (mag_diff_sort[0] <=0): #if the obj. w/ largest mag. diff. doesn't make the mag. cut, but the other does, append match_1_list w/ info for eastern object, and nearest (by mag.) western object
          match_1_list.append( "%19s %d %19s %9.5f %9.5f %6.3f %6.3f %9.3f %9.3f %9.3f %9.3f" %(e_id_m_sort[0],0, w_id[i], mean_ra, mean_dec,e_mag_m_sort[0],w_mag[i],e_x_m_sort[0],e_y_m_sort[0],w_x[i],w_y[i]) )
          #single_e_w.append(e_id[i] + " " + w_id_m_sort[0]) #used to keep track of 1:1 matches for comparison between e to w and w to e matching
        if (mag_diff_sort[0] > 1.0):
          match_0_list.append("%19s %d %19s %9.5f %9.5f %6s %6.3f %9s %9s %9.3f %9.3f" %('NULL',0, w_id[i], w_ra[i], w_dec[i],'NULL',w_mag[i],'NULL','NULL',w_x[i],w_y[i]) )

      if (np.sum(rcond)==3): #if 3 west objects match to one east object
        asep = angsep[rcond]
        e_id_m = e_id_[rcond]
        e_x_m = e_x_[rcond]
        e_y_m = e_y_[rcond]
        e_ra_m = e_ra_[rcond]
        e_dec_m = e_dec_[rcond]
        e_file_line_m = e_file_line_[rcond]
        order = np.argsort(asep)
        asep_sort = asep[order]
        e_id_m_sort = e_id_m[order]
        e_x_m_sort = e_x_m[order]
        e_y_m_sort = e_y_m[order]
        e_ra_m_sort = e_ra_m[order]
        e_dec_m_sort = e_dec_m[order]
        e_mag_m = e_mag_[rcond]
        e_mag_m_sort = e_mag_m[order]
        mag_diff = abs(w_mag[i] - e_mag_m)
        mag_diff_sort = mag_diff[order]
        if (mag_diff_sort[2] <= 1.0):
          #single_e_w.append(e_id[i] + " " + w_id_m_sort[0])
          match_3_list.append( "%19s %d %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f" \
          %(w_id[i],1, w_x[i],e_y[i],w_ra[i], w_dec[i], w_mag[i], e_id_m_sort[0],e_x_m_sort[0],e_y_m_sort[0], e_ra_m_sort[0], e_dec_m_sort[0], e_mag_m_sort[0], e_id_m_sort[1], e_x_m_sort[1],e_y_m_sort[1],e_ra_m_sort[1], e_dec_m_sort[1], e_mag_m_sort[1],e_id_m_sort[2], e_x_m_sort[2],e_y_m_sort[2],e_ra_m_sort[2], e_dec_m_sort[2], e_mag_m_sort[2] ))
        if (mag_diff_sort[2] > 1.0) & (mag_diff_sort[1] <=1.0):
          #single_e_w.append(e_id[i] + " " + w_id_m_sort[0])
          match_2_list.append( "%19s %d %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f" \
          %(w_id[i],1, w_x[i],w_y[i],w_ra[i], w_dec[i], w_mag[i], e_id_m_sort[0],e_x_m_sort[0],e_y_m_sort[0], e_ra_m_sort[0], e_dec_m_sort[0], e_mag_m_sort[0], e_id_m_sort[1], e_x_m_sort[1],e_y_m_sort[1],e_ra_m_sort[1], e_dec_m_sort[1], e_mag_m_sort[1] ))
        if (mag_diff_sort[0] <= 1.0):
          #single_e_w.append(e_id[i] + " " + w_id_m_sort[0])
          match_1_list.append( "%19s %d %19s %9.5f %9.5f %6.3f %6.3f %9.3f %9.3f %9.3f %9.3f" %(e_id_m_sort[0],0, w_id[i], mean_ra, mean_dec,e_mag_m_sort[0],w_mag[i],e_x_m_sort[0],e_y_m_sort[0],w_x[i],w_y[i]) )
        if (mag_diff_sort[0] > 1.0):
          match_0_list.append("%19s %d %19s %9.5f %9.5f %6s %6.3f %9s %9s %9.3f %9.3f" %('NULL',0, w_id[i], w_ra[i], w_dec[i],'NULL',w_mag[i],'NULL','NULL',w_x[i],w_y[i]) )


      if (np.sum(rcond)==4): #if 3 west objects match to one east object
        asep = angsep[rcond]
        e_id_m = e_id_[rcond]
        e_x_m = e_x_[rcond]
        e_y_m = e_y_[rcond]
        e_ra_m = e_ra_[rcond]
        e_dec_m = e_dec_[rcond]
        e_file_line_m = e_file_line_[rcond]
        order = np.argsort(asep)
        asep_sort = asep[order]
        e_id_m_sort = e_id_m[order]
        e_x_m_sort = e_x_m[order]
        e_y_m_sort = e_y_m[order]
        e_ra_m_sort = e_ra_m[order]
        e_dec_m_sort = e_dec_m[order]
        e_mag_m = e_mag_[rcond]
        e_mag_m_sort = e_mag_m[order]
        mag_diff = abs(w_mag[i] - e_mag_m)
        mag_diff_sort = mag_diff[order]

        if (mag_diff_sort[3] <= 1.0):
          #single_e_w.append(e_id[i] + " " + w_id_m_sort[0])
          match_4_list.append(  "%19s %d %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f" \
          %(w_id[i],1, w_x[i],w_y[i],w_ra[i], w_dec[i], w_mag[i], e_id_m_sort[0],e_x_m_sort[0],e_y_m_sort[0], e_ra_m_sort[0], e_dec_m_sort[0], e_mag_m_sort[0], e_id_m_sort[1], e_x_m_sort[1],e_y_m_sort[1],e_ra_m_sort[1], e_dec_m_sort[1], e_mag_m_sort[1],e_id_m_sort[2], e_x_m_sort[2],e_y_m_sort[2],e_ra_m_sort[2], e_dec_m_sort[2], e_mag_m_sort[2],e_id_m_sort[3], e_x_m_sort[3],e_y_m_sort[3],e_ra_m_sort[3], e_dec_m_sort[3], e_mag_m_sort[3] ))
        if (mag_diff_sort[3] > 1.0) & (mag_diff_sort[2] <=1.0):
          #single_e_w.append(e_id[i] + " " + w_id_m_sort[0])
          match_3_list.append( "%19s %d %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f" \
          %(w_id[i],1, w_x[i],e_y[i],w_ra[i], w_dec[i], w_mag[i], e_id_m_sort[0],e_x_m_sort[0],e_y_m_sort[0], e_ra_m_sort[0], e_dec_m_sort[0], e_mag_m_sort[0], e_id_m_sort[1], e_x_m_sort[1],e_y_m_sort[1],e_ra_m_sort[1], e_dec_m_sort[1], e_mag_m_sort[1],e_id_m_sort[2], e_x_m_sort[2],e_y_m_sort[2],e_ra_m_sort[2], e_dec_m_sort[2], e_mag_m_sort[2] ))
        if (mag_diff_sort[2] > 1.0) & (mag_diff_sort[1] <=1.0):
          #single_e_w.append(e_id[i] + " " + w_id_m_sort[0])
          match_2_list.append( "%19s %d %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f %19s %9.3f %9.3f %9.5f %9.5f %6.3f" \
          %(w_id[i],1, w_x[i],w_y[i],w_ra[i], w_dec[i], w_mag[i], e_id_m_sort[0],e_x_m_sort[0],e_y_m_sort[0], e_ra_m_sort[0], e_dec_m_sort[0], e_mag_m_sort[0], e_id_m_sort[1], e_x_m_sort[1],e_y_m_sort[1],e_ra_m_sort[1], e_dec_m_sort[1], e_mag_m_sort[1] ))
        if (mag_diff_sort[1] > 1.0) & (mag_diff_sort[0] <= 1.0):
          match_1_list.append( "%19s %d %19s %9.5f %9.5f %6.3f %6.3f %9.3f %9.3f %9.3f %9.3f" %(e_id_m_sort[0],0, w_id[i], mean_ra, mean_dec,e_mag_m_sort[0],w_mag[i],e_x_m_sort[0],e_y_m_sort[0],w_x[i],w_y[i]) )
          #single_e_w.append(e_id[i] + " " + w_id_m_sort[0])
        if (mag_diff_sort[0] > 1.0):
          match_0_list.append("%19s %d %19s %9.5f %9.5f %6s %6.3f %9s %9s %9.3f %9.3f" %('NULL',0, w_id[i], w_ra[i], w_dec[i],'NULL',w_mag[i],'NULL','NULL',w_x[i],w_y[i]) )


#output columns in ./kelt_e_w_match/ are: e_id, # of coord matches, w_id, ra (avg.), dec (avg.), e_mag, w_mag, e_x, e_y, w_x, w_y
  np.savetxt('/media/sf_Astro/k_cat/single_w_e/{0}_we_cat.dat'.format(field),(match_0_list + match_1_list),fmt='%s')
  np.savetxt('/media/sf_Astro/k_cat/multi_w_e/{0}_we_cat.dat'.format(field),(match_2_list + match_3_list + match_4_list),fmt='%s')
  np.savetxt('/media/sf_Astro/k_cat/no_match_w_e/{0}_w_cat.dat'.format(field),match_0_list,fmt='%s') #do this so that later, with another script, I can make a single file with all east objects, and either 0 or 1 west matches, then do the same for west first... 
  print field + ' : ' + '%.1f' %((time.clock() - start_time)/60.) + " minutes"


## Changes from v7: changed all the paths to files to be correct for use in Unix box
