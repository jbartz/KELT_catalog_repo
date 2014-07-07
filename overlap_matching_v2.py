#!/usr/bin/python
import os,sys
import numpy as np
import time
import matplotlib.pyplot as plt
from math import *
from os import listdir
from os.path import isfile, join

#treat id labels as strings, but specify their size as the largest id

#in northern fields, only expect adjacent fields to overlap

#create script/subroutine that takes two fields and determines the maximum and min RA and Dec, and then
#checks whether there is overlap between those coordinates between the two fields (both RA and dec overlap)
#will I need to hardcode the field size to do this?
#2 fields, arbitrarly placed on sky... can we determine overlap without knowing size of fields?
#just use 26 deg. as the field size- easier
#if one field has RA positions within > 334 degrees, watch out for RA -> 0

#output should be 1+2_overlap, 2+3_overlap, etc...

#for output coordinates for field 1, east west match, use the mean of (ra_1 and ra_2), and (dec_1, and dec_2) as RA, Dec
#precision to a tenth of an arcsecond


#once all overlap files are output, then compare those to eachother...

#for canonical output file:
#   key point is a single line for every object (things with a match)
#    for cases where 1 and 2 match: East_1, West_1, East_2, West_2, East_3, West_3, etc..., also: 
#  in the long term, build this with the option where this file will not necessarily be rectangular

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

kelt_cat_files = [ f for f in listdir('/media/sf_Astro/k_cat/single_canon/') if isfile(join('/media/sf_Astro/k_cat/single_canon/',f))]
for file_no in range(0,len(kelt_cat_files)-1):
#for file_no in range(1,2,1):
  
  
### Use this block when looking at overlap between 1 -> 2
  file1 = kelt_cat_files[file_no]
  file2 = kelt_cat_files[file_no+1]
### Use this block when looking at overlap between 1 -> 2  
  

### Use this block when looking at overlap between 2 -> 1
  #file2 = kelt_cat_files[file_no]
  #file1 = kelt_cat_files[file_no+1]
### Use this block when looking at overlap between 2 -> 1
  
  
  field1 = file1[:3]
  field2 = file2[:3]
  #f = open('./kelt_e_w_match/{0}'.format(file2),mode='r')
  #file2_line = f.readlines() #is a list of lines from west_file
  #ar_file2_line = np.array(file2_line) #is an array of lines from west_file
  ##print 'read all lines in file 2...'
  #f.close()  
  
#######Use this to get data for most adjacent files (like N01 and N02, N02 and N03,...)
####### The single_canon folder contains files for each field with unique entries, as opposed to 
####### files in single_all, which have duplicate lines
  data1 = np.genfromtxt('/media/sf_Astro/k_cat/single_canon/{0}'.format(file1), dtype='|S44')
  data2 = np.genfromtxt('/media/sf_Astro/k_cat/single_canon/{0}'.format(file2), dtype='|S44')

#######Use this to get data for non-adjacent files (like N01 and N13)
### be careful with the order. for 1 -> 2 matching, data1=N01, data2=N02. reverse for 2 -> 1
  #data1 = np.genfromtxt('/media/sf_Astro/k_cat/single_canon/N01_all.dat', dtype='|S44')
  #data2 = np.genfromtxt('/media/sf_Astro/k_cat/single_canon/N13_all.dat', dtype='|S44')

  #tacks on east or west, and field number: lc_xxxxx.dat -> e_02_lc_xxxxx.dat. Takes very little time to run this part.

#for some reason sorting like this seems to fuck up the east_1 RA...
#  data1_sort = data1[data1[:,3].argsort()] #sorts by RA
#  data2_sort = data2[data2[:,3].argsort()] #sorts by RA 
  data1_sort = data1
  data2_sort = data2

#When this runs, we only see 2:2 matches (1_e, 1_w):(2_e, 2_w) and 1:1 matches (1_e, null):(2_e,null)
  eid2 = np.array([ row[0] for row in data2_sort ])
  eid1 = np.array([ row[0] for row in data1_sort ])
  wid2 = np.array([ row[2] for row in data2_sort ])
  wid1 = np.array([ row[2] for row in data1_sort ])
  dec2 = np.array([ row[4] for row in data2_sort ]).astype(np.float)
  dec1 = np.array([ row[4] for row in data1_sort ]).astype(np.float)
  ra2 = np.array([ row[3] for row in data2_sort ]).astype(np.float)
  ra1 = np.array([ row[3] for row in data1_sort ]).astype(np.float)
  mag2e = np.array([ row[5] for row in data2_sort ])
  mag1e = np.array([ row[5] for row in data1_sort ])
  mag2w = np.array([ row[6] for row in data2_sort ])
  mag1w = np.array([ row[6] for row in data1_sort ])
  e_x1 = np.array([ row[7] for row in data1_sort ])
  e_y1 = np.array([ row[8] for row in data1_sort ])
  w_x1 = np.array([ row[9] for row in data1_sort ])
  w_y1 = np.array([ row[10] for row in data1_sort ])
  
  e_x2 = np.array([ row[7] for row in data2_sort ])
  e_y2 = np.array([ row[8] for row in data2_sort ])
  w_x2 = np.array([ row[9] for row in data2_sort ])
  w_y2 = np.array([ row[10] for row in data2_sort ])
  
  
  if ra2.max() > ra1.max():
    ra_min = ra2.min()
    ra_max = ra1.max()
  else: 
    ra_min = ra1.min()
    ra_max = ra2.max()  
    
  #if dec2.max() > dec1.max():
  #  dec_min = dec2.min()
  #  dec_max = dec1.max()
  #else: 
  #  dec_min = dec1.min()
  #  dec_max = dec2.max()     
  #gets the coordinate range that we want to search through


  #limits our search to objects only within the overlap region between two fields
  coord1_cond = (ra1 >= ra_min) & (ra1 <= ra_max) # & (dec1 >= dec_min) & (dec1 <= dec_max) 
  coord2_cond = (ra2 >= ra_min) & (ra2 <= ra_max) # & (dec2 >= dec_min) & (dec2 <= dec_max) 
  
  eid2 = eid2[coord2_cond]
  eid1 = eid1[coord1_cond]
  wid2 = wid2[coord2_cond]
  wid1 = wid1[coord1_cond]
  dec2 = dec2[coord2_cond]
  dec1 = dec1[coord1_cond]
  ra2 = ra2[coord2_cond]
  ra1 = ra1[coord1_cond]
  mag2e = mag2e[coord2_cond]
  mag1e = mag1e[coord1_cond] 
  mag2w = mag2w[coord2_cond]
  mag1w = mag1w[coord1_cond]
  e_x1 = e_x1[coord1_cond] 
  e_y1 = e_y1[coord1_cond]
  w_x1 = w_x1[coord1_cond] 
  w_y1 = w_y1[coord1_cond]
  
  e_x2 = e_x2[coord2_cond] 
  e_y2 = e_y2[coord2_cond]
  w_x2 = w_x2[coord2_cond] 
  w_y2 = w_y2[coord2_cond]
  
  match_0_list = []
  match_1_list = []
  match_2_list = []
  match_3_list = []
  match_4_list = []
  rad = 0.0147 #matching radius in degrees. 0.0147 deg = 52.9" = 2.3 kelt pixels

  for (i,cand) in enumerate(eid1):
    lowdec = dec1[i] - rad
    highdec = dec1[i] + rad
    cond = (dec2 >= lowdec) & (dec2 <= highdec)
    if (np.sum(cond)>=1):
      wid2_ = wid2[cond]
      eid2_ = eid2[cond]
      ra2_ = ra2[cond]
      dec2_ = dec2[cond]
      file2_line_ = data2_sort[cond]
      mag2e_ = mag2e[cond]
      mag2w_ = mag2w[cond]
      e_x2_ = e_x2[cond]
      e_y2_ = e_y2[cond]
      w_x2_ = w_x2[cond]
      w_y2_ = w_y2[cond]
      angsep = dAngSep(ra1[i], dec1[i],ra2_, dec2_)
      rcond = (angsep<rad)
#Do we note objects with 0 cross matches?      if (np.sum(rcond)==0):
#        match_0_list.append(' '.join(data1_sort[i]))
      if (np.sum(rcond)==1):
        wid2_m = wid2_[rcond]
        eid2_m = eid2_[rcond]
        asep = angsep[rcond]
        ra2_m = ra2_[rcond]
        dec2_m = dec2_[rcond]
        e_x2_m = e_x2_[rcond]
        e_y2_m = e_y2_[rcond]
        w_x2_m = w_x2_[rcond]
        w_y2_m = w_y2_[rcond]
        
#        file2_line_m = file2_line_[rcond] #is this used?
        order = np.argsort(asep)
        asep_sort = asep[order]
        wid2_m_sort = wid2_m[order]
        eid2_m_sort = eid2_m[order]
        ra2_m_sort = ra2_m[order]
        dec2_m_sort = dec2_m[order]
        mag2e_m = mag2e_[rcond]
        mag2e_m_sort = mag2e_m[order]
        mag2w_m = mag2w_[rcond]
        mag2w_m_sort = mag2w_m[order]
        
        dec2_m_sort = dec2_m[order]
        e_x2_m_sort = e_x2_m[order]
        e_y2_m_sort = e_y2_m[order]
        w_x2_m_sort = w_x2_m[order]
        w_y2_m_sort = w_y2_m[order]
#        if (abs(mag1e[i] - mag2e_m_sort[0]) <= 2.0) :
        mean_ra = 0.5 * (ra1[i] + ra2_m_sort[0])
        mean_dec = 0.5 * (dec1[i] + dec2_m_sort[0])
#        match_1_list.append( "%19s %19s %19s %19s %i %9.5f %9.5f %6s %6s %6s %6s %8s %8s %8s %8s %8s %8s %8s %8s " %(eid1[i],wid1[i],eid2_m_sort[0],wid2_m_sort[0],np.sum(rcond), mean_ra, mean_dec,mag1e[i],mag1w[i],mag2e_m_sort[0],mag2w_m_sort[0], e_x1[i],e_y1[i],w_x1[i], w_y1[i],e_x2_m_sort[0],e_y2_m_sort[0],w_x2_m_sort[0],w_y2_m_sort[0]) )
        match_1_list.append( "%i %19s %19s %9.5f %9.5f %6s %6s %8s %8s %8s %8s %19s %19s %9.5f %9.5f %6s %6s %8s %8s %8s %8s" %(np.sum(rcond),eid1[i],wid1[i], ra1[i], dec1[i], mag1e[i],mag1w[i], e_x1[i],e_y1[i],w_x1[i], w_y1[i],eid2_m_sort[0],wid2_m_sort[0],ra2_m_sort[0],dec2_m_sort[0],mag2e_m_sort[0],mag2w_m_sort[0],e_x2_m_sort[0],e_y2_m_sort[0],w_x2_m_sort[0],w_y2_m_sort[0]) )

        #np.savetxt('./kelt_overlap_match/N01_N13_ew_cat.dat',match_1_list,fmt='%s')
        #put in some code to change the last entry in this line from '0' -> '1', to show that this object has overlap matches. This needs to make it
        #to the final canonical file, so we can see if an object has other id's
        #data1
#  np.savetxt('./kelt_overlap_match/{0}_ew_cat.dat'.format(field1 + '_' + field2),match_1_list,fmt='%s')

      if (np.sum(rcond)==2):
        wid2_m = wid2_[rcond]
        eid2_m = eid2_[rcond]
        asep = angsep[rcond]
        ra2_m = ra2_[rcond]
        dec2_m = dec2_[rcond]
        e_x2_m = e_x2_[rcond]
        e_y2_m = e_y2_[rcond]
        w_x2_m = w_x2_[rcond]
        w_y2_m = w_y2_[rcond]
        
#        file2_line_m = file2_line_[rcond] #is this used?
        order = np.argsort(asep)
        asep_sort = asep[order]
        wid2_m_sort = wid2_m[order]
        eid2_m_sort = eid2_m[order]
        ra2_m_sort = ra2_m[order]
        dec2_m_sort = dec2_m[order]
        mag2e_m = mag2e_[rcond]
        mag2e_m_sort = mag2e_m[order]
        mag2w_m = mag2w_[rcond]
        mag2w_m_sort = mag2w_m[order]
        
        dec2_m_sort = dec2_m[order]
        e_x2_m_sort = e_x2_m[order]
        e_y2_m_sort = e_y2_m[order]
        w_x2_m_sort = w_x2_m[order]
        w_y2_m_sort = w_y2_m[order]
#        if (abs(mag1e[i] - mag2e_m_sort[0]) <= 2.0) :
        mean_ra = 0.5 * (ra1[i] + ra2_m_sort[0])
        mean_dec = 0.5 * (dec1[i] + dec2_m_sort[0])
        match_2_list.append( "%i %19s %19s %9.5f %9.5f %6s %6s %8s %8s %8s %8s %19s %19s %9.5f %9.5f %6s %6s %8s %8s %8s %8s %19s %19s %9.5f %9.5f %6s %6s %8s %8s %8s %8s" %(np.sum(rcond),eid1[i],wid1[i], ra1[i], dec1[i], mag1e[i],mag1w[i], e_x1[i],e_y1[i],w_x1[i], w_y1[i],eid2_m_sort[0],wid2_m_sort[0],ra2_m_sort[0],dec2_m_sort[0],mag2e_m_sort[0],mag2w_m_sort[0],e_x2_m_sort[0],e_y2_m_sort[0],w_x2_m_sort[0],w_y2_m_sort[1],eid2_m_sort[1],wid2_m_sort[1],ra2_m_sort[1],dec2_m_sort[1],mag2e_m_sort[1],mag2w_m_sort[1],e_x2_m_sort[1],e_y2_m_sort[1],w_x2_m_sort[1],w_y2_m_sort[1]) )

      if (np.sum(rcond)==3):
        wid2_m = wid2_[rcond]
        eid2_m = eid2_[rcond]
        asep = angsep[rcond]
        ra2_m = ra2_[rcond]
        dec2_m = dec2_[rcond]
        e_x2_m = e_x2_[rcond]
        e_y2_m = e_y2_[rcond]
        w_x2_m = w_x2_[rcond]
        w_y2_m = w_y2_[rcond]
        
#        file2_line_m = file2_line_[rcond] #is this used?
        order = np.argsort(asep)
        asep_sort = asep[order]
        wid2_m_sort = wid2_m[order]
        eid2_m_sort = eid2_m[order]
        ra2_m_sort = ra2_m[order]
        dec2_m_sort = dec2_m[order]
        mag2e_m = mag2e_[rcond]
        mag2e_m_sort = mag2e_m[order]
        mag2w_m = mag2w_[rcond]
        mag2w_m_sort = mag2w_m[order]
        
        dec2_m_sort = dec2_m[order]
        e_x2_m_sort = e_x2_m[order]
        e_y2_m_sort = e_y2_m[order]
        w_x2_m_sort = w_x2_m[order]
        w_y2_m_sort = w_y2_m[order]
#        if (abs(mag1e[i] - mag2e_m_sort[0]) <= 2.0) :
        mean_ra = 0.5 * (ra1[i] + ra2_m_sort[0])
        mean_dec = 0.5 * (dec1[i] + dec2_m_sort[0])
        match_3_list.append( "%i %19s %19s %9.5f %9.5f %6s %6s %8s %8s %8s %8s %19s %19s %9.5f %9.5f %6s %6s %8s %8s %8s %8s %19s %19s %9.5f %9.5f %6s %6s %8s %8s %8s %8s %19s %19s %9.5f %9.5f %6s %6s %8s %8s %8s %8s" %(np.sum(rcond),eid1[i],wid1[i], ra1[i], dec1[i], mag1e[i],mag1w[i], e_x1[i],e_y1[i],w_x1[i], w_y1[i],eid2_m_sort[0],wid2_m_sort[0],ra2_m_sort[0],dec2_m_sort[0],mag2e_m_sort[0],mag2w_m_sort[0],e_x2_m_sort[0],e_y2_m_sort[0],w_x2_m_sort[0],w_y2_m_sort[1],eid2_m_sort[1],wid2_m_sort[1],ra2_m_sort[1],dec2_m_sort[1],mag2e_m_sort[1],mag2w_m_sort[1],e_x2_m_sort[1],e_y2_m_sort[1],w_x2_m_sort[1],w_y2_m_sort[1], eid2_m_sort[2],wid2_m_sort[2],ra2_m_sort[2],dec2_m_sort[2],mag2e_m_sort[2],mag2w_m_sort[2],e_x2_m_sort[2],e_y2_m_sort[2],w_x2_m_sort[2],w_y2_m_sort[2]) )

      if (np.sum(rcond)==4):
        wid2_m = wid2_[rcond]
        eid2_m = eid2_[rcond]
        asep = angsep[rcond]
        ra2_m = ra2_[rcond]
        dec2_m = dec2_[rcond]
        e_x2_m = e_x2_[rcond]
        e_y2_m = e_y2_[rcond]
        w_x2_m = w_x2_[rcond]
        w_y2_m = w_y2_[rcond]
        
#        file2_line_m = file2_line_[rcond] #is this used?
        order = np.argsort(asep)
        asep_sort = asep[order]
        wid2_m_sort = wid2_m[order]
        eid2_m_sort = eid2_m[order]
        ra2_m_sort = ra2_m[order]
        dec2_m_sort = dec2_m[order]
        mag2e_m = mag2e_[rcond]
        mag2e_m_sort = mag2e_m[order]
        mag2w_m = mag2w_[rcond]
        mag2w_m_sort = mag2w_m[order]
        
        dec2_m_sort = dec2_m[order]
        e_x2_m_sort = e_x2_m[order]
        e_y2_m_sort = e_y2_m[order]
        w_x2_m_sort = w_x2_m[order]
        w_y2_m_sort = w_y2_m[order]
#        if (abs(mag1e[i] - mag2e_m_sort[0]) <= 2.0) :
        mean_ra = 0.5 * (ra1[i] + ra2_m_sort[0])
        mean_dec = 0.5 * (dec1[i] + dec2_m_sort[0])
        print 'Damn! ' + str(field1) + ' ' + str(field2)
        match_4_list.append( "%i %19s %19s %9.5f %9.5f %6s %6s %8s %8s %8s %8s %19s %19s %9.5f %9.5f %6s %6s %8s %8s %8s %8s %19s %19s %9.5f %9.5f %6s %6s %8s %8s %8s %8s %19s %19s %9.5f %9.5f %6s %6s %8s %8s %8s %8s %19s %19s %9.5f %9.5f %6s %6s %8s %8s %8s %8s" %(np.sum(rcond),eid1[i],wid1[i], ra1[i], dec1[i], mag1e[i],mag1w[i], e_x1[i],e_y1[i],w_x1[i], w_y1[i],eid2_m_sort[0],wid2_m_sort[0],ra2_m_sort[0],dec2_m_sort[0],mag2e_m_sort[0],mag2w_m_sort[0],e_x2_m_sort[0],e_y2_m_sort[0],w_x2_m_sort[0],w_y2_m_sort[1],eid2_m_sort[1],wid2_m_sort[1],ra2_m_sort[1],dec2_m_sort[1],mag2e_m_sort[1],mag2w_m_sort[1],e_x2_m_sort[1],e_y2_m_sort[1],w_x2_m_sort[1],w_y2_m_sort[1], eid2_m_sort[2],wid2_m_sort[2],ra2_m_sort[2],dec2_m_sort[2],mag2e_m_sort[2],mag2w_m_sort[2],e_x2_m_sort[2],e_y2_m_sort[2],w_x2_m_sort[2],w_y2_m_sort[3],eid2_m_sort[3],wid2_m_sort[3],ra2_m_sort[3],dec2_m_sort[3],mag2e_m_sort[3],mag2w_m_sort[3],e_x2_m_sort[3],e_y2_m_sort[3],w_x2_m_sort[3],w_y2_m_sort[3]) )


#use this line when doing 1 -> 2 overlap matching
  np.savetxt('/media/sf_Astro/k_cat/overlap_1_2/{0}_overlap.dat'.format(field1 + '_' + field2),(match_1_list + match_2_list + match_3_list + match_4_list),fmt='%s')
#  np.savetxt('/media/sf_Astro/k_cat/overlap_2_1/{0}_overlap.dat'.format(field2 + '_' + field1),(match_1_list + match_2_list + match_3_list + match_4_list),fmt='%s')


#  np.savetxt('./new/kelt_overlap_match_all/{0}_ew_cat.dat'.format('N01_N13'),(match_1_list + match_2_list + match_3_list + match_4_list),fmt='%s')


#use this line when doing N01 -> N13 overlap
#  np.savetxt('/media/sf_Astro/k_cat/overlap_1_2/{0}_overlap.dat'.format('N01_N13'),(match_1_list + match_2_list + match_3_list + match_4_list),fmt='%s')

#use this line when doing N13 -> N01 overlap
#  np.savetxt('/media/sf_Astro/k_cat/overlap_2_1/{0}_overlap.dat'.format('N13_N01'),(match_1_list + match_2_list + match_3_list + match_4_list),fmt='%s')

#  np.savetxt('./new/kelt_overlap_multi_reverse/{0}_overlap.dat'.format(field1 + '_' + field2),(match_2_list + match_3_list + match_4_list),fmt='%s')
  #np.savetxt('./kelt_e_w_multi_match/{0}_ew_cat.dat'.format(field),(match_2_list + match_3_list + match_4_list),fmt='%s')
  #np.savetxt('./kelt_e_w_no_matches/{0}_e_cat.dat'.format(field),match_0_list,fmt='%s')  
  print field1 + ' ' + field2 + ' : ' + '%.1f' %((time.clock() - start_time)/60.) + " minutes"
'''      
      if (np.sum(rcond)==2): #if 2 west objects match to one east object
        asep = angsep[rcond]
        ra2_m = ra2_[rcond]
        dec2_m = dec2_[rcond]
        file2_line_m = file2_line_[rcond]
        order = np.argsort(asep)
        asep_sort = asep[order]
        ra2_m_sort = ra2_m[order]
        dec2_m_sort = dec2_m[order]
        mag2_m = mag2_[rcond]
        mag2_m_sort = mag2_m[order]
        if (abs(mag1[i] - mag2_m_sort[0]) <= 1.0) & (abs(mag1[i] - mag2_m_sort[1]) <= 1.0):
          match_2_list.append( "%s %d %9.5f %8.5f %2.3f %s %s" %(id1[i],np.sum(rcond), ra1[i], dec1[i],\
          mag1[i], ' '.join(map(str,file2_line_m[0])),' '.join(map(str,file2_line_m[1]))) )
          match_1_list.append( "%s %d %9.5f %8.5f %2.3f %s" %(id1[i],np.sum(rcond), ra1[i], dec1[i],\
          mag1[i], ' '.join(map(str,file2_line_m[0]))) )
        if (abs(mag1[i] - mag2_m_sort[0]) <= 1.0) & (abs(mag1[i] - mag2_m_sort[1]) > 1.0):
          match_1_list.append( "%s %d %9.5f %8.5f %2.3f %s" %(id1[i],np.sum(rcond), ra1[i], dec1[i],\
          mag1[i], ' '.join(map(str,file2_line_m[0]))) )

      if (np.sum(rcond)==3): #if 3 west objects match to one east object
        asep = angsep[rcond]
        ra2_m = ra2_[rcond]
        dec2_m = dec2_[rcond]
        file2_line_m = file2_line_[rcond]
        order = np.argsort(asep)
        asep_sort = asep[order]
        ra2_m_sort = ra2_m[order]
        dec2_m_sort = dec2_m[order]
        mag2_m = mag2_[rcond]
        mag2_m_sort = mag2_m[order]
        match_3_list.append("%s %d %9.5f %8.5f %2.3f %s %s %s" %(id1[i],np.sum(rcond), ra1[i], dec1[i],\
        mag1[i], ' '.join(map(str,file2_line_m[0])),' '.join(map(str,file2_line_m[1])),' '.join(map(str,file2_line_m[2]))) )
      if (np.sum(rcond)==4): #if 3 west objects match to one east object
        asep = angsep[rcond]
        ra2_m = ra2_[rcond]
        dec2_m = dec2_[rcond]
        file2_line_m = file2_line_[rcond]
        order = np.argsort(asep)
        asep_sort = asep[order]
        ra2_m_sort = ra2_m[order]
        dec2_m_sort = dec2_m[order]
        mag2_m = mag2_[rcond]
        mag2_m_sort = mag2_m[order]
        match_4_list.append("%s %d %9.5f %8.5f %2.3f %s %s %s %s" %(id1[i],np.sum(rcond), ra1[i], dec1[i],\
        mag1[i],' '.join(map(str,file2_line_m[0])),' '.join(map(str,file2_line_m[1])),' '.join(map(str,file2_line_m[2])),\
        ' '.join(map(str,file2_line_m[3]))) )
'''


