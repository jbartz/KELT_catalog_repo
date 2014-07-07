#!/usr/bin/python
import os,sys
import numpy as np
import time
import matplotlib.pyplot as plt
from math import *
from os import listdir
from os.path import isfile, join
start_time = time.clock()

#better to do this for each canonical file seperately, I guess?
#so like, 
#for file_name in file_list:
#  info = file_name[whatever]
#canon_list = np.genfromtxt('./canon_all.dat',dtype=str)

#only want to read the 'rectangular' part of the overlap_all.dat file
col_range = range(21)
#overlap_list = np.genfromtxt('/media/sf_Astro/k_cat/overlap_all.dat',usecols=(col_range),dtype=str)
overlap_list = np.genfromtxt('/media/sf_Astro/k_cat/overlap_all.dat',usecols=(col_range),dtype=str)
overlap_list = overlap_list[:,1:]

fields=[]
files_canon = [ f for f in listdir('/media/sf_Astro/k_cat/single_canon/') if isfile(join('/media/sf_Astro/k_cat/single_canon/',f))]
for filename in files_canon:
  field = filename[:3]
  exec "file_%s = np.genfromtxt('/media/sf_Astro/k_cat/single_canon/{0}'.format(filename), dtype='|S44')" %field
  exec "remove_f_%s = []" %field
  exec "keep_f_%s = []" %field
  fields.append(field)
  f = open(('/media/sf_Astro/k_cat/single_canon/' + filename),mode='r')
  exec('file_%s_lines = f.readlines()' %field)
  #print 'read all lines in file 2...'
  f.close()
#  
#f = open('./new/canon_all.dat',mode='r')
#canon_list_actual_list = f.readlines()
##print 'read all lines in file 2...'
#f.close()

#this only looks at the two nearest (in coordinates) objects, and ignores any higher order matching between fields
#This can be improved by considering, for example, if one N01 object matches to two N02 objects, and
#the N01 object is farther from the edge of the detector, we would want to delete both N02 objects 
#from the canonical list??
e1id_over = overlap_list[:,0]
w1id_over = overlap_list[:,1]
e2id_over = overlap_list[:,10]
w2id_over = overlap_list[:,11]

#east and west id of the first field
#all objects in these canon lists belong to the same field
#e1id_canon = canon_list[:,0]
#w1id_canon = canon_list[:,2]
#cond = (dec2 >= lowdec) & (dec2 <= highdec)
#canon_no_overlap = [x for x in overlap_list if not x in a]

#when going through overlap list, iterate over each line, and read the contents of line to e1_id, w1_id, e2_id, w2_id, etc
#if there exist two id's for field 1: choose east id. Otherwise, choose the non-'null' id
#do the same for field 2 with choosing identifiers
#if there are both east and west id's in the line, choose the one closest to the edge. This is the file to search for and remove from the canonical list
removed_lines=[]
removed_lines2=[]
test_list=[]
start_time = time.clock()
counter=0
#for each line in the overlap_list, see which objects it matches to in the canon_list
for i, line in enumerate(overlap_list):
  e1id = line[0]
  w1id = line[1]
  e2id = line[10]
  w2id = line[11]
  
  
  #####FIX THIS!!!!!!!!!!!!!!!!!!!
  e1x = line[6]
  w1x = line[8]
  e2x = line[16]
  w2x = line[18]
  
  #in canon_file, we want to delete the duplicate closest to the edge, so make sure to get that in memory
  if e1id == 'NULL':
    c_id1 = w1id
    c_1x = w1x.astype(np.float)
    if w1id == 'NULL':
      print 'Error! both field 1 objects are "null" in line ' + i
  if e1id != 'NULL':
    c_id1 = e1id
    c_1x = e1x.astype(np.float)
  if e2id == 'NULL':
    c_id2 = w2id
    c_2x = w2x.astype(np.float)
    if w2id == 'NULL':
      print 'Error! both field 2 objects are "null" in line ' + i
  if e2id != 'NULL':
    c_id2 = e2id
    c_2x = e2x.astype(np.float)

#look into how many pixels!!! replace the value used here with accurate number
  if c_1x >= 2034:
    edge_dist1 = abs(4068 - c_1x)
  else:
    edge_dist1 = c_1x
  if c_2x >= 2034:
    edge_dist2 = abs(4068 - c_2x)
  else:
    edge_dist2 = c_2x
  if edge_dist1 >= edge_dist2:
    c_id = c_id2
    keep_id = c_id1
  else:
    c_id = c_id1
    keep_id = c_id2
  field = c_id[2:5]
  keep_field = keep_id[2:5]
  eval("remove_f_%s.append(c_id)" %field)
  eval("keep_f_%s.append(keep_id)" %keep_field)
    #I want this 'keep' list, so that I can go through the final, duplicate-removed, canonical file, and
    #somehow use sed to change the last entry in each line this hits from '0' -> '1', indicating
    #that the object has other matches and should be looked up
for j in fields:
  exec "this = '/media/sf_Astro/k_cat/sed_remove/sed_remove%s' " %j
  exec "that = '/media/sf_Astro/k_cat/sed_keep/sed_keep%s' " %j
  with open( this, 'a') as outputfile1:
    for line in eval("remove_f_%s" %j):
      outputfile1.write("/" + line + "/d" + "\n")
  with open( that, 'a') as outputfile2:
    for line in eval("keep_f_%s" %j):
      outputfile2.write("/" + line + "/s/0$/1/" + "\n")
outputfile1.close()
outputfile2.close()

#file 1 no dupes : removed all the lines in the 'to be deleted' file
#file 2 matched to all but 14 lines in 'to be deleted' file... why the difference?