#!/usr/bin/python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*-coding:Utf-8 -*

"""
  syntheticpopulation.py
  Auteur : Alexis Petit, PhD Student, Namur University
  Date : 2017 03 22

  This scrip analyses a population of space debris in the GEO region.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os
from datetime import datetime
from scipy.stats import norm


def read_simulation(name_file):
  """ Read the snapshot of a simulation of the space debris environment.

  INPUT
  -----

  name_file: string
    Name of the file containing the population.

  RETURN
  ------

  population: dataframe
    The population of space debris.

  """

  f = open(name_file,'r')
  data = f.readlines()
  f.close()

  yy = int(data[1].strip().split()[0])
  mm = int(data[1].strip().split()[1])
  dd = int(data[1].strip().split()[2])
  hh = int(data[1].strip().split()[3])
  mi = int(data[1].strip().split()[4])
  ss = int(data[1].strip().split()[5])
  date = datetime(yy,mm,dd,hh,mi,ss)

  headers = data[2].strip().split()
  population = []
  for row in data[3:]:
    if 'GEOTEST' not in row:
      population.append(row.strip().split())
  population = pd.DataFrame(population,columns=headers)    

  category = {}
  for i in range(0,len(population)-1,1):
    name = population.loc[i]['NAME']
    if name not in category:
        category[name] = 1
    else:
        category[name] += 1
  print category

  plt.title(date.strftime('%y/%m/%d'),fontsize=14)
  plt.grid(True)
  plt.rc('xtick', labelsize=12) 
  plt.rc('ytick', labelsize=12) 
  plt.ticklabel_format(style='sci',useOffset=False)
  plt.xlabel('RAAN [deg]',fontsize=14)
  plt.ylabel('i [deg]',fontsize=14)

  for name_selected in category:
    flag = population['NAME']==name_selected
    raan = population['RAAN[deg]'][flag]
    inc = population['i[deg]'][flag]
    plt.plot(raan,inc,'o',markersize=2,label=name_selected.replace('_',' '))   
  plt.xlim([0,360])
  plt.ylim([0,20])
  plt.legend(loc='upper left',fontsize=14)
  plt.savefig('simulation')
  plt.close()
  
  return population
  
  
def read_controls(name_file,var_list):  

  f = open(name_file,'r')
  data = f.readlines()
  f.close()

  headers = data[2].strip().split()
  population = []
  for row in data[3:]:
    if 'GEOTEST' not in row:
      population.append(row.strip().split())
  population = pd.DataFrame(population,columns=headers)
     
  fig = plt.figure(figsize=(12,8))

  for k,var in enumerate(var_list):

    flag = population['NAME']=='EKRAN_2_DEB'
    var_data = [float(i) for i in population[var][flag].tolist()]
    
    ax = fig.add_subplot(2,3,k+1)
    ax.grid(True)
    ax.set_title(var)
    ax.hist(var_data,bins=25)
  plt.savefig('controls')


def discretize(variable,n_dim):
  """ Discretize a serie of data.

  INPUT
  -----

  variable: dataframe
  n_dim: integer

  RETURN
  ------

  variable_discretized:

  """

  variable = variable.convert_objects(convert_numeric=True)
  out = pd.cut(variable,n_dim)
  counts = pd.value_counts(out)
  variable_discretized = counts.reindex(out.cat.categories)
  
  return variable_discretized
  
 
def pop_2_cross_table(population,var_list,n_dim):
    """ Create the initial cross table.

    INPUT
    -----

    population: dataframe
    var_list: list
    n_dim: integer

    RETURN
    ------

    cross_table: matrix-array 
      The cross table generated from the initial population.

    """

    bounds = {}
    limits = {}
    total = {}
    for var in var_list:
      population[var] = population[var].convert_objects(convert_numeric=True)
      bounds[var],limits[var] = pd.cut(population[var],n_dim,retbins=True,labels=False)
      counts = pd.value_counts(bounds[var])
      counts = [counts.loc[i] for i in range(0,len(counts),1)]
      total[var] = counts
      
    cross_table = np.zeros((n_dim,n_dim)) 
    for i in range(0,len(limits['a[m]'])-1,1):
      for j in range(0,len(limits['RAAN[deg]'])-1,1):
        flag = (bounds['a[m]']==i) & (bounds['RAAN[deg]']==j)
        cross_table[i][j] = len(population[flag])
      
    return cross_table,total


def ipf_process(cross_table,total,controls,var_list):
    """ Apply the IPF process.

    """
  
    new_cross_table = np.zeros((n_dim,n_dim)) 

    return new_cross_table


def cross_table_2_pop(cross_table,var_list):

  population = []

  return population


# FIRST STEP: read the population coming from our simulation
  
population = read_simulation('simulation1.txt')
flag = population['NAME']=='EKRAN_2_DEB'
population = population[flag]

# SECOND STEP: create the cross-table

var_list = ['a[m]','RAAN[deg]']
n_dim = 5
cross_table,total = pop_2_cross_table(population,var_list,n_dim)
print 'Cross-table'
print cross_table

# THIRD STEP: calculate the controls
controls = read_controls('simulation1.txt',var_list)

# FOURTH STEP: apply the IPF process
new_cross_table = ipf_process(cross_table,total,controls,var_list)

# FITH STEP: compute the synthetic population
new_population = cross_table_2_pop(cross_table,var_list)
