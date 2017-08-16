#!/usr/bin/python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*-coding:Utf-8 -*

"""
  syntheticpopulation.py
  Auteur : Alexis Petit, PhD Student, Namur University
  Date : 2017 03 22

  This scrip analyses a population of space debris in the GEO region.
"""

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os
from datetime import datetime
from scipy.stats import norm,lognorm
from random import random


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
  population['A/M[m2/kg]'] = population['BC[m2/kg]'].astype('float64')/2.2
  
  category = {}
  for i in range(0,len(population)-1,1):
    name = population.loc[i]['NAME']
    if name not in category:
        category[name] = 1
    else:
        category[name] += 1

  print
  print '> Categories in the simulation:'      
  print '> ',category
  print

  fig, ax = plt.subplots()
  ax.set_title(date.strftime('%y/%m/%d'),fontsize=14)
  ax.set_xlabel('RAAN [deg]',fontsize=14)
  ax.set_ylabel('i [deg]',fontsize=14)
  ax.spines["top"].set_visible(False)
  ax.spines["right"].set_visible(False)

  for name_selected in category:
    flag = population['NAME']==name_selected
    raan = population['RAAN[deg]'][flag]
    inc = population['i[deg]'][flag]
    ax.plot(raan,inc,'o',markersize=2,label=name_selected.replace('_',' '))   
  ax.set_xlim([0,360])
  ax.set_ylim([0,20])
  plt.legend(loc='upper left',fontsize=14)
  plt.savefig('set_simulation_1')
  plt.close()
  
  return population
  
  
def read_controls(name_file,var_list,limits,nb_obj):
  """ Infer contraints from data

  INPUTS
  ------

  name_file: string 
    Name of the files used as contraints
  var_list: list
    Variables used [a,i,RAAN,...]
  limits: list

  RETURN
  ------

  controls: dataframe
    Frequencies used as contraints

  """

  #controls empty
  frequencies = {}
  laws = {}

  #read the files containing the contraints
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
  population['A/M[m2/kg]'] = population['BC[m2/kg]'].astype('float64')/2.2
  
  category = {}
  for i in range(0,len(population)-1,1):
    name = population.loc[i]['NAME']
    if name not in category:
        category[name] = 1
    else:
        category[name] += 1

  fig, ax = plt.subplots()
  ax.set_title(date.strftime('%y/%m/%d'),fontsize=14)
  ax.set_xlabel('RAAN [deg]',fontsize=14)
  ax.set_ylabel('i [deg]',fontsize=14)
  ax.spines["top"].set_visible(False)
  ax.spines["right"].set_visible(False)

  for name_selected in category:
    flag = population['NAME']==name_selected
    raan = population['RAAN[deg]'][flag]
    inc = population['i[deg]'][flag]
    ax.plot(raan,inc,'o',markersize=2,label=name_selected.replace('_',' '))   
  ax.set_xlim([0,360])
  ax.set_ylim([0,20])
  plt.legend(loc='upper left',fontsize=14)
  plt.savefig('set_simulation_2')
  plt.close()  

  print '> Plot the distribution (constraints)' 
  fig = plt.figure(figsize=(10,8))

  for k,var in enumerate(var_list):

    #we select only the objects with the tag 'EKRAN_2_DEB'
    flag = population['NAME']=='EKRAN_2_DEB'
    var_data = [float(i) for i in population[var][flag].tolist()]

    ax = fig.add_subplot(2,2,k+1)
    ax.grid(True)
    ax.set_title(var)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    #ax.hist(var_data,bins=25)
    xmin =  limits[var][0]
    xmax =  limits[var][-1]
    x = np.linspace(xmin, xmax, 100)
    if var == 'A/M[m2/kg]':
      #param = lognorm.fit(var_data)
      #laws[var] = param
      #p = lognorm.pdf(x,param[0])
      mu,std = norm.fit(var_data)
      laws[var] = [mu,std]
      p = norm.pdf(x,mu,std)
    else:
      mu,std = norm.fit(var_data)
      laws[var] = [mu,std]
      p = norm.pdf(x,mu,std)
    ax.plot(x,p,'k',linewidth=2)
    ax.set_xlim([xmin,xmax])

    frequencies[var] = create_controls(var,xmin,xmax,mu,std,nb_obj,limits[var])
    #print 'freq',frequencies[var]
    
  plt.savefig('controls')
  plt.close()
  
  return frequencies,laws


def comparison_clouds(name1,name2,limit_raan,limit_inc):

  #read the files containing the contraints
  f = open(name1,'r')
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
  population['A/M[m2/kg]'] = population['BC[m2/kg]']
  
  category = {}
  for i in range(0,len(population)-1,1):
    name = population.loc[i]['NAME']
    if name not in category:
        category[name] = 1
    else:
        category[name] += 1

  fig, ax = plt.subplots()
  ax.set_title(date.strftime('%y/%m/%d'),fontsize=14)
  ax.set_xlabel('RAAN [deg]',fontsize=14)
  ax.set_ylabel('i [deg]',fontsize=14)
  ax.spines["top"].set_visible(False)
  ax.spines["right"].set_visible(False)

  flag = population['NAME']=='EKRAN_2_DEB'
  raan = population['RAAN[deg]'][flag]
  inc = population['i[deg]'][flag]
  ax.plot(raan,inc,'o',markersize=2,label='Simulation 1')

  #read the files containing the contraints
  f = open(name2,'r')
  data = f.readlines()
  f.close()

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

  flag = population['NAME']=='EKRAN_2_DEB'
  raan = population['RAAN[deg]'][flag]
  inc = population['i[deg]'][flag]
  ax.plot(raan,inc,'o',markersize=2,label='Simulation 2')
    
  ax.set_xlim(limit_raan)
  ax.set_ylim(limit_inc)
  plt.legend(loc='upper left',fontsize=14)
  plt.savefig('comparison_simulation')
  plt.close()    


def create_controls(var,vmin,vmax,mu,std,nb_obj,limits):
  """ From statistical data we compute frequencies

  """
 
  counter = 0
  pop = []
  while(len(pop)<nb_obj):
    if var=='A/M[m2/kg]':
      x = np.random.normal(loc=mu,scale=std)
    else:
      x = np.random.normal(loc=mu,scale=std)
    if x>vmin and x<vmax:
      pop.append(x)
      
  frequencies = np.zeros(len(limits)-1)
  for i in range(0,len(limits)-1,1):
    for var in pop:
      if (var>limits[i]) and (var<limits[i+1]):
        frequencies[i] += 1      
    
  return frequencies.tolist()

 
def pop_2_cross_table(population,var_list,n_dim):
    """ Create the initial cross table.

    INPUT
    -----

    population: dataframe
    var_list: list
    n_dim: integer

    RETURN
    ------

    cross_table: float-matrix 
      The cross table generated from the initial population.
    frequencies: dataframe
      Count the frequencies for each variables.

    """

    fig = plt.figure(figsize=(10,8))
    
    bounds = {}
    limits = {}
    frequencies = {}
    k=0
    for var in var_list:
      population[var] = population[var].convert_objects(convert_numeric=True)
      bounds[var],limits[var] = pd.cut(population[var],n_dim[var],retbins=True,labels=False)
      counts = pd.value_counts(bounds[var])
      counts = [counts.loc[i] for i in range(0,len(counts),1)]
      frequencies[var] = counts

      ax = fig.add_subplot(2,2,k+1)
      #ax.spines["top"].set_visible(False)
      #ax.spines["right"].set_visible(False)
      data_max = len(population[var])
          
      if var == 'A/M[m2/kg]':
        mu,sigma = norm.fit(population[var])
        #weights = np.ones_like(population[var])/float(len(population[var]))
        n, bins, patches = plt.hist(population[var],50,facecolor='green',alpha=0.75,normed=1)#,weights=weights)
        p = mlab.normpdf(bins,mu,sigma)#*data_max
        #plt.plot(bins,p,'r--',linewidth=2)
        param = lognorm.fit(population[var])
        xmin =  min(population[var])
        xmax =  max(population[var])
        x = np.linspace(xmin,xmax,50)
        p = lognorm.pdf(x,param[0])
        plt.plot(x,p,'r--',linewidth=2)
      else:
        mu,sigma = norm.fit(population[var])
        weights = np.ones_like(population[var])/float(len(population[var]))
        n, bins, patches = plt.hist(population[var],50,facecolor='green',alpha=0.75,normed=1)#,weights=weights)
        p = mlab.normpdf(bins,mu,sigma)#*data_max
        plt.plot(bins,p,'r--',linewidth=2)
        
      for i in range(0,n_dim[var]+1,1):
        x = min(population[var])+((max(population[var])-min(population[var])))*i/n_dim[var]
        plt.axvline(x,color='k',linestyle='--')

      k += 1
      #ax.grid(True)
      ax.set_title(var)
      
    plt.savefig('initial_frequencies')
    plt.close()
    #sys.exit()
    
    cross_table = np.zeros((n_dim['a[m]'],n_dim['i[deg]'],n_dim['RAAN[deg]'],n_dim['A/M[m2/kg]'])) 
    for i in range(0,len(limits['a[m]'])-1,1):
      for j in range(0,len(limits['i[deg]'])-1,1):
        for k in range(0,len(limits['RAAN[deg]'])-1,1):
          for l in range(0,len(limits['A/M[m2/kg]'])-1,1):
            flag = (bounds['a[m]']==i) & (bounds['i[deg]']==j) & (bounds['RAAN[deg]']==k) & (bounds['A/M[m2/kg]']==l)
            cross_table[i,j,k,l] = len(population[flag])
      
    return cross_table,frequencies,limits


def ipf_process(cross_table,frequencies,controls):
    """ Apply the IPF process.

    INPUTS
    ------

    cross_table: float-matrix
    frequencies: dataframe
    controls: dataframe
    var_list: list

    RETURN
    ------

    new_cross_table: float-matrix

    """
      
    epsilon = 10E-8
    distance = 1
    list_distance = []
    counter = 0

    while (epsilon<distance) and (counter<100):

      cross_table_saved = cross_table.copy()
      
      cross_table,total = update_variable_cross_table(cross_table,frequencies,controls) 
      
      distance = 0
      for i in range(0,n_dim['a[m]'],1):
        for j in range(0,n_dim['i[deg]'],1):
          for k in range(0,n_dim['RAAN[deg]'],1):
            for l in range(0,n_dim['A/M[m2/kg]'],1):
              distance += np.abs(cross_table[i,j,k,l]-cross_table_saved[i,j,k,l])

      list_distance.append(distance)          
      counter += 1

    #print cross_table

    fig, ax = plt.subplots()
    ax.grid(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    #plt.rc('xtick', labelsize=12) 
    #plt.rc('ytick', labelsize=12) 
    #plt.ticklabel_format(style='sci',useOffset=False)
    ax.set_xlabel('Iteration',fontsize=14)
    ax.set_ylabel('Distance',fontsize=14)
    ax.set_xlim([0,len(list_distance)-1])
    ax.set_yscale('log')
    ax.plot(list_distance)   
    plt.savefig('convergence')
    plt.close()   
      
    return cross_table


def update_variable_cross_table(cross_table,frequency,controls):
  """ Fitting of the cross table

  INPUTS
  ------

  cross_table: float-matrix
  frequency: dataframe
  controls: dataframes

  RETURN
  ------

  cross_table: float-matrix
  new_frequency: dataframe
      
  """
 
  new_frequency = frequency.copy()
  for item in new_frequency:
    dim = len(new_frequency[item])
    new_frequency[item] = np.zeros(dim).tolist()
  for i in range(0,len(frequency['a[m]']),1):
    for j in range(0,len(frequency['i[deg]']),1):
      for k in range(0,len(frequency['RAAN[deg]']),1):
        for l in range(0,len(frequency['A/M[m2/kg]']),1):
          if frequency['a[m]'][i]==0:
            cross_table[i,j,k,l]=0
          else:
            cross_table[i,j,k,l] = cross_table[i,j,k,l]*controls['a[m]'][i]/frequency['a[m]'][i]
          new_frequency['a[m]'][i] += cross_table[i,j,k,l]
          new_frequency['i[deg]'][j] += cross_table[i,j,k,l]
          new_frequency['RAAN[deg]'][k] += cross_table[i,j,k,l]
          new_frequency['A/M[m2/kg]'][l] += cross_table[i,j,k,l]
  frequency = new_frequency.copy()
  #print cross_table
  
  new_frequency = frequency.copy()
  for item in new_frequency:
    dim = len(new_frequency[item])
    new_frequency[item] = np.zeros(dim).tolist()  
  for i in range(0,len(frequency['a[m]']),1):
    for j in range(0,len(frequency['i[deg]']),1):
      for k in range(0,len(frequency['RAAN[deg]']),1):      
        for l in range(0,len(frequency['A/M[m2/kg]']),1):
          if frequency['i[deg]'][j]==0:
            cross_table[i,j,k,l] = 0
          else:
            cross_table[i,j,k,l] = cross_table[i,j,k,l]*controls['i[deg]'][j]/frequency['i[deg]'][j]
          new_frequency['a[m]'][i] += cross_table[i,j,k,l]
          new_frequency['i[deg]'][j] += cross_table[i,j,k,l]      
          new_frequency['RAAN[deg]'][k] += cross_table[i,j,k,l]
          new_frequency['A/M[m2/kg]'][l] += cross_table[i,j,k,l]
  frequency = new_frequency.copy() 
  #print cross_table
  
  new_frequency = frequency.copy()
  for item in new_frequency:
    dim = len(new_frequency[item])
    new_frequency[item] = np.zeros(dim).tolist()
  for i in range(0,len(frequency['a[m]']),1):
    for j in range(0,len(frequency['i[deg]']),1):
      for k in range(0,len(frequency['RAAN[deg]']),1):
        for l in range(0,len(frequency['A/M[m2/kg]']),1):
          if frequency['RAAN[deg]'][k]==0:
            cross_table[i,j,k,l] = 0
          else:
            cross_table[i,j,k,l] = cross_table[i,j,k,l]*controls['RAAN[deg]'][k]/frequency['RAAN[deg]'][k]
          new_frequency['a[m]'][i] += cross_table[i,j,k,l]
          new_frequency['i[deg]'][j] += cross_table[i,j,k,l]      
          new_frequency['RAAN[deg]'][k] += cross_table[i,j,k,l]
          new_frequency['A/M[m2/kg]'][l] += cross_table[i,j,k,l]
  frequency = new_frequency.copy()
  #print cross_table

  new_frequency = frequency.copy()
  for item in new_frequency:
    dim = len(new_frequency[item])
    new_frequency[item] = np.zeros(dim).tolist()
  for i in range(0,len(frequency['a[m]']),1):
    for j in range(0,len(frequency['i[deg]']),1):
      for k in range(0,len(frequency['RAAN[deg]']),1):
        for l in range(0,len(frequency['A/M[m2/kg]']),1):
          if frequency['A/M[m2/kg]'][l]==0:
            cross_table[i,j,k,l] = 0
          else:
            cross_table[i,j,k,l] = cross_table[i,j,k,l]*controls['A/M[m2/kg]'][l]/frequency['A/M[m2/kg]'][l]
          new_frequency['a[m]'][i] += cross_table[i,j,k,l]
          new_frequency['i[deg]'][j] += cross_table[i,j,k,l]      
          new_frequency['RAAN[deg]'][k] += cross_table[i,j,k,l]
          new_frequency['A/M[m2/kg]'][l] += cross_table[i,j,k,l]
  frequency = new_frequency.copy()

  
  return cross_table,new_frequency

def trs(cross_table,n_dim,nb_obj):

  new_cross_table = np.zeros((n_dim['a[m]'],n_dim['i[deg]'],n_dim['RAAN[deg]'],n_dim['A/M[m2/kg]'])) 
  weight = np.zeros((n_dim['a[m]'],n_dim['i[deg]'],n_dim['RAAN[deg]'],n_dim['A/M[m2/kg]']))
  sum_weight = 0

  total = 0
  for i in range(0,n_dim['a[m]'],1):
    for j in range(0,n_dim['i[deg]'],1):
      for k in range(0,n_dim['RAAN[deg]'],1):
        for l in range(0,n_dim['A/M[m2/kg]'],1):
          new_cross_table[i,j,k,l] += int(cross_table[i,j,k,l])
          total += new_cross_table[i,j,k,l] 
          weight[i,j,k,l] = cross_table[i,j,k,l] - new_cross_table[i,j,k,l]
          sum_weight += weight[i,j,k,l]

  bound = np.zeros((n_dim['a[m]'],n_dim['i[deg]'],n_dim['RAAN[deg]'],n_dim['A/M[m2/kg]']))         
  sum_bound = 0      
  for i in range(0,n_dim['a[m]'],1):
    for j in range(0,n_dim['i[deg]'],1):
      for k in range(0,n_dim['RAAN[deg]'],1):
        for l in range(0,n_dim['A/M[m2/kg]'],1):
          weight[i,j,k,l] = weight[i,j,k,l]/sum_weight
          sum_bound += weight[i,j,k,l]
          bound[i,j,k,l] = sum_bound

  while (total<nb_obj):
    flag = False
    for i in range(0,n_dim['a[m]'],1):
      for j in range(0,n_dim['i[deg]'],1):
        for k in range(0,n_dim['RAAN[deg]'],1):
          for l in range(0,n_dim['A/M[m2/kg]'],1):
            if (bound[i,j,k,l] > random()) and new_cross_table[i,j,k,l]!=0:
               new_cross_table[i,j,k,l] +=1
               total += 1
               flag = True
            if flag:
              break
          if flag:
            break
        if flag:
          break
      if flag:
        break
      
  return new_cross_table
  

def cross_table_2_pop(cross_table,frequencies,limits,laws,var_list):
  """ Convert the cross table to a population of objects

  """

  population = []
  for i in range(0,len(frequencies['a[m]']),1):
    for j in range(0,len(frequencies['i[deg]']),1):
      for k in range(0,len(frequencies['RAAN[deg]']),1):
        for l in range(0,len(frequencies['A/M[m2/kg]']),1):
          counter = 0
          while counter<cross_table[i,j,k,l]:
            sma = limits['a[m]'][i] + (limits['a[m]'][i+1]-limits['a[m]'][i])*random()
            while(True):
              inc = np.random.normal(laws['i[deg]'][0],laws['i[deg]'][1],1)
              if (inc>=limits['i[deg]'][j] and inc<limits['i[deg]'][j+1]):
                break
              #x = random()
              #inc = limits['i[deg]'][j] + (limits['i[deg]'][j+1]-limits['i[deg]'][j])*random()
              #p = norm(laws['i[deg]'][0], laws['i[deg]'][1]).pdf(inc)
              #if x<p:
              #  break             
            ecc = 0.1*random()
            while(True):
              x = random()
              raan = limits['RAAN[deg]'][k] + (limits['RAAN[deg]'][k+1]-limits['RAAN[deg]'][k])*random()
              p = norm(laws['RAAN[deg]'][0], laws['RAAN[deg]'][1]).pdf(raan)
              if x<p:
                break
            omega = 360.*random()
            ma = 360.*random()
            bc = limits['A/M[m2/kg]'][l] + (limits['A/M[m2/kg]'][l+1]-limits['A/M[m2/kg]'][l])*random()
            population.append([sma,ecc,inc,raan,omega,ma,bc])
            counter += 1
  
  #['a[m]', 'e[-]', 'i[deg]', 'RAAN[deg]', 'Omega[deg]', 'MA[deg]', 'A/M[m2/kg]', 'S[m]', 'IDNORAD', 'NAME']      
  headers = ['a[m]', 'e[-]', 'i[deg]', 'RAAN[deg]', 'Omega[deg]', 'MA[deg]', 'A/M[m2/kg]']
  population = pd.DataFrame(population,columns=headers) 
        
  return population


if __name__ == "__main__": 

  # FIRST STEP: read the population coming from our simulation
  print
  print 'Read the data of the simulation'
  population = read_simulation('simulation1.txt')
  flag = population['NAME']=='EKRAN_2_DEB'
  population = population[flag]
  
  # SECOND STEP: create the cross-table
  print
  print 'Compute the cross-table'
  var_list = ['a[m]','i[deg]','RAAN[deg]','A/M[m2/kg]']
  n_dim = {}
  n_dim['a[m]'] = 4
  n_dim['i[deg]'] = 7
  n_dim['RAAN[deg]'] = 7
  n_dim['A/M[m2/kg]'] = 3

  cross_table,frequencies,limits = pop_2_cross_table(population,var_list,n_dim)
  print 'Cross-table'
  print cross_table
  print 'Frequencies'
  print frequencies
  print 'Limits'
  print limits
  print

  # THIRD STEP: calculate the controls
  nb_obj = 460#*2
  print
  print 'Compute the contraints'
  controls,laws = read_controls('simulation2.txt',var_list,limits,nb_obj)
  print '> ',controls
  print 

  # FOURTH STEP: apply the IPF process
  new_cross_table = ipf_process(cross_table,frequencies,controls)
  print 'New cross-table'
  print cross_table
  
  # TRS process
  cross_table = trs(new_cross_table,n_dim,nb_obj)
  #cross_table = new_cross_table
  
  #check total
  total = 0
  for i in range(0,n_dim['a[m]'],1):
    for j in range(0,n_dim['i[deg]'],1):
      for k in range(0,n_dim['RAAN[deg]'],1):
        for l in range(0,n_dim['A/M[m2/kg]'],1):
          total += cross_table[i,j,k,l]
  print 'New total = ',total

  # FITH STEP: compute the synthetic population
  new_population = cross_table_2_pop(cross_table,frequencies,limits,laws,var_list)
  raan = new_population['RAAN[deg]']
  inc = new_population['i[deg]']

  fig, ax = plt.subplots()
  ax.set_xlabel('RAAN [deg]',fontsize=14)
  ax.set_ylabel('i [deg]',fontsize=14)
  ax.spines["top"].set_visible(False)
  ax.spines["right"].set_visible(False)  
  ax.plot(population['RAAN[deg]'],population['i[deg]'],'x',markersize=2,label='Simulation ('+str(len(population))+' fragments)',alpha=0.7)
  ax.plot(raan,inc,'o',markersize=2,label='Synthetic population ('+str(len(raan))+' fragments)',alpha=0.7)
  for j in range(0,n_dim['i[deg]'],1):
    x = limits['i[deg]'][j]
    ax.axhline(x, color='k', linestyle='--')
  x = limits['i[deg]'][j+1]
  ax.axhline(x, color='k', linestyle='--')
  for k in range(0,n_dim['RAAN[deg]'],1):
    x = limits['RAAN[deg]'][k]
    ax.axvline(x, color='k', linestyle='--')
  x = limits['RAAN[deg]'][k+1]
  ax.axvline(x, color='k', linestyle='--')
  ax.set_xlim([limits['RAAN[deg]'][0],limits['RAAN[deg]'][k+1]])
  ax.set_ylim([limits['i[deg]'][0],limits['i[deg]'][j+1]])
  plt.legend(loc='upper left',fontsize=12)
  plt.savefig('synthetic_population')
  plt.close()

  limit_raan = [limits['RAAN[deg]'][0],limits['RAAN[deg]'][k+1]]
  limit_inc = [limits['i[deg]'][0],limits['i[deg]'][j+1]]
  comparison_clouds('simulation1.txt','simulation2.txt',limit_raan,limit_inc)
