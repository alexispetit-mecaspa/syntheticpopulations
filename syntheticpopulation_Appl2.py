#!/usr/bin/python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*-coding:Utf-8 -*

"""
  syntheticpopulation.py
  Auteur : Alexis Petit, PhD Student, Namur University
  Date : 2017 03 22
  Adaptation du code : Morgane Dumont, PhD Student, Namur University
  Date : 2018 02 26

  This scrip analyses a population of space debris in the GEO region.
"""

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import linear_model
import numpy as np
from scipy import stats
import sys
import math
import os
from datetime import datetime
from scipy.stats import norm,lognorm, ks_2samp
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
  
  
def read_controls(name_file,var_list,limits,nb_obj,change_n_obj):
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
      #mu,std = norm.fit(var_data)
      #laws[var] = [mu,std]
      #p = norm.pdf(x,mu,std)
      shape, loc, scale = lognorm.fit(var_data, floc=0)
      std = np.log(shape)
      mu = np.log(scale)
      p = lognorm.pdf(x,std,loc=0, scale = scale)
    else:
      mu,std = norm.fit(var_data)
      laws[var] = [mu,std]
      p = norm.pdf(x,mu,std)
    ax.plot(x,p,'k',linewidth=2)
    ax.set_xlim([xmin,xmax])

    if change_n_obj:
	# to fit the approximated laws:
    	frequencies[var] = create_controls(var,xmin,xmax,mu,std,nb_obj,limits[var])
    else:
	# to fit exactly the same structure as simu2
	frequencies[var]=create_cut(var,var_data,limits[var])

    #print 'freq ',var,frequencies[var]
    
  plt.savefig('controls')
  plt.close()
  flag = population['NAME']=='EKRAN_2_DEB'
  population = population[flag]
  return frequencies,laws, population


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
  ax.plot(raan,inc,'o',markersize=4,label='Simulation 1')

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
  ax.plot(raan,inc,'o',markersize=4,color='darkorange',label='Simulation 2')
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

  ax.set_xlim(limit_raan)
  ax.set_ylim(limit_inc)
  plt.legend(loc='upper left',fontsize=14)
  plt.savefig('comparison_simulation')
  plt.close()    

def create_cut(var,var_data,limits):
  """ aggregate var_data within classes defined by limits[var]
  """
  frequencies = np.zeros(len(limits)-1)
  for data in var_data:
	found = False
	index = 0
	if data<=limits[len(limits)-1]:
		while found==False:
			if data <= limits[index+1]:
				frequencies[index]+=1
				found=True
			else: 
				index+=1
	else:
		frequencies[len(limits)-2]+=1
  return frequencies.tolist()	
	
def def_bounds(population,var_list,n_dim, quantile_cut):
    fig = plt.figure(figsize=(10,8))
    bounds = {}
    limits = {}
    frequencies = {}
    k=0
    for var in var_list:
      #print var
      
      
      population[var] = pd.to_numeric(population[var])
      if quantile_cut == True:
      	bounds[var],limits[var] = pd.qcut(population[var],n_dim[var],retbins=True,labels=False)
      else: 
      	bounds[var],limits[var] = pd.cut(population[var],n_dim[var],retbins=True,labels=False)
      #print population[var]
      #print limits[var]
      #print bounds[var]
      #print frequencies[var]
      countsTemp = pd.value_counts(bounds[var])
      
      print 'count avant', countsTemp
      counts = [0 for i in range(0,n_dim[var])]
      for i in  range(0,len(counts),1):
	if i in countsTemp.index:
		counts[i]= countsTemp.loc[i]
      print 'count apres', counts

      #print counts
      frequencies[var] = counts
      #print frequencies[var]	

      ax = fig.add_subplot(2,2,k+1)
      #ax.spines["top"].set_visible(False)
      #ax.spines["right"].set_visible(False)
      data_max = len(population[var])
          
      if var == 'A/M[m2/kg]':
        n, bins, patches = plt.hist(population[var],50,facecolor='green',alpha=0.75,normed=True)#,weights=weights)
        shape, loc, scale = lognorm.fit(population[var], floc=0)
        xmin =  min(population[var])
        xmax =  max(population[var])
        x = np.linspace(xmin,xmax,50)
        p = lognorm.pdf(x,shape,loc=0, scale = scale)
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
      
    plt.savefig('initial_frequencies_testMo')
    plt.close()
    return frequencies,limits,bounds


def create_controls(var,vmin,vmax,mu,std,nb_obj,limits):
  """ From statistical data we compute frequencies

  """
 
  counter = 0
  pop = []
  while(len(pop)<nb_obj):
    if var=='A/M[m2/kg]':
      x = np.random.lognormal(mean=mu,sigma=std)
      #print x
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

 
def pop_2_cross_table(population,var_list,n_dim,quantile_cut):
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


    frequencies, limits, bounds = def_bounds(population,var_list,n_dim, quantile_cut)
      
    #sys.exit()
    
    cross_table = np.zeros((n_dim['a[m]'],n_dim['i[deg]'],n_dim['RAAN[deg]'],n_dim['A/M[m2/kg]'])) 
    for i in range(0,len(limits['a[m]'])-1,1):
      for j in range(0,len(limits['i[deg]'])-1,1):
        for k in range(0,len(limits['RAAN[deg]'])-1,1):
          for l in range(0,len(limits['A/M[m2/kg]'])-1,1):
            flag = (bounds['a[m]']==i) & (bounds['i[deg]']==j) & (bounds['RAAN[deg]']==k) & (bounds['A/M[m2/kg]']==l)
            cross_table[i,j,k,l] = len(population[flag])

    return cross_table,frequencies,limits

def delete_zero_cell_problem(cross_table,frequencies, n_dim):
	for i in range(0,len(frequencies['a[m]']),1):
		#print i
		#print frequencies
		for j in range(0,len(frequencies['i[deg]']),1):
			for k in range(0,len(frequencies['RAAN[deg]']),1):
				for l in range(0,len(frequencies['A/M[m2/kg]']),1):
					if cross_table[i,j,k,l] == 0:
						cross_table[i,j,k,l] = 0.001
						frequencies['a[m]'][i] += cross_table[i,j,k,l]
						frequencies['i[deg]'][j] += cross_table[i,j,k,l]
						frequencies['RAAN[deg]'][k] += cross_table[i,j,k,l]
						frequencies['A/M[m2/kg]'][l] += cross_table[i,j,k,l]
	#print cross_table
	return cross_table, frequencies


def ipf_process(cross_table,frequencies,controls, n_dim):
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

    #cross_table, frequencies = delete_zero_cell_problem(cross_table,frequencies, n_dim)
      
    epsilon = 1E-13 # 10^-13 because eps machine is 10^-16 for each cell of contingency table (+-1000)
    print n_dim['a[m]']*n_dim['i[deg]']*n_dim['RAAN[deg]']*n_dim['A/M[m2/kg]']
    distance = 1
    list_distance = []
    counter = 0

    while (epsilon<distance) and (counter<200):
      #print "iteration",counter+1
      cross_table_saved = cross_table.copy()
      
      cross_table, frequencies = update_variable_cross_table(cross_table,frequencies,controls) 
      total = frequencies
      
      distance = 0
      for i in range(0,n_dim['a[m]'],1):
        for j in range(0,n_dim['i[deg]'],1):
          for k in range(0,n_dim['RAAN[deg]'],1):
            for l in range(0,n_dim['A/M[m2/kg]'],1):
              distance += np.abs(cross_table[i,j,k,l]-cross_table_saved[i,j,k,l])
              
      #print distance
      list_distance.append(distance)          
      counter += 1

    #print cross_table
    print 'fin IPF', frequencies

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
          if frequency['i[deg]'][j]==0:
            cross_table[i,j,k,l] = 0
          else:
            cross_table[i,j,k,l] = cross_table[i,j,k,l]*controls['i[deg]'][j]/frequency['i[deg]'][j]
          new_frequency['a[m]'][i] += cross_table[i,j,k,l]
          new_frequency['i[deg]'][j] += cross_table[i,j,k,l]      
          new_frequency['RAAN[deg]'][k] += cross_table[i,j,k,l]
          new_frequency['A/M[m2/kg]'][l] += cross_table[i,j,k,l]
  frequency = new_frequency.copy() 
 

  new_frequency = frequency.copy()
  for item in new_frequency:
    dim = len(new_frequency[item])
    new_frequency[item] = np.zeros(dim).tolist()


  for i in range(0,len(frequency['a[m]']),1):
    for j in range(0,len(frequency['i[deg]']),1):
      for k in range(0,len(frequency['RAAN[deg]']),1):
        for l in range(0,len(frequency['A/M[m2/kg]']),1):
          if frequency['RAAN[deg]'][k]!=0:
			cross_table[i,j,k,l] = cross_table[i,j,k,l]*controls['RAAN[deg]'][k]/frequency['RAAN[deg]'][k]
          new_frequency['a[m]'][i] += cross_table[i,j,k,l]
          new_frequency['i[deg]'][j] += cross_table[i,j,k,l]      
          new_frequency['RAAN[deg]'][k] += cross_table[i,j,k,l]
          new_frequency['A/M[m2/kg]'][l] += cross_table[i,j,k,l]
  frequency = new_frequency.copy()

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
  print 'laws :',laws
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
            ecc = 0.1*random()
            while(True):
              raan = np.random.normal(laws['RAAN[deg]'][0],laws['RAAN[deg]'][1],1)
              if (raan>=limits['RAAN[deg]'][k] and raan<limits['RAAN[deg]'][k+1]):
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
  raan_init = pd.to_numeric(population['RAAN[deg]'])
  inc_init = pd.to_numeric(population['i[deg]'])
  print 'initial correlation between RAAN and i for sim 1 ', np.corrcoef(raan_init,inc_init)
  #print len(population)
  
  # SECOND STEP: create the cross-table
  print
  print 'Compute the cross-table'
  var_list = ['a[m]','i[deg]','RAAN[deg]','A/M[m2/kg]']
  n_dim = {}
  n_dim['a[m]'] = 1
  n_dim['i[deg]'] = 10
  n_dim['RAAN[deg]'] = 10
  n_dim['A/M[m2/kg]'] = 1#3
  quantile_cut = False
  change_n_obj = False
  nb_obj = 920#*2

  cross_table,frequencies,limits = pop_2_cross_table(population,var_list,n_dim,quantile_cut)
  #print 'Cross-table'
  #print cross_table
  print 'Frequencies'
  print frequencies
  print 'Limits'
  print limits
  print

  # THIRD STEP: calculate the controls
  print
  print 'Compute the contraints'
  controls,laws, population_constraint = read_controls('simulation2.txt',var_list,limits,nb_obj,change_n_obj)
  raan_init_cons = pd.to_numeric(population_constraint['RAAN[deg]'])
  inc_init_cons = pd.to_numeric(population_constraint['i[deg]'])
  print 'initial correlation between RAAN and i for sim 2 ', np.corrcoef(raan_init_cons,inc_init_cons)
  print '> ',controls
  print 
  nb_obj = sum(controls['a[m]'])
  print nb_obj

  # FOURTH STEP: apply the IPF process
  new_cross_table = ipf_process(cross_table,frequencies,controls, n_dim)
  print 'New cross-table'
  #print cross_table
  
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

  print type(raan_init)
  print 'new correlation between RAAN and i', np.corrcoef([a[0] for a in raan],[a[0] for a in inc])
  raan_approx = np.random.normal(laws['RAAN[deg]'][0],laws['RAAN[deg]'][1],9600)
  inc_approx = np.random.normal(laws['i[deg]'][0],laws['i[deg]'][1],9600)
  print 'distrib raan - approx :', ks_2samp(raan_approx, [a[0] for a in raan])
  print 'distrib inc - approx :', ks_2samp(inc_approx, [a[0] for a in inc])

  ##############################################################
  # linear reg
  slope_init, intercept_init, r_value, p_value, std_err = stats.linregress(raan_init,inc_init)
  print('Coefficients linear reg sim 1: \n', slope_init, intercept_init)
  print('r_value : \n',r_value)
  slope_init_cons, intercept_init_cons, r_value, p_value, std_err = stats.linregress(raan_init_cons,inc_init_cons)
  print('Coefficients linear reg sim 1: \n', slope_init_cons, intercept_init_cons)
  print('r_value : \n',r_value)
  slope_synthPop, intercept_synthPop, r_value, p_value, std_err = stats.linregress([a[0] for a in raan],[a[0] for a in inc])
  print('Coefficients linear reg synthetic pop: \n', slope_synthPop, intercept_synthPop)
  print('r_value : \n',r_value)


  fig, ax = plt.subplots()
  ax.set_xlabel('RAAN [deg]',fontsize=14)
  ax.set_ylabel('i [deg]',fontsize=14)
  ax.spines["top"].set_visible(False)
  ax.spines["right"].set_visible(False)  
  ax.plot(population['RAAN[deg]'],population['i[deg]'],'o',markersize=4,label='Simulation ('+str(len(population))+' fragments)',alpha=0.7)
  ax.plot(raan,inc,'o',markersize=4,color='g',label='Synthetic population ('+str(len(raan))+' fragments)',alpha=0.7)
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

################ figure with everything and regression
  fig, ax = plt.subplots()
  ax.set_xlabel('RAAN [deg]',fontsize=14)
  ax.set_ylabel('i [deg]',fontsize=14)
  ax.spines["top"].set_visible(False)
  ax.spines["right"].set_visible(False)  
  ax.plot(population['RAAN[deg]'],population['i[deg]'],'o',markersize=4,color='blue',label='Simulation 1 ('+str(len(population))+' fragments)',alpha=0.7)
  X = np.array([min(raan_init), max(raan_init)])
  plt.plot(X, X*slope_init + intercept_init, 'b')
  ax.plot(population_constraint['RAAN[deg]'],population_constraint['i[deg]'],'o',markersize=4,color='darkorange',label='Simulation 2 ('+str(len(population_constraint))+' fragments)',alpha=0.7)
  plt.plot(X, X*slope_init_cons + intercept_init_cons ,color='darkorange')
  ax.plot(raan,inc,'o',markersize=4,label='Synthetic population ('+str(len(raan))+' fragments)',alpha=0.7,color='g')
  plt.plot(X, X*slope_synthPop + intercept_synthPop ,color='g')
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
  plt.savefig('synthetic_population_all')
  plt.close()


  plt.subplot(1, 2, 1)
  plt.hist([a[0] for a in raan], histtype='stepfilled',color="orange", alpha=0.4, normed=True, bins=40, label="Synthetic population")
  plt.hist(raan_init, histtype='stepfilled', alpha=0.4, normed=True, bins=40, label="Simulation 1")
  plt.hist(raan_init_cons, histtype='stepfilled', alpha=0.4, normed=True, bins=40, label="Simulation 2")
  plt.xlabel('RAAN [deg]')
  plt.ylabel('Density')
  plt.legend(loc='upper left',fontsize=12)
  plt.ylim((0,1))
  plt.subplot(1, 2, 2)
  plt.hist([a[0] for a in inc], histtype='stepfilled',color="orange", alpha=0.4, normed=True, bins=40, label="Synthetic population")
  plt.hist(inc_init, histtype='stepfilled', alpha=0.4, normed=True, bins=40, label="Simulation 1")
  plt.hist(inc_init_cons, histtype='stepfilled', alpha=0.4, normed=True, bins=40, label="Simulation 2")
  plt.xlabel('i[deg]')
  plt.ylim((0,1.4))
  plt.savefig('Distributions_pop_appl2')
  plt.close()


  limit_raan = [limits['RAAN[deg]'][0],limits['RAAN[deg]'][k+1]]
  limit_inc = [limits['i[deg]'][0],limits['i[deg]'][j+1]]
  comparison_clouds('simulation1.txt','simulation2.txt',limit_raan,limit_inc)
