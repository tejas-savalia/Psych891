#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy
import matplotlib.pyplot as plt
import csv


# noise() generates standard normal noise.
# 
# resp_time(m, x, c) is the linearly decreasing response time function with noise. Linear as observed in the Serial Reacion Time Task
# 
# 
# dist_covered() returns the distance from the next stimulus which corresponds to the boundary separation parameter. This is equivalent to log(x) with noise added. 
# 
# 
# Data is generated and stored in the data.csv file -- this code need not be run.

# In[2]:


def noise():
    return numpy.random.normal()


# In[3]:


def resp_time(m, x, c):
    y = m*x + c
    return -y + noise()


# In[4]:


trials = 104
m = -1
c = 4
rt = list()


# In[5]:


for i in range(trials, c, -1):
    rt.append(resp_time(m, i, c))


# In[14]:


def dist_covered(x):
    return 5*numpy.log(x) + noise()


# In[15]:


dist = list()


# In[16]:


for i in range(c, trials):
    dist.append(dist_covered(i))


# In[17]:


plt.plot(dist)


# In[10]:


plt.plot(rt)


# In[18]:


with open('data.csv', mode = 'w') as data:
    data_writer = csv.writer(data)
    data_writer.writerow(rt)
    data_writer.writerow(dist)


# In[ ]:




