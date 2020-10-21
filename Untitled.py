#!/usr/bin/env python
# coding: utf-8

# In[20]:


get_ipython().run_line_magic('matplotlib', 'nbagg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import csv
import pandas as pd


# In[4]:


fig, ax = plt.subplots()
x = np.linspace(norm.ppf(0.01),
                norm.ppf(0.99), 100)

ax.plot(x, norm.pdf(x),
       'r-', lw=5, alpha=0.6)

ax.plot(x, np.full(len(x), 0.2),
       'b-', lw=1)

fig.show()


# In[24]:


tsv = pd.read_csv("/home/balberti/Documents/S5/protéo/tp1/tp-proteomics/data/TCL_wt1.tsv",delimiter="\t")
tsv=pd.DataFrame(tsv)
tsv=pd.DataFrame.dropna


# In[ ]:





# In[ ]:




