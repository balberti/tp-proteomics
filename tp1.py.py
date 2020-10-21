#!/usr/bin/env python
# coding: utf-8

# In[2]:


get_ipython().run_line_magic('matplotlib', 'nbagg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import csv
import pandas as pd


# In[3]:


fig, ax = plt.subplots()
x = np.linspace(norm.ppf(0.01),
                norm.ppf(0.99), 100)

ax.plot(x, norm.pdf(x),
       'r-', lw=5, alpha=0.6)

ax.plot(x, np.full(len(x), 0.2),
       'b-', lw=1)

fig.show()


# In[4]:


tsv = pd.read_csv("/home/balberti/Documents/S5/protéo/tp1/tp-proteomics/data/TCL_wt1.tsv",delimiter="\t")
tsv=pd.DataFrame(tsv)
tsv=pd.DataFrame.dropna(tsv)
tsv.iloc[:, 3:7]=tsv.iloc[:, 3:7].astype(float)
tsv


# In[5]:


log2_CAR=tsv["Log2 Corrected Abundance Ratio"]
plt.figure()
plt.xlabel("Ratio d'abondance corrigée (base Log2)")
plt.ylabel("Nombre de ratio")
plt.title("Histogramme de la variable Log2 Corrected Abundance Ratio")
plt.hist(log2_CAR)


# In[6]:


mu=np.mean(log2_CAR)
mu


# In[7]:


S_2=np.var(log2_CAR)*np.sqrt(len(log2_CAR)/(len(log2_CAR)+1))
S_2


# In[8]:


fig, ax = plt.subplots()
hist = ax.hist(log2_CAR, bins=100) # draw histogram
x = np.linspace(min(log2_CAR), max(log2_CAR), 100) # generate PDF domain points
dx = hist[1][1] - hist[1][0] # Get single value bar height
scale = len(log2_CAR)*dx # scale accordingly
ax.plot(x, norm.pdf(x, mu, np.sqrt(S_2))*scale) # compute theoritical PDF and draw it

ax.set_xlabel("Log2 Corrected Abundance Ratio")
ax.set_ylabel("protein counts")
ax.set_title("Histogramme de log2 corrected abundance ratio et sa densité de probabilité")
plt.show()


# In[9]:


log10_pval=tsv["-LOG10 Adj.P-val"]
fig, ax = plt.subplots()
ax.scatter(log2_CAR,log10_pval)
ax.set_xlabel("log2 corrected Abundance Ratio")
ax.set_ylabel("-log10 adjusted P-value")
ax.set_title("Volcano plot")

log10_lim=3
log2_lim=mu+np.sqrt(S_2) #norm mu+sigma
ax.plot(ax.get_xlim(),(log10_lim,log10_lim),color="black",linestyle="dashed")
ax.plot((log2_lim,log2_lim),ax.get_ylim(),color="black",linestyle="dashed")
ax.plot()
plt.show()


# seuil : Alpha = 1*10-3  -log10 de 0.0001 (a gauche) / mu+sigma =-0.17 (a droite) ou -0.28 (tout ce qui est a mu + 1 ecart type est sur abondant)
# Donc en gros on prend que ce qui est dans le carré en haut à droite

# In[10]:


df_ORA=tsv.loc[(tsv['Log2 Corrected Abundance Ratio']>mu+np.sqrt(S_2)) & (tsv['-LOG10 Adj.P-val']>-np.log10(0.001))]
df_ORA["Accession"]


# In[11]:


from xml.etree.ElementTree import parse, dump
# Parse the E.Coli proteome XML Document
tree = parse('data/uniprot-proteome_UP000000625.xml')
root = tree.getroot()
ns = '{http://uniprot.org/uniprot}' # MANDATORY PREFIX FOR ANY SEARCH within document
# Store all entries aka proteins in a list of xml nodes
proteins = root.findall(ns + 'entry')
# Display the xml subtree of the first protein 
# dump(proteins[0])

# Find the xml subtree of a protein with accession "P31224"
identifiant=df_ORA["Accession"]
for iden in identifiant:
    for entry in proteins:
        accessions = entry.findall(ns+"accession")
        for acc in accessions:
            if acc.text == iden:
                dump(entry)
                break


# J'ai récupéré ces résultats et utilisé sublim text pour ne récupérer que les identifiants Go d'intérêts. C'est une opération manuelle mais je n'arrivais pas à le faire avec un code. Il était plus important pour moi de contourner le problème pour continuer le TP.

# In[12]:


go_id=open("data/Go_id.txt",'r')
liste=[]
for line in go_id:
    liste.append(line)
liste_seule=set(liste)
print(len(liste))
len(liste_seule)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




