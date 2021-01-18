import os
import numpy as np
import pandas as pd
import webbrowser
os.chdir('/Users/elise/Desktop/GitHub/Hubness_sc/Data_bulk/ARCHS4/')
os.getcwd()

# import list of suspicious GSE (with a sample dropout >= 60)
suspicious_list = np.array(pd.read_csv('GSE_list_dropout>=60.csv', usecols=[1]).T).tolist()[0]
all_list = np.array(pd.read_csv('GSE_list.csv', usecols=[1]).T).tolist()[0]
semi_all = np.array(pd.read_csv('GSM_list_dropout>=60.csv', usecols=[1]).T).tolist()[0]
len(suspicious_list)

# SCRAPING - merely open web browser page actually
urls = []
for i in semi_all:
    urls.append('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+i)
for url in urls:
    webbrowser.open(url)

