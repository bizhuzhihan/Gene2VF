#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import re
from urllib import request
import subprocess

# Genes
path = "/Users/zhuzhihan/Desktop/AMR/data/genes.xlsx"
profile = pd.read_excel(path)

# VFDB
import xlrd

VFs = xlrd.open_workbook("/Users/zhuzhihan/Downloads/VFs.xls", formatting_info=True)
vf = VFs.sheet_by_index(0)
vf_name = []
key = []
url = []
m = []
fun = []
bacteria = []
for row in range(2, 576):
    rowValues = vf.row_values(row, start_colx=0, end_colx=8)
    vf_name.append(rowValues[0])
    link = vf.hyperlink_map.get((row, 0))
    url.append('(No URL)' if link is None else link.url_or_path)
    m.append(rowValues[6])
    fun.append(rowValues[5])
    key.append(rowValues[7])
    bacteria.append(rowValues[2])
from urllib import request
import re

html = []

def reducing(d):
    d = re.sub(r'<.*?>', "", d)
    d = re.sub(r'\\t\\t\\t', "", d)
    d = re.sub(r'\\r\\n', "", d)
    return d

for i in range(len(url)):
    print(i)
    with request.urlopen(url[i]) as f:
        data = f.read()
        data = str(data)
        
        name = '<a name="'+vf_name[i]
        begin = re.search(name, data, flags=0)
        
        if (begin!=None):
            data = data[begin.span()[1]:]
            end = re.search(r'</P>', data, flags=0)
            data = data[:end.span()[0]]
            
            f = re.search(r'<b>Related genes:</b>', data, flags=0)
            if (f!=None):
                r = re.search(r'<b>Keywords:</b>', data, flags=0)
                if (r!=None):
                    html.append(reducing(data[f.span()[1]:r.span()[0]])) 
                else:
                    r = re.search(r'<b>Functions:</b>', data, flags=0)
                    html.append(reducing(data[f.span()[1]:r.span()[0]])) 
            else:
                html.append('None') 
        else:
            html.append('None') 

df = pd.DataFrame()
df['vf_name']=vf_name
df['bacteria']=bacteria
df['key']=key
df['gene']=html
df['mechanism']=m
df['function']=fun
df['link']=url
df.to_csv('/Users/zhuzhihan/Desktop/vfdb.csv')

# Merge
new_gene = []
new_count = []
new_VF = []
new_Bacteria = []
new_key = []
new_mechanism = []
new_function = []
new_link = []
for i in range(len(profile['gene'])):
    for j in range(len(df['gene'])):
        gg = re.search(profile['gene'][i], df['gene'][j], flags=0)
        if (gg!=None):
            new_gene.append(profile['gene'][i])
            new_count.append(profile['count'][i])
            new_VF.append(df['vf_name'][j])
            new_Bacteria.append(df['bacteria'][j])
            new_key.append(df['key'][j])
            new_mechanism.append(df['mechanism'][j])
            new_function.append(df['function'][j])
            new_link.append(df['link'][j])
df_new = pd.DataFrame()
df_new['gene']=new_gene
df_new['count']=new_count
df_new['vf']=new_VF
df_new['bacteria']=new_Bacteria
df_new['key']=new_key
df_new['mechanism']=new_mechanism
df_new['function']=new_function
df_new['link']=new_link
df_new.to_csv('/Users/zhuzhihan/Desktop/AMR/data/genes_vf_bacteria.csv')

# Filter the genes related to multiple VFs
genesi = df_new.groupby('gene').count()
genesi = genesi[genesi['vf']==1]
genesi['gene'] = genesi.index
genesi = genesi.drop(columns = ['count','vf','bacteria','key','mechanism','function','link'])
genesi.index = range(len(genesi['gene']))
genesi = genesi.sort_values('gene')
genesi = pd.merge(genesi, df_new, on='gene')
genem = df_new.groupby('gene').count()
genem = genem[genem['vf']>1]
genem['gene'] = genem.index
genem = genem.drop(columns = ['count','vf','bacteria','key','mechanism','function','link'])
genem.index = range(len(genem['gene']))
genem = pd.merge(genem, df_new, on='gene')

# Download the genomes from NCBI
b = pd.read_csv('/Users/zhuzhihan/Desktop/Bacteria.csv')
species = []
html = []
for i in range(len(b['Bacteria'])):
    print(i)
    url = "https://www.ncbi.nlm.nih.gov/genome/?term=" + b['Bacteria'][i]
    with request.urlopen(url) as f:
        data = f.read()
        data = str(data)
        marker = 'Download sequences in FASTA format for <a href="'
        begin = re.search(marker, data, flags=0)
        if (begin!=None):
            data = data[begin.span()[1]:]
            end = re.search(r'">genome</a>', data, flags=0)
            data = data[:end.span()[0]]
            species.append(b['Bacteria'][i])
            html.append(data) 
for i in range(0,len(html)):
    print(i)
    path = '/Users/zhuzhihan/Desktop/AMR/data/genome/' + species[i] + '.fna.gz'
    request.urlretrieve(html[i], path)

# Calculate the dist from KP
dist = []
sp = []
spr = []
for i in genem['bacteria']:
    j = i
    spr.append(j)
    i = re.sub(r'Klebsiella pneumoniae', r'Klebsiella pneumonia', i, count=0, flags=0)
    end = re.search(r" \(", i, flags=0)
    if(end!=None):
        sp.append(re.sub(r' ', r'+', i, count=0, flags=0)[:end.span()[0]])
    else:
        sp.append(re.sub(r' ', r'+', i, count=0, flags=0))
 
s = []
for i in set(sp):
    path = '/Users/zhuzhihan/Desktop/AMR/data/genome/' + i + '.fna'
    a = subprocess.run(['/Users/zhuzhihan/Desktop/AMR/data/mash', 'dist', path, '/Users/zhuzhihan/Desktop/AMR/data/genome/Klebsiella+pneumonia.fna'], capture_output=True)
    a = str(a.stdout)
    a = a.split(sep=r'\t')[2]
    dist.append(a)
    s.append(i)

Dist = pd.DataFrame()
Dist['bacteria'] = s
Dist['dist'] = dist
Dist.to_csv('/Users/zhuzhihan/Desktop/AMR/data/Dist.csv')
SP = pd.DataFrame()
SP['bacteria'] = sp
SP['bacteria2'] = spr
SP = (pd.merge(SP, Dist, on='bacteria'))
SP['bacteria'] = SP['bacteria2']
SP = SP.drop(columns = 'bacteria2')
SP = SP.drop_duplicates()
genem = (pd.merge( genem,SP, on='bacteria'))
genem = genem.sort_values('gene')
genem.index = range(len(genem['gene']))
genem.to_csv('/Users/zhuzhihan/Desktop/AMR/data/1vsM.csv')

df = genem
df_result=pd.DataFrame(columns=df.columns)
gene_list = df['gene'].unique()
for i in gene_list:
    df_sort = df[df['gene'] == i].sort_values(by = 'dist')
    df_result.loc[i] = df_sort.iloc[0]
df_result.index = range(len(df_result['gene']))
df_result.to_csv('/Users/zhuzhihan/Desktop/AMR/data/closest.csv')

profile = pd.concat([genesi,df_result.drop(columns = 'dist')])
profile = profile.sort_values('gene')
profile.index = range(len(profile['gene']))
profile.to_csv('/Users/zhuzhihan/Desktop/AMR/data/profile.csv')