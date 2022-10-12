import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm
from functools import reduce
import os
import itertools

organs = ['brain', 'colon', 'heart', 'kidney', 'liver', 'lung','soleus','spleen','testis']
min_cell_number = 300
threshold_number = 4
universe = pd.read_csv('data/universe_annotated.csv')
deg = pd.read_csv('data/deg_annotated.csv')
output = {}


def do_corr(df, filename, threshold_number):
    f = f'./results/corr/{filename}.csv'
    if os.path.isfile(f):
        print('correlation is been reading from the file.')
        corr = pd.read_csv(f, index_col=0)
    else:
        corr = df.corr(method='pearson')
        corr.index.names = ['genes']
    stack = corr.stack().reset_index()
    stack.columns = ['gene1', 'gene2', 'score']
    stack = stack[stack.gene1 != stack.gene2]
    df = stack[stack.gene1.astype(str) < stack.gene2.astype(str)]
    threshold =  threshold_number * stack['score'].std() + stack['score'].mean() 
    return df, corr, threshold


def filter_corr_matrix(df, rows, columns):
    # delete recurrent scores i.e. ABC XYZ, XYZ ABC
    df = df.query(f'gene1 in {rows} and gene2 in {columns} or gene1 in {columns} and gene2 in {rows}').sort_values(by=['score'])
    df = df[df.gene1.astype(str) < df.gene2.astype(str)]
    return df

def read_raw_data(organs, min_cell_number):
    gene_counts = {} 
    min_gene_list = {}
    dfs = {} 

    for organ in organs:
        dfs[organ] = pd.read_csv(f'data/{organ}/Data.csv',index_col=0)
        # ham veri icerisinde hangi genlerin kac farkli hucrede eksprese edildigini sayiyoruz.
        # bu veriyi en az N sayida hucrede eksprese edilmis veri seklinde kullanacagiz.
        gene_counts[organ] = dfs[organ][dfs[organ]!=0].count(axis=1).to_frame(name='count')
        # her bir organ icin min sayiya gore filtrelenen gen listesi (exp verisi yok)
        min_gene_list[organ] = gene_counts[organ][gene_counts[organ]['count']>min_cell_number].index.to_list()
    
    return dfs, min_gene_list


def annotate(min_gene_list, universe, deg):
    universe_list = {}
    deg_list = {}
    # her bir organ icin filtrelenen genleri tabloya atiyoruz.
    x = pd.DataFrame.from_dict(min_gene_list,'index').fillna(value=np.NaN).T
    for o in x.columns:
        # burada filtrelenen liste ile annotate edilen listeleri karsilastirip olmayanlari eliyoruz
        universe_list[o] = list(set(x[o].dropna()) & set( universe['name']))
        deg_list[o] = list(set(x[o].dropna()) & set( deg['name']))

    return universe_list, deg_list



itterlist = itertools.product([200,300], [3,4])

multiple_output = {}
for i in itterlist:

    min_cell_number, threshold_number = i[0], i[1]
    dfs, min_gene_list = read_raw_data(organs, min_cell_number)
    universe_list, deg_list = annotate(min_gene_list, universe, deg)

    # her bir organ icin
    # 3 farkli matrix secerek corr yapacagiz.
    for o in organs:

        #if o is 'brain':

        # hangi genlerin kullanilacgini belirtiyoruz.
        cases = {
            'DEG_vs_nonDEG': [deg_list[o], list(set(universe_list[o])-set(deg_list[o])) ],
            'nonDEG_vs_nonDEG': [list(set(universe_list[o])-set(deg_list[o])), 
            list(set(universe_list[o])-set(deg_list[o]))],
            'DEG_vs_DEG': [ deg_list[o], deg_list[o] ],

        }

        # hangi organ icin corr yapilacaksa o organa ait ham veriyi cekiyoruz,
        # icinde 0 olan ama sadece istedigimiz genler var.
        df = dfs[o].loc[ universe_list[o] ].T

        # corr yapiyoruz.
        corr_df, corr_matrix, threshold = do_corr(df, o, threshold_number)
        #corr_matrix.to_csv(f'results/corr/{o}.csv')

        m = f'{o}_{min_cell_number}_{threshold_number}'
        single_output = {}
        for k,v in cases.items():
            fcorr_df = filter_corr_matrix(corr_df, v[0], v[1])
            fcorr_thresholded = fcorr_df.query('score > @threshold')

            abs_edge = fcorr_df.index.size
            signal_edge = fcorr_thresholded.index.size
            density = signal_edge/abs_edge
            single_output[k] = {
                'abs_edge': abs_edge,
                'signal_edge': signal_edge,
                'density':   density,
                'threshold': threshold,
                'abs_mean' : fcorr_df['score'].mean(),
                'abs_median' : fcorr_df['score'].median(),
                'abs_std' : fcorr_df['score'].std(),
            }

        multiple_output[m] = single_output


t = pd.DataFrame.from_dict(multiple_output, 'index')
for i in ['DEG_vs_DEG', 'nonDEG_vs_nonDEG', 'DEG_vs_nonDEG']:
    exceldfs = t[i].values
    exceldfnames = t.index
    startrow = 0
    with pd.ExcelWriter(f'{i}.xlsx') as writer:
        for k, exceldf in enumerate(exceldfs):
            exceldf = pd.DataFrame([exceldf])
            exceldf = exceldf.rename(index={0: exceldfnames[k]})
            if k is 0:
                exceldf.to_excel(writer, engine="xlsxwriter", startrow=startrow)
                startrow += (exceldf.shape[0] + 1)
            else:
                exceldf.to_excel(writer, engine="xlsxwriter", startrow=startrow, header=False)
                startrow += (exceldf.shape[0] )

