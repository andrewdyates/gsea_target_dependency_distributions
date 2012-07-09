#!/usr/bin/python
"""Print boxplots of targets vs non-targets.

python script.py gsea_fname=~/Desktop/c3.tft.v3.0.symbols.gmt dependency_json=1 tabfile=/Users/qq/Dropbox/biostat/eqtl_data/GSE2034/GSE2034.GPL96.eQTL.normed.tab
"""
import json
import sys
from util import *
from py_symmetric_matrix import *
import numpy as np
import itertools

def main(gsea_fname=None, dependency_json=None, tabfile=None):
  assert gsea_fname and dependency_json and tabfile

  varlist = []
  for row in name_iter(open(tabfile), varlist): pass
  cleaned_varlist = [clean(s) for s in varlist]

  # first get names
  fp = open(gsea_fname)
  gsea_list = set()
  for line in fp:
    row = line.strip('\n').split('\t')
    for name in row[2:]:
      gsea_list.add(clean(name))
  fp.close()

  shared_genes = gsea_list & set(cleaned_varlist)
  shared_genes_list = list(shared_genes)
  print "Num genes in GSEA list", len(gsea_list)
  print "Num genes in Tab file", len(cleaned_varlist)
  print "Num of intersecting genes", len(shared_genes)

  # load pairs into graph
  B = NamedSymmetricMatrix(store_diagonal=False, dtype=np.bool, var_list=shared_genes_list)

  fp = open(gsea_fname)
  for line in fp:
    row = set([clean(s) for s in line.strip('\n').split('\t')[2:]])
    row = row & shared_genes
    for x, y in itertools.combinations(row, 2):
      B.set(x,y,1)

  # load a dependency matrix
  #D = json.load(open(dependency_json))
  #D["dependencies"]["mic"]
  

  
if __name__ == "__main__":
  main(**dict([s.split("=") for s in sys.argv[1:]])) 
