#!/usr/bin/python
"""Print boxplots of co-occurences in targets file.

python $HOME/gsea_target_dependency_distributions/script.py gsea_fname=$HOME/gsea_target_dependency_distributions/c3.tft.v3.0.symbols.gmt dependency_json=$HOME/gse2034/gse2034_promising.json tabfile=$HOME/gse2034/GSE2034.GPL96.eQTL.normed.tab outdir=$HOME/gse2034/gsea_boxplots

WARNING: THIS IS VERY INEFFICIENT: USE INDEXED ADJ MATRICES
"""
import json
import random
import os
import sys
from util import *
from py_symmetric_matrix import *
import numpy as np
import itertools
import matplotlib 
matplotlib.use('agg') # required for use on OSC servers
import matplotlib.pyplot as plt


# per Dr. Kun
GENE_LIST = ["E2F1", "ETS2", "ESR1"]


def main(gsea_fname=None, dependency_json=None, tabfile=None, outdir=None):
  assert gsea_fname and dependency_json and tabfile and outdir

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
  if os.path.exists("%s.B.npy" % gsea_fname):
    print "Loading '%s.B.npy'" % gsea_fname
    cached = np.load("%s.B.npy" % gsea_fname)
  else:
    cached = None
  B = NamedSymmetricMatrix(store_diagonal=False, dtype=np.bool, var_list=shared_genes_list, matrix=cached)

  if cached is None:
    fp = open(gsea_fname)
    n_set = 0
    n_dupe = 0
    for line in fp:
      row = set([clean(s) for s in line.strip('\n').split('\t')[2:]])
      row = row & shared_genes
      for x, y in itertools.combinations(row, 2):
        if not B.get(x,y):
          B.set(x,y,1)
          n_set += 1
        else:
          n_dupe += 1
    print "Set %d interactions (%d dupes)" % (n_set, n_dupe)
    np.save("%s.B.npy" % gsea_fname, B._m)
    print "Saved '%s.B.npy' % gsea_fname to file"

  # load a dependency matrix
  D = json.load(open(dependency_json))
  for dep_name, d in D["dependencies"].items():
    #if dep_name != 'pcc': continue # hack for now

    M = np.load(os.path.join(d['dir'], d['values_file']))
    Mask = np.load(os.path.join(d['dir'], d['bool_file']))

    Q = NamedSymmetricMatrix(store_diagonal=False, dtype=np.float, var_list=cleaned_varlist, matrix=M)
    Q_Mask = NamedSymmetricMatrix(store_diagonal=False, dtype=np.bool, var_list=cleaned_varlist, matrix=M)

    for gene in GENE_LIST:
      interact_scores = []
      random_scores = []
      sample_scores = []
      for other in shared_genes_list:
        if other != gene and B.get(gene, other) and Q_Mask.get(gene, other):
          interact_scores.append(Q.get(gene,other))
      n = len(interact_scores)
      # random
      for other in random.sample(shared_genes_list, n):
        if other != gene and Q_Mask.get(gene, other):
          random_scores.append(Q.get(gene,other))
      # sample
      for other in random.sample(cleaned_varlist, 1000):
        if other != gene and Q_Mask.get(gene, other):
          sample_scores.append(Q.get(gene,other))

      print "Scores for %s." % gene
      print "%d interactions. %d random from intersecting gene list. Missing values skipped." % (len(interact_scores),len(random_scores))
      print "%d samples of %d from dependency matrix." % (len(sample_scores), len(cleaned_varlist))
      print "interact:", np.mean(interact_scores), np.std(interact_scores), np.max(interact_scores)
      print "random:", np.mean(random_scores), np.std(random_scores), np.max(random_scores)
      print "sample:", np.mean(sample_scores), np.std(sample_scores), np.max(sample_scores)

      # BOXPLOTS
      plt.clf(); plt.cla()
      plt.boxplot((interact_scores, random_scores, sample_scores))
      plt.savefig(os.path.join(outdir,"%s_%s_gse2034_boxplot.png" % (gene, dep_name)))
    

  
if __name__ == "__main__":
  main(**dict([s.split("=") for s in sys.argv[1:]])) 
