#!/usr/bin/python
"""Print boxplots of co-occurences in targets file.

python $HOME/gsea_target_dependency_distributions/script.py gsea_fname=$HOME/gsea_target_dependency_distributions/c3.tft.v3.0.symbols.gmt dependency_json=$HOME/gse2034/gse2034_best.json tabfile=$HOME/gse2034/GSE2034.GPL96.eQTL.normed.tab outdir=$HOME/gse2034/gsea_boxplots_all
"""
import json
from pprint import pprint
import os
import sys
from util import *
from py_symmetric_matrix import *
import numpy as np
import matplotlib
from scipy import stats
matplotlib.use('agg') # required for use on OSC servers
import matplotlib.pyplot as plt


# per Dr. Kun
GENE_LIST = ["E2F1", "ETS2", "ESR1"]
RX_SYM = re.compile("[^$]*\$([^_]+)_(.*)") # (gene name, target id)


def main(gsea_fname=None, dependency_json=None, tabfile=None, outdir=None, list_all=False):
  assert gsea_fname and dependency_json and tabfile and outdir

  # Load list of variables in study data file.
  varlist = []
  for row in name_iter(open(tabfile), varlist): pass
  cleaned_varlist = [clean(s) for s in varlist]
  cleaned_varset = set(cleaned_varlist) # for lookups

  # Load targets from gsea file
  targets = {} # {gene: {id: set(genes)}}
  gsea_gene_set = set()
  for line in open(gsea_fname):
    row = line.split('\t')
    m = RX_SYM.match(row[0])
    if not m:
      continue
    gene, target_group = clean(m.group(1)), m.group(2)
    gsea_gene_set.add(gene)
    if gene in cleaned_varset:
      d = targets.setdefault(gene, {})
      d[target_group] = set(map(clean, row[2:]))
      
  print "Loaded %d genes; kept %d of these in study." % (len(gsea_gene_set), len(targets))

  # Load dependencies
  D = json.load(open(dependency_json))
  for dep_name, d in D["dependencies"].items():
    print "Loading %s..." % dep_name
    M = np.load(os.path.join(d['dir'], d['values_file']))
    B = np.load(os.path.join(d['dir'], d['bool_file']))
    d['Q'] = NamedSymmetricMatrix(store_diagonal=False, dtype=np.float, var_list=cleaned_varlist, matrix=M)
    d['Q_Mask'] = NamedSymmetricMatrix(store_diagonal=False, dtype=np.bool, var_list=cleaned_varlist, matrix=B)


  if not list_all:
    LIST = GENE_LIST
  else:
    LIST = targets.keys()
  R = {}
    
  for gene in LIST:

    if gene not in cleaned_varset:
      print "WARNING: %s not in cleaned_varset." % gene
      continue
    if gene not in targets:
      print "WARNING: %s not in targets." % gene
      continue
    R[gene] = {}

    for dep_name, d in D["dependencies"].items():
      R[gene][dep_name] = {}
      Q, Q_Mask = d['Q'], d['Q_Mask']

      # Get all other scores
      cleaned_varset.remove(gene) # remove self from list of genes to consider
      all_scores = [ \
        Q.get(gene,q) for q in cleaned_varset if \
          q != gene and Q_Mask.get(gene, q)]
      cleaned_varset.add(gene) # remove self from list of genes to consider

      # Get per-target-set scores
      for target_set_name, target_set in targets[gene].items():
        # this could be done by views if performance were an issue
        shared_targets = target_set & cleaned_varset - set([gene])
        target_scores = [ \
          Q.get(gene,q) for q in shared_targets if \
             Q_Mask.get(gene, q)]
        ptest = stats.ranksums(target_scores, all_scores)

        R[gene][dep_name] = {'z': ptest[0], 'p': ptest[1]}
        R[gene][dep_name]['target'] = compile_stats(target_scores)
        R[gene][dep_name]['all'] = compile_stats(all_scores)
        print "Scores for %s, set %s." % (gene, target_set_name)
        print "Ranksum ptest (Wilcoxon): z=%f p=%f" % ptest
        print "%d target interactions. %d targets in study considered." % (len(target_set), len(shared_targets))
        print "targets only: "; pprint(R[gene][dep_name]['target'])
        print "all genes: "; pprint(R[gene][dep_name]['all'])
  
        # BOXPLOTS
        id_set = (gene, target_set_name, dep_name)
        plt.clf(); plt.cla()
        plt.boxplot((target_scores, all_scores))
        plt.title(" ".join(id_set) + " targets vs all %.6f" % ptest[1]*100)
        fname = os.path.join(outdir,"%s.png" % ("_".join(id_set)))
        plt.savefig(fname)
        R[gene][dep_name]['plot_fname'] = fname
        print "saved boxplot %s" % fname
        print 

  json.dump(R, open(os.path.join(outdir, "results.json"), "w"))

def compile_stats(a):
  return {
    'mean': np.mean(a),
    'std': np.std(a),
    'min': np.min(a),
    'max': np.max(a),
    }
      
if __name__ == "__main__":
  main(**dict([s.split("=") for s in sys.argv[1:]])) 
