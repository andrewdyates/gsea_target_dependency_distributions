#!/usr/bin/python
"""Print boxplots of co-occurences in targets file.

python $HOME/gsea_target_dependency_distributions/script.py gsea_fname=$HOME/gsea_target_dependency_distributions/c3.tft.v3.0.symbols.gmt dependency_json=$HOME/gse2034/gse2034_promising.json tabfile=$HOME/gse2034/GSE2034.GPL96.eQTL.normed.tab outdir=$HOME/gse2034/gsea_boxplots
"""
import json
import os
import sys
from util import *
from py_symmetric_matrix import *
import numpy as np
import matplotlib 
matplotlib.use('agg') # required for use on OSC servers
import matplotlib.pyplot as plt


# per Dr. Kun
GENE_LIST = ["E2F1", "ETS2", "ESR1"]
RX_SYM = re.compile("[^$]*\$([^_]+)_(.*)") # (gene name, target id)


def main(gsea_fname=None, dependency_json=None, tabfile=None, outdir=None):
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
    print gene, target_group
    d = targets.setdefault('gene', {})
    d[target_group] = map(clean, row[2:]) 
    gsea_gene_set.add(gene)

  shared_set = cleaned_varset & gsea_gene_set
  print "Loaded %d genes with targets, %d of these in study." % (len(gsea_gene_set), len(shared_set))

  # Load dependencies
  D = json.load(open(dependency_json))
  for dep_name, d in D["dependencies"].items():
    print "Loading %s..." % dep_name
    M = np.load(os.path.join(d['dir'], d['values_file']))
    B = np.load(os.path.join(d['dir'], d['bool_file']))
    D["dependencies"]['Q'] = NamedSymmetricMatrix(store_diagonal=False, dtype=np.float, var_list=cleaned_varlist, matrix=M)
    D["dependencies"]['Q_Mask'] = NamedSymmetricMatrix(store_diagonal=False, dtype=np.bool, var_list=cleaned_varlist, matrix=B)


  for gsea_gene in GENE_LIST:

    if gsea_gene not in shared_set:
      print "WARNING: gsea_gene not in shared set! Skipping..."
      continue

    for dep_name, d in D["dependencies"].items():
      Q, Q_Mask = d['Q'], d['Q_Mask']

      # Get all other scores
      all_scores = np.zeros(len(cleaned_varlist))
      for i, q in enumerate(cleaned_varset):
        if q != gsea_gene and Q_Mask.get(gsea_gene, q):
          all_scores[i] = Q.get(gsea_gene, q)

      # Get per-target-set scores
      for target_set_name, target_list in targets[gsea_gene].items():
        shared_targets = shared_set - set(target_list)
        idxs = [Q.get_idx(q, gsea_gene) for q in shared_targets if Q_Mask.get(q, gsea_gene)]
        target_scores = all_scores.take(idxs)
  
        print "Scores for %s, set %s." % (gsea_gene, target_set_name)
        
        print "%d target interactions. %d targets in study considered." % (len(target_list), len(shared_targets))
        print "targets only: ", compile_stats(target_scores)
        print "all genes: ", compile_stats(all_scores)
  
        # BOXPLOTS
        id_set = (gsea_gene, target_set_name, dep_name)
        plt.clf(); plt.cla()
        plt.boxplot((target_scores, all_scores))
        plt.title(" ".join(id_set) + " targets vs all")
        fname = os.path.join(outdir,"%s.png" % "_".join(id_set))
        plt.savefig(fname)
        print "saved boxplot %s" % fname
        print "target scores: ", target_scores
        print "first 20 all scores: ", all_scores[:20]
        print 
    


def compile_stats(a):
  return {
    'mean': np.mean(a),
    'std': np.std(a),
    'min': np.min(a),
    'max': np.max(a),
    }
      
if __name__ == "__main__":
  main(**dict([s.split("=") for s in sys.argv[1:]])) 
