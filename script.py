#!/usr/bin/python
"""Print boxplots of targets vs non-targets.

python script.py gsea_fname=~/Desktop/c3.tft.v3.0.symbols.gmt dependency_json=1 tabfile=/Users/qq/Dropbox/biostat/eqtl_data/GSE2034/GSE2034.GPL96.eQTL.tab
"""

import sys
from util import *
from py_symmetric_matrix import *

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
    
  print len(gsea_list)

  
if __name__ == "__main__":
  main(**dict([s.split("=") for s in sys.argv[1:]])) 
