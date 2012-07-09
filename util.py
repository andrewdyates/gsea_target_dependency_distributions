#!/usr/bin/python
import datetime
import re

GENE_RENAMES = {
  'MO25ALPHA': 'MO25A',
  'MO25BETA': 'MO25B',
  'STRADBETA': 'STRADB',
  'STRADALPHA': 'STRADA',
  '1DEC': 'DEC1',
  '1MAR': 'MARCH1',
  '10MAR': 'MARCH10',
  '11MAR': 'MARCH11',
  '2MAR': 'MARCH2',
  '3MAR': 'MARCH3',
  '4MAR': 'MARCH4',
  '5MAR': 'MARCH5',
  '6MAR': 'MARCH6',
  '7MAR': 'MARCH7',
  '8MAR': 'MARCH8',
  '9MAR': 'MARCH9',
  '1SEP': 'SEPT1',
  '10SEP': 'SEPT10',
  '11SEP': 'SEPT11',
  '12SEP': 'SEPT12',
  '13SEP': 'SEPT13',
  '14SEP': 'SEPT14',
  '2SEP': 'SEPT2',
  '3SEP': 'SEPT3',
  '4SEP': 'SEPT4',
  '5SEP': 'SEPT5',
  '6SEP': 'SEPT6',
  '7SEP': 'SEPT7',
  '8SEP': 'SEPT8',
  '9SEP': 'SEPT9',
}

def timestamp(date=None):
  if not date:
    date = datetime.datetime.now()
  return date.isoformat(' ')

def name_iter(fp, varlist):
  """Load a labeled row matrix as a line iterator."""
  for line in fp:
    if line[0] in ('#', '\n'): continue
    name,c,row = line.partition('\t')
    varlist.append(name)
    yield row

def clean(s):
  """Standardize gene names to be all caps, alphanumeric."""
  s = re.sub('[^a-zA-Z0-9]', '', s.upper())
  if s in GENE_RENAMES:
    s = GENE_RENAMES[s]
  return s
