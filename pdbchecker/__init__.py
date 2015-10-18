from ftplib import FTP
import tempfile
import gzip
import MDAnalysis as mda
import json
import numpy
import numpy as np

THRESHOLD = 4.5   # Threshold for identifiction of a gap in residues

def get_PDB_universe(pdbcode):
  # Get file from PDB
  filename = 'pdb'+str(pdbcode)+'.ent.gz'
  ftp = FTP('ftp.wwpdb.org')
  ftp.login()
  ftp.cwd('pub/pdb/data/structures/all/pdb')
  gzipfile = tempfile.NamedTemporaryFile(suffix='.gz')
  ftp.retrbinary('RETR {}'.format(filename), gzipfile.write)
  ftp.quit()

  # unzip PDB file
  gzipfile.seek(0)
  with gzip.GzipFile(pdbcode, 'rb', 0, gzipfile) as unzipper:
    pdbcontent = unzipper.read()
  gzipfile.close()

  with open('temp.pdb', mode='w') as test:
    test.write(pdbcontent)

  # Is there no way to create a Universe directly from a text object?
  fh = tempfile.NamedTemporaryFile(suffix='.pdb', mode='w')
  fh.write(pdbcontent)
  fh.flush()

  # create universe
  u = mda.Universe(fh.name)

  # clean up and return
  fh.close()
  return u

def get_apolar(u):
    types = json.load(open('types.json'))
    selector = ' or '.join(['(resname {} and name {})'.format(pair[0], pair[1]) for pair in types['Apolar']])
    apolar = u.atoms.select_atoms(selector)
    return apolar

def get_polar_mixt(u):
    types = json.load(open('types.json'))
    selector = ' or '.join(['(resname {} and name {})'.format(pair[0], pair[1]) for pair in types['Polar Mixt']])
    selection = u.atoms.select_atoms(selector)
    return selection
def get_polar_donor_ionic_positive(u):
    types = json.load(open('types.json'))
    selector = ' or '.join(['(resname {} and name {})'.format(pair[0], pair[1]) for pair in types['Polar Donnor, Ionic Positive']])
    selection = u.atoms.select_atoms(selector)
    return selection
def get_polar_acceptor_ionic_negative(u):
    types = json.load(open('types.json'))
    selector = ' or '.join(['(resname {} and name {})'.format(pair[0], pair[1]) for pair in types['Polar Acceptor, Ionic Negative']])
    selection = u.atoms.select_atoms(selector)
    return selection
def get_polar_donor(u):
    types = json.load(open('types.json'))
    selector = ' or '.join(['(resname {} and name {})'.format(pair[0], pair[1]) for pair in types['Polar Donor']])
    selection = u.atoms.select_atoms(selector)
    return selection
def get_polar_acceptor(u):
    types = json.load(open('types.json'))
    selector = ' or '.join(['(resname {} and name {})'.format(pair[0], pair[1]) for pair in types['Polar Acceptor']])
    selection = u.atoms.select_atoms(selector)
    return selection

def get_gaps(u):
  ca = u.atoms.select_atoms('name CA')

  for c1 in ca.atoms:
      found = False
      for c2 in ca.atoms:
          if c2.residue.resnum-c1.residue.resnum > 0 and not(found):
              found = True
              p1 = c1.position
              p2 = c2.position
              dist_ca = np.linalg.norm(p1-p2)
              if (dist_ca > THRESHOLD):
                  print "Residue Gap between residues:", c1.residue.resnum, c2.residue.resnum, dist_ca
              
