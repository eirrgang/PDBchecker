import bioservices
import tempfile
import MDAnalysis as mda
import json

def get_PDB_universe(molecule):
  s = bioservices.PDB()
  res = s.get_file("1aki", 'pdb')

  fh = tempfile.NamedTemporaryFile(suffix='.pdb')
  fh.write(res)

  u = mda.Universe(fh.name)
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
