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
