import bioservices
import tempfile
import MDAnalysis as mda

def get_PDB_universe(molecule):
  s = bioservices.PDB()
  res = s.get_file("1aki", 'pdb')

  fh = tempfile.NamedTemporaryFile(suffix='.pdb')
  fh.write(res)

  u = mda.Universe(fh.name)
  return u
