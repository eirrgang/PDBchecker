# PDBchecker
## Get a universe for a PDB file from the web database
    import pdbchecker
    u = pdbchecker.get_PDB_universe('1aki')

## Identify atoms of certain types
    pdbchecker.get_apolar(u)
    pdbchecker.get_polar_acceptor(u)
    pdbchecker.get_polar_acceptor_ionic_negative(u)
    pdbchecker.get_polar_donor(u)
    pdbchecker.get_polar_donor_ionic_positive(u)
    pdbchecker.get_polar_mixt(u)

## Identify gaps between residues
