import pymatgen.io.vasp.outputs as pivo
vr = pivo.Vasprun("vasprun.xml")
print(vr.eigenvalue_band_properties)
