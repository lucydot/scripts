#!/Applications/Julia-0.5.app/Contents/Resources/julia/bin/julia

using JuliaPhonons
using JLD
using PyCall

# read data from POSCAR and mesh.yaml using JuliaPhonons module. If already read then load from file.
if isfile("jld_phonon_variables.jld") == false
  println("Reading poscar and mesh.yaml.....")
  poscar_data = read_POSCAR(open("POSCAR"))
  mesh = read_meshyaml(open("mesh.yaml"),poscar_data)
  println("Saving variables as .jld......") # save variables for further analysis with JuliaPhonons
  save("./jld_phonon_variables.jld","mesh",mesh,"poscar_data",poscar_data)
else
  println("loading poscar and mesh.yaml data.....")
  poscar_data = load("jld_phonon_variables.jld")["poscar_data"]
  mesh = load("jld_phonon_variables.jld")["mesh"]
end

# Calculate IPR data from every phonon mode. Store data in dict IPR_dictionary[mode][key1/atom][key2]
IPR_dictionary = Dict()
println("calculating IPR data for all phonon modes.....")
for x = 1:length(mesh[1]) #iterate through all phonon modes 
    mode_dictionary = Dict()
    mode_data = decompose_eigenmode_atom_contributions(poscar_data, mesh[2][x], mesh[1][x])
    mode_dictionary["IPR-energy"] = mode_data[1]
    mode_dictionary["IPR-distance"] = mode_data[2]
    mode_dictionary["sume"] = mode_data[3]
    mode_dictionary["sumd"] = mode_data[4]
    mode_dictionary["freq"] = mode_data[6]
    IPR_dictionary[x] = merge(mode_dictionary,mode_data[5]) #mode_data[5] is a dictionary whose keys correspond to atom number. This contains another nested dictionary with keys "efrac","dfrac" and "PRE"
end

println("converting julia dictionary to python dictionary.....")
# Use Pycall to convert Julia dict into Python dict
@pyimport pickle
filehandle = open("IPR_dictionary.p","w")
pickle.dump(IPR_dictionary,filehandle)


