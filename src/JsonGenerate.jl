#import Pkg ; Pkg.activate("/home/m3g/Documentos/polyisobutileno/mddf_kb_1000ns/enviroments/Json")
#using PDBTools, ComplexMixtures, CSV, DataFrames, Statistics, Dates

export d_load, d_name

s_vect_s = [" ", " "]

function d_load(pdb_file::String, trajectory_file::String,)
    if !isfile(pdb_file)
        error("Arquivo PDB não encontrado: $pdb_file")
    end
    if !isfile(trajectory_file)
        error("Arquivo de trajetória não encontrado: $trajectory_file")
    end
    println("Carregando arquivos:")
    println("  Arquivo PDB: $pdb_file")
    println("  Arquivo de trajetória: $trajectory_file")
    (s_vect_s[1] = pdb_file)
    (s_vect_s[2] = trajectory_file)
    return s_vect_s

end

function d_name(replica_temp::String, simulation_time::String, bulk::String, 
    output_name::String)
    println("Gerando nome do arquivo de saída:")
    println("  Arquivo de saída: coc$(replica_temp)_$(output_name)_$(simulation_time)_$(bulk).json")
    return "coc$(replica_temp)_$(output_name)_$(simulation_time)_$(bulk).json"
end


function json_generate(; 
    pdb_file::String, 
    trajectory_file::String, 
    replica_temp::String, 
    simulation_time::String, 
    bulk::String, 
    output_name::String)
    
    # Verifica se os arquivos fornecidos existem
    if !isfile(pdb_file)
        error("Arquivo PDB não encontrado: $pdb_file")
    end
    if !isfile(trajectory_file)
        error("Arquivo de trajetória não encontrado: $trajectory_file")
    end

    println("Carregando arquivos:")
    println("  Arquivo PDB: $pdb_file")
    println("  Arquivo de trajetória: $trajectory_file")

    # Carrega o arquivo PDB
    atoms = readPDB(pdb_file)

    # Seleciona o polímero e os solventes
    polymer = PDBTools.select(atoms, "segname POLI")
    coc = PDBTools.select(atoms, "resname COC")

    # Configura o soluto (polímero)
    solute = AtomSelection(polymer, nmols=1)

    # Configura o solvente (número de átomos: oct = 26, coc = 24)
    solvent_coc = AtomSelection(coc, natomspermol=24)

    # Configura a estrutura de trajetória
    trajectory = Trajectory(trajectory_file, solute, solvent_coc)

    # Executa o cálculo de mddf e salva os resultados
    results = mddf(trajectory, Options(bulk_range=(15.0, 20.0)))
    output_file = "files_json/coc$(replica_temp)_$(output_name)_$(simulation_time)_$(bulk).json"
    save(results, output_file)

    println("Resultados salvos em: $output_file")
end










































































































































#==

# Argumentos fornecidos pelo Bash
pdb_file        = ARGS[1]  # Caminho absoluto ou relativo para o arquivo PDB
trajectory_file = ARGS[2]  # Caminho absoluto ou relativo para o arquivo de trajetória
replica_temp    = ARGS[3]  # Replica e temperatura
simulation_time = ARGS[4]  # Tempo de simulação
bulk            = ARGS[5]  # bulk1520
output_name     = ARGS[6]  # Nome do arquivo de saída

function json_generate(pdb_file::String, trajectory_file::String, replica_temp::String, simulation_time::String, bulk::String, output_name::String)
    # Verifica se os arquivos fornecidos existem
    if !isfile(pdb_file)
        error("Arquivo PDB não encontrado: $pdb_file")
    end
    if !isfile(trajectory_file)
        error("Arquivo de trajetória não encontrado: $trajectory_file")
    end

    println("Carregando arquivos:")
    println("  Arquivo PDB: $pdb_file")
    println("  Arquivo de trajetória: $trajectory_file")

    # Load PDB file of the system
    atoms = readPDB(pdb_file)

    # Select the polymer and the solvents
    polymer = PDBTools.select(atoms, "segname POLI")
    coc = PDBTools.select(atoms, "resname COC")

    # Setup solute (polymer)
    solute = AtomSelection(polymer, nmols=1)

    # Setup solvent (number of atoms of oct = 26, coc = 24)
    solvent_coc = AtomSelection(coc, natomspermol=24)

    # Setup the Trajectory structure
    trajectory = Trajectory(trajectory_file, solute, solvent_coc)

    # Run mddf calculation and save results
    results = mddf(trajectory, Options(bulk_range=(15.0, 20.0)))
    output_file = "files_json/coc$(replica_temp)_$(output_name)_$(simulation_time)_$(bulk).json"
    save(results, output_file)

    println("Resultados salvos em: $output_file")
end
# Verifica se o número correto de argumentos foi fornecido
if length(ARGS) != 6
    println("Uso: julia JsonGenerate.jl <pdb_file> <trajectory_file> <replica_temp> <simulation_time> <bulk> <output_name>")
    exit(1)
end


# Verificar se os arquivos fornecidos existem
if !isfile(pdb_file)
    error("Arquivo PDB não encontrado: $pdb_file")
end
if !isfile(trajectory_file)
    error("Arquivo de trajetória não encontrado: $trajectory_file")
end

println("Carregando arquivos:")
println("  Arquivo PDB: $pdb_file")
println("  Arquivo de trajetória: $trajectory_file")

# Load PDB file of the system
atoms = readPDB(pdb_file)

# Select the polymer and the solvents
polymer = PDBTools.select(atoms, "segname POLI")
coc = PDBTools.select(atoms, "resname COC")

# Setup solute (polymer)
solute = AtomSelection(polymer, nmols=1)

# Setup solvent (number of atoms of oct = 26, coc = 24)
solvent_coc = AtomSelection(coc, natomspermol=24)

# Setup the Trajectory structure
trajectory = Trajectory(trajectory_file, solute, solvent_coc)

# Run mddf calculation and save results
results = mddf(trajectory, Options(bulk_range=(15.0, 20.0)))
output_file = "files_json/coc$(replica_temp)_$(output_name)_$(simulation_time)_$(bulk).json"
save(results, output_file)

println("Resultados salvos em: $output_file")

==#