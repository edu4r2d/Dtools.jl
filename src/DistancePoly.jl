# module PolymerDistance

export distancia_entre_atomospdb

using PDBTools
using LinearAlgebra
using Printf

"""
    distancia_entre_atomospdb(pdb_file::String,
                              resnum1::Int, atom_name1::String,
                              resnum2::Int, atom_name2::String)

Calcula a distância entre dois átomos especificados em um arquivo PDB e salva o resultado em um arquivo `.txt`.

# Parâmetros
- `pdb_file`: Caminho para o arquivo `.pdb`
- `resnum1`: Número do resíduo do primeiro átomo
- `atom_name1`: Nome do primeiro átomo (ex: "C1")
- `resnum2`: Número do resíduo do segundo átomo
- `atom_name2`: Nome do segundo átomo (ex: "CZ")
"""
function distancia_entre_atomospdb(pdb_file::String,
                                   resnum1::Int, atom_name1::String,
                                   resnum2::Int, atom_name2::String)

    atoms = read_pdb(pdb_file)

    idx1 = findfirst(at -> at.resnum == resnum1 && at.name == atom_name1, atoms)
    idx2 = findfirst(at -> at.resnum == resnum2 && at.name == atom_name2, atoms)

    if idx1 === nothing || idx2 === nothing
        error("Um ou ambos os átomos não foram encontrados no arquivo PDB.")
    end

    atom1 = atoms[idx1]
    atom2 = atoms[idx2]

    coord1 = [atom1.x, atom1.y, atom1.z]
    coord2 = [atom2.x, atom2.y, atom2.z]

    distancia = Float64(norm(coord1 - coord2))
    distancia_rounded = round(distancia, digits=3)

    println("Distância entre $atom_name1 (resíduo $resnum1) e $atom_name2 (resíduo $resnum2): $distancia_rounded Å")

    # Nome do arquivo de saída
    saida_txt = replace(pdb_file, r"\.pdb$" => "_distancia.txt")
    open(saida_txt, "w") do io
        println(io, "Arquivo: $(basename(pdb_file))")
        println(io, "Átomo inicial: $atom_name1 (resíduo $resnum1)")
        println(io, "Átomo final: $atom_name2 (resíduo $resnum2)")
        @printf(io, "Distância: %.3f Å\n", distancia_rounded)
    end

    println("Resultado salvo em: $saida_txt")

    return distancia_rounded
end

# end # module


#==
include("PolymerDistance.jl")
using .PolymerDistance

# Chamada da função
dist = distancia_entre_atomospdb("gamma_p100.pdb", 1, "C1", 100, "CZ")
==#