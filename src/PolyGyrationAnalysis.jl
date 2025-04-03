module PolyGyrationAnalysis
export run_analysis, plot_radius_of_gyration, save_radius_of_gyration
using PDBTools
using MolSimToolkit
using StaticArrays
using ProgressMeter
using LinearAlgebra
# using Statistics
using Plots 
using KernelDensity

    # Definir o arquivo de log
    const LOG_FILE = "debug_log.txt"

    # Função para registrar mensagens no log
    function log_message(message)
        open(LOG_FILE, "a") do file
            write(file, message * "\n")
        end
    end

    # Função para obter todos os números de resíduos (monômeros) presentes no polímero
    function get_all_monomers(poli_data)
        return unique(atom.resnum for atom in poli_data)
    end

    # Função para calcular o centro de massa
    function center_of_mass(simulation, indices::AbstractVector{Int}, coor::Vector{SVector{3, Float64}})
        totmass = 0.0
        cm = SVector{3, Float64}(0.0, 0.0, 0.0)
        for i in indices
            m = atomic_mass(atoms(simulation)[i])
            cm += coor[i] * m
            totmass += m
        end
        return cm / totmass
    end

    # Função para calcular o raio de giro
    function radius_of_gyration(simulation, indices::AbstractVector{Int}, coor::Vector{SVector{3, Float64}})
        cm = center_of_mass(simulation, indices, coor)
        dist_squared_sum = 0.0
        total_mass = 0.0
        for i in indices
            m = atomic_mass(atoms(simulation)[i])
            dist_squared_sum += m * norm(coor[i] - cm)^2
            total_mass += m
        end
        return sqrt(dist_squared_sum / total_mass)
    end

    # Função para converter posições do frame para SVectors
    function convert_positions_to_svector(positions::FramePositions)
        return [SVector{3, Float64}(pos.x, pos.y, pos.z) for pos in positions]
    end

    # Função para calcular o raio de giro para cada frame
    function calculate_radius_of_gyration_per_frame(simulation, monomer_indices, num_frames)
        Rg_list = Float64[]
        p = Progress(num_frames, desc="Calculando raio de giro")
        for (iframe, frame) in enumerate(simulation)
            if iframe > num_frames
                break
            end
            next!(p)
            coor_raw = positions(frame)
            coor = convert_positions_to_svector(coor_raw)
            # Coletar todos os índices dos monômeros selecionados
            all_monomer_indices = vcat(values(monomer_indices)...)
            Rg = radius_of_gyration(simulation, all_monomer_indices, coor)
            #println("Frame $iframe: Raio de Giro = $Rg Å")
            push!(Rg_list, Rg)
        end
        return Rg_list
    end

    # Função para salvar os dados do raio de giro em um arquivo .dat
    function save_radius_of_gyration(Rg_list, filename::String="radius_of_gyration.dat")
        open(filename, "w") do file
            for (i, Rg) in enumerate(Rg_list)
                write(file, "$i $Rg\n")
            end
        end
        println("Dados do raio de giro salvos em '$filename'")
    end

    # Função para plotar Radius of Gyration e KDE
    function plot_radius_of_gyration(Rg_list)
        # Preparar dados para o plot
        data = hcat(collect(1:length(Rg_list)), Rg_list)
        
        # Gráfico de Radius of Gyration vs Frames
        p1 = plot(data[:, 1] / 100, data[:, 2],
            label="Radius of Gyration",
            line=:solid,
            linecolor=:blue,
            title="Radius of Gyration vs Frames",
            xlabel="Frames/100",
            ylabel="Radius of gyration (Å)",
            xticks=0.0:5.0:(maximum(data[:, 1] / 100) + 5),
            yticks=10:5.0:(maximum(data[:, 2]) + 5),
            xlim=[0.0, maximum(data[:, 1] / 100) + 5],
            ylim=[10.0, maximum(data[:, 2]) + 5],
            framestyle=:box,
            grid=false,
            legend=:topright
        )
        savefig(p1, "plot_radius_of_gyration_meu_script.png")
        println("Gráfico de Radius of Gyration salvo como 'plot_radius_of_gyration_meu_script.png'")
        
        # Gráfico KDE
        kde = KernelDensity.kde(data[:, 2])
        x = range(minimum(data[:, 2]) - 5, stop=maximum(data[:, 2]) + 5, length=length(data[:, 2]))
        
        p2 = plot(x, z -> pdf(kde, z),
            label="KDE",
            line=:solid,
            linecolor=:red,
            title="KDE of Radius of Gyration",
            xlabel="Radius of gyration (Å)",
            ylabel="Probability density",
            xticks=(minimum(data[:, 2]) - 5):5.0:(maximum(data[:, 2]) + 5),
            ylim=[0.0, maximum(pdf(kde, x)) + 0.04],
            framestyle=:box,
            grid=true
        )
        savefig(p2, "plot_kde_radius_of_gyration.png")
        println("Gráfico KDE salvo como 'plot_kde_radius_of_gyration.png'")
    end

# Função que realiza toda a análise
# Parâmetros:
#   pdb_file   : nome do arquivo PDB
#   traj_file  : nome do arquivo da trajetória
#   num_frames : número de frames a analisar (padrão: 15000)
#   seg_select : string para seleção dos átomos, por exemplo, "segname POLI" (padrão)

    function run_analysis(pdb_file::String, traj_file::String, num_frames::Int=15000, seg_select::String="segname POLI")
        log_message("Carregando arquivos PDB: $pdb_file e Trajetória: $traj_file")
        simulation = Simulation(pdb_file, traj_file)
        first_frame!(simulation)
        log_message("Simulação inicializada com sucesso.")
        
        # Ler os dados do PDB e selecionar os átomos de interesse
        pdb_data = readPDB(pdb_file)
        poli_data = select(pdb_data, seg_select)
        
        # Selecionar todos os monômeros presentes na seleção
        monomer_targets = get_all_monomers(poli_data)
        log_message("Nenhum monômero fornecido. Usando todos os monômeros: $monomer_targets")
        
        # Criar dicionário com os índices dos átomos para cada monômero
        monomer_indices = Dict{Int, Vector{Int}}()
        for resnum in monomer_targets
            atom_indices = findall(atom -> atom.resnum == resnum, poli_data)
            monomer_indices[resnum] = atom_indices
        end
        
        log_message("Iniciando cálculo do raio de giro.")
        Rg_list = calculate_radius_of_gyration_per_frame(simulation, monomer_indices, num_frames)
        log_message("Cálculo do raio de giro concluído.")
        return Rg_list
    end

end # module