# Funçõs exportadas para o cálculo do raio de giro de polímeros
export calc_rg, plot_rg, sv_rg

# Definir o arquivo de log
const LOG_FILE = "debug_Gyrayion_log.txt"

# Função para registrar mensagens no log
function log_message(message)
    open(LOG_FILE, "a") do file
        write(file, message * "\n")
    end
end

function check_file(file_path::AbstractString)
    if !isfile(file_path)
        error("Arquivo não encontrado: $file_path")
    end
end

# /home/m3g/.julia/dev/DTools

# Função para obter todos os números de resíduos (monômeros) presentes no polímero
function get_all_monomers(poli_data)
    # Obter todos os números de resíduos (monômeros) presentes no polímero
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

@testitem "RG" begin

    @test 1 == 1

end

# Função para salvar os dados do raio de giro em um arquivo .dat
function sv_rg(Rg_list, filename::String)
    open(filename, "w") do file
        for (i, Rg) in enumerate(Rg_list)
            write(file, "$i $Rg\n")
        end
    end
    println("Dados do raio de giro salvos em '$filename'")
end

function plot_rg(Rg_list, filename::String)
    
    plot_font = "Computer Modern" 
    default(fontfamily=plot_font, framestyle=:box)
    
    # Preparar dados para o plot
    data = hcat(collect(1:length(Rg_list)), Rg_list)
    
    # Tamanho a fonte dos eixos
    fsize = 12

    # Gráfico de Radius of Gyration vs Frames (p1)
    p1 = plot(data[:, 1] / 100, data[:, 2],
        label="Radius of Gyration",
        line=:solid,
        linecolor=:red,
        xtickfontsize=fsize,
        ytickfontsize=fsize,
        xlabelfontsize=fsize,
        ylabelfontsize=fsize,
        title="Radius of Gyration vs Frames",
        xlabel="Frames/100",
        ylabel="Radius of gyration (Å)",
        xticks=0.0:20.0:(maximum(data[:, 1] / 100) + 5),
        yticks=10:5.0:(maximum(data[:, 2]) + 5),
        xlim=[0.0, maximum(data[:, 1] / 100) + 1],
        ylim=[10.0, maximum(data[:, 2]) + 5],
        framestyle=:box,
        grid=true,
        legend=:topright
    )

    # Gráfico KDE (p2)
    kde = KernelDensity.kde(data[:, 2])
    x = range(minimum(data[:, 2]) - 5, stop=maximum(data[:, 2]) + 5, length=length(data[:, 2]))
    p2 = plot(x, z -> pdf(kde, z),
        label="KDE",
        lw = 2,
        line=:solid,
        linecolor=:red,
        xtickfontsize=fsize,
        ytickfontsize=fsize,
        xlabelfontsize=fsize,
        ylabelfontsize=fsize,
        title="KDE of Radius of Gyration",
        xlabel="Radius of gyration (Å)",
        ylabel="Probability density",
        xticks=(minimum(data[:, 2]) - 5):10.0:(maximum(data[:, 2]) + 5),
        ylim=[0.0, maximum(pdf(kde, x)) + 0.04],
        framestyle=:box,
        grid=true
    )

    # Cálculo automático das coordenadas para a anotação no p1
    # Os limites são definidos conforme os parâmetros xlim e ylim do p1
    xlims1 = (0.0, maximum(data[:, 1] / 100) + 1)
    ylims1 = (10.0, maximum(data[:, 2]) + 5)
    xA = xlims1[1] + 0.05 * (xlims1[2] - xlims1[1])
    yA = ylims1[2] - 0.05 * (ylims1[2] - ylims1[1])
    annotate!(p1, (xA, yA, text("A", (fsize + 8), :black, plot_font)))

    # Cálculo automático das coordenadas para a anotação no p2
    # Usamos os limites definidos pelo range x e os limites do eixo y
    xlims2 = (minimum(x), maximum(x))
    ylims2 = (0.0, maximum(pdf(kde, x)) + 0.04)
    xB = xlims2[1] + 0.05 * (xlims2[2] - xlims2[1])
    yB = ylims2[2] - 0.05 * (ylims2[2] - ylims2[1])
    annotate!(p2, (xB, yB, text("B", (fsize + 8), :black, plot_font)))

    # Combina os dois plots lado a lado com as configurações de layout
    combined = plot(p1, p2,
        layout=(1, 2),
        size=(1200, 600),
        legend=false,
        margin=6Plots.Measures.mm,
        dpi=400
    )

    savefig(combined, "layout_$(filename).png")
end

function calc_rg(pdb_file::String, traj_file::String; num_frames::Union{Nothing,Int}=nothing, seg_select::String)
    log_message("Carregando arquivos PDB: $pdb_file e Trajetória: $traj_file")
    simulation = Simulation(pdb_file, traj_file)
    first_frame!(simulation)
    log_message("Simulação inicializada com sucesso.")

    # Se num_frames não foi passado, usa todos os frames da simulação
    if isnothing(num_frames)
        num_frames = length(simulation)
    end
    
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