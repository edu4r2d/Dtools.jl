export calc_median

# Constante para o arquivo de log
const LOG_FILE = "debug_Median_log.txt"

# Função para registrar mensagens no log
function log_message(message::AbstractString)
    open(LOG_FILE, "a") do file
        write(file, message * "\n")
    end
end

# Função para verificar a existência de um arquivo
function check_file(file_path::AbstractString)
    if !isfile(file_path)
        error("Arquivo não encontrado: $file_path")
    end
end

# Função para carregar a simulação e selecionar os dados do polímero
function load_simulation(pdb_file::AbstractString, traj_file::AbstractString, seg_select::AbstractString)
    log_message("Carregando arquivos: PDB: $pdb_file, Trajetória: $traj_file")
    
    # Verificar existência dos arquivos
    check_file(pdb_file)
    check_file(traj_file)
    
    println("Carregando arquivos:")
    println("  PDB: $pdb_file")
    println("  Trajetória: $traj_file")
    
    # Inicializar a simulação e selecionar o primeiro frame
    simulation = Simulation(pdb_file, traj_file)
    first_frame!(simulation)
    
    # Carregar dados do PDB e selecionar os átomos do polímero
    pdb_data = readPDB(pdb_file)
    poli_data = select(pdb_data, "$seg_select")
    
    return simulation, poli_data
end

# Função para preparar o dicionário de índices dos átomos para cada monômero
function prepare_monomer_indices(poli_data)
    # Se não forem fornecidos monômeros específicos, usar todos os encontrados
    monomer_targets = unique([atom.resnum for atom in poli_data])
    log_message("Monômeros detectados: $monomer_targets")
    
    monomer_indices = Dict{Int, Vector{Int}}()
    for resnum in monomer_targets
        atom_indices = findall(atom -> atom.resnum == resnum, poli_data)
        if !isempty(atom_indices)
            monomer_indices[resnum] = atom_indices
        end
    end
    
    if isempty(monomer_indices)
        error("Nenhum monômero encontrado! Verifique a seleção.")
    end
    
    return monomer_indices
end

# Função para iterar sobre os frames e calcular as distâncias
function analyze_distances(simulation, polymer_indices, monomer_indices, num_frames::Union{Int, Nothing}=nothing)
    
    if isnothing(num_frames)
        num_frames = length(simulation)
    end
    
    distances_all_frames = []
    
    # Inicializar a barra de progresso
    p = Progress(num_frames, desc="Processando frames")
    
    for (iframe, frame) in enumerate(simulation)
        if iframe > num_frames
            break
        end
        
        next!(p)
        coor = positions(frame)
        
        # Calcular o centro de massa do polímero
        cm_polymer = PDBTools.center_of_mass(polymer_indices, simulation, coor; iref=nothing)
        
        distances_this_frame = Float64[]
        for (_, atom_indices) in monomer_indices
            # Calcular o centro de massa do monômero
            cm_residue = PDBTools.center_of_mass(atom_indices, simulation, coor; iref=nothing)
            # Calcular a distância e registrar
            distance = norm(cm_residue - cm_polymer)
            push!(distances_this_frame, distance)
        end
        
        push!(distances_all_frames, distances_this_frame)
    end
    
    return distances_all_frames
end

# ======================================================
#           ANÁLISE E PLOTAGEM
# ======================================================

# Função para calcular as medianas e plotar os resultados
function analyze_and_plot(distances_all_frames)
    medians = [isempty(frame_distances) ? NaN : median(frame_distances) for frame_distances in distances_all_frames]
    
    plt = plot(medians,
        title="Mediana das Distâncias por Frame",
        xlabel="Frames",
        ylabel="Mediana da Distância",
        legend=false)
    display(plt)
    
    return medians
end

# ======================================================
#           SALVAMENTO DOS RESULTADOS
# ======================================================

# Função para salvar os dados de mediana em um arquivo .dat
function save_to_dat_file(medians, filename::AbstractString)
    open(filename, "w") do file
        for (i, mediana) in enumerate(medians)
            write(file, "$i $mediana\n")
        end
    end
    println("Médias salvas em: $filename")
end

# ======================================================
#           FUNÇÃO PRINCIPAL EXPORTADA
# ======================================================

"""
    calc_median(pdb_file, traj_file, output_file; num_frames=150000)

Função principal que orquestra o carregamento dos dados, processamento dos frames, cálculo das medianas e salvamento dos resultados.
- `pdb_file`: Caminho absoluto para o arquivo PDB.
- `traj_file`: Caminho absoluto para o arquivo de trajetória.
- `output_file`: Caminho absoluto para o arquivo de saída (.dat).
- `num_frames`: Número de frames a serem processados (padrão: 150000).
"""
function calc_median(pdb_file::AbstractString, traj_file::AbstractString, output_file::AbstractString; num_frames::Union{Nothing, Int}=nothing, seg_select::AbstractString)
    # Carregar a simulação e os dados do PDB
    simulation, poli_data = load_simulation(pdb_file, traj_file, seg_select)
    
    println("Saída: $output_file")
    
    # Preparar os índices dos monômeros
    monomer_indices = prepare_monomer_indices(poli_data)

    # Dentro da função calc_median, extraia o nome do segmento:
    seg_words = split(seg_select)
    if length(seg_words) < 2
        error("O argumento seg_select deve conter pelo menos duas palavras, por exemplo: 'segname POLI'")
    end
    seg_name = seg_words[2]

    
    # Obter os índices dos átomos do polímero na simulação
    polymer_indices = findall(atom -> atom.segname == "$seg_name", atoms(simulation))

     # Se num_frames não for fornecido, usa todos os frames disponíveis
     if isnothing(num_frames)
        num_frames = length(simulation)
    end
    
    # Processar os frames para calcular as distâncias
    distances_all_frames = analyze_distances(simulation, polymer_indices, monomer_indices, num_frames)
    
    # Calcular as medianas e plotar os resultados
    medians = analyze_and_plot(distances_all_frames)
    
    # Salvar os dados no arquivo de saída
    save_to_dat_file(medians, output_file)
    
    return medians
end