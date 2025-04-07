#using CSV, DataFrames, Statistics, Dates

export frames_counter

"""
    main(folders::Vector{String}; file_folder::String, file_analise, multiplier::Real)

Função principal que itera sobre as pastas especificadas em `folders`. Em cada pasta, a função
procura na subpasta `file_folder` o arquivo de análise e realiza o processamento dos dados.
O parâmetro `file_analise` pode ser:
  - Uma _String_ com o nome do arquivo.
  - Uma _Tuple_ (prefix, extensão) para buscar um arquivo cujo nome inicie com `prefix` e termine com `extensão`.
O `multiplier` é usado para calcular o limiar (threshold) a partir do desvio padrão.
"""
function frames_counter(folders::Vector{String}; file_folder::String, file_analise, multiplier::Real)
    base_path = pwd()

    println("Iniciando a análise de frames...")

    for folder in folders
        # Caminho para a pasta principal
        path_main = joinpath(base_path, folder)
        if !isdir(path_main)
            println("❌ Diretório não encontrado: $folder. Pulando para o próximo.")
            continue
        end

        # Caminho para a subpasta onde o arquivo de análise está localizado
        sub_path = joinpath(path_main, file_folder)
        if !isdir(sub_path)
            println("❌ Diretório '$file_folder' não encontrado em $folder. Pulando para o próximo.")
            continue
        end

        # Determinando o nome do arquivo de análise com base no parâmetro recebido
        chosen_file = nothing
        if isa(file_analise, Tuple{String,String})
            prefix, ext = file_analise
            # Busca por um arquivo que inicia com `prefix` e termina com `ext`
            candidates = filter(fname -> startswith(fname, prefix) && endswith(fname, ext), readdir(sub_path))
            if isempty(candidates)
                println("⚠️ Nenhum arquivo encontrado em $sub_path que inicie com \"$prefix\" e termine com \"$ext\". Pulando para o próximo.")
                continue
            end
            chosen_file = candidates[1]  # Utiliza o primeiro arquivo encontrado
        elseif isa(file_analise, String)
            chosen_file = file_analise
        else
            println("⚠️ Parâmetro file_analise inválido em $folder. Pulando para o próximo.")
            continue
        end

        file_path = joinpath(sub_path, chosen_file)
        if !isfile(file_path)
            println("⚠️ Arquivo '$chosen_file' não encontrado em $sub_path. Pulando para o próximo.")
            continue
        end

        # Lendo os dados do arquivo
        s1 = CSV.read(file_path, DataFrame; header=false)

        # Renomeando as colunas para facilitar a manipulação
        rename!(s1, Dict(:Column1 => :Step, :Column2 => :Radius_of_Gyration))

        # Criando a coluna "Frame"
        s1.Frame = collect(0:size(s1, 1)-1)

        # Calculando estatísticas
        mean_val = mean(s1.Radius_of_Gyration)
        std_val  = std(s1.Radius_of_Gyration)
        threshold = multiplier * std_val

        # Selecionando os frames com base no threshold
        extended_frames = collect(s1.Frame[s1.Radius_of_Gyration .>= (mean_val + threshold)])
        compact_frames  = collect(s1.Frame[s1.Radius_of_Gyration .<= (mean_val - threshold)])

        # Criando o arquivo de log na mesma subpasta
        log_file = "frames_selection_info.log"
        log_path = joinpath(sub_path, log_file)
        open(log_path, "w") do file
            println(file, """📄 Log de Seleção de Frames - $folder/$file_folder

🕒 Data e Hora da Execução: $(now())

📂 Arquivo de análise: $chosen_file

📊 Estatísticas:
   - Média do raio de giro: $mean_val
   - Desvio padrão: $std_val
   - Limiar (Threshold): $threshold

🎯 Frames Selecionados:
   - Total de frames estendidos: $(length(extended_frames))
   - Total de frames compactos : $(length(compact_frames))

   - Frames Estendidos: $(join(extended_frames, ", "))
   - Frames Compactos : $(join(compact_frames, ", "))

✅ Execução concluída.
""")
        end

        # Exibindo um resumo dos resultados no terminal
        println("\n📂 Pasta: $folder")
        println("   Subpasta: $file_folder")
        println("   Arquivo: $chosen_file")
        println("📊 Estatísticas:")
        println("   - Média do raio de giro: $mean_val")
        println("   - Desvio padrão: $std_val")
        println("   - Limiar (Threshold): $threshold")
        println("   - Total de frames estendidos: $(length(extended_frames))")
        println("   - Total de frames compactos : $(length(compact_frames))")
    end
    println("✅ Análise de frames concluída.")

end

# Exemplo de uso:
# folders = ["0", "20", "40", "60", "80", "100"]
# Aqui passamos uma tupla para file_analise: buscar um arquivo que comece com "median" e termine com ".dat"
# frames_counter(folders; file_folder="analises", file_analise=("radius", ".dat"), multiplier=1.55)
