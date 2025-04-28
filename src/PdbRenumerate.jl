export pdbRen

using Printf
function pdbRen(pdb_path::String, atomos_para_remover::Vector{Int}, output_path::String)
    linhas = readlines(pdb_path)
    novas_linhas = String[]
    novo_serial = 1

    for linha in linhas
        if startswith(linha, "ATOM") || startswith(linha, "HETATM")
            serial_str = strip(linha[7:11])
            try
                serial = parse(Int, serial_str)
                if serial ∈ atomos_para_remover
                    continue  # pula o átomo a ser removido
                end

                # Substitui o número do átomo com alinhamento
                linha_nova = string(linha[1:6], @sprintf("%5d", novo_serial), linha[12:end])
                push!(novas_linhas, linha_nova)
                novo_serial += 1
            catch
                push!(novas_linhas, linha)  # se falhar, mantém a linha
            end
        else
            push!(novas_linhas, linha)  # outras linhas são preservadas
        end
    end

    open(output_path, "w") do io
        for linha in novas_linhas
            println(io, linha)
        end
    end

    println("✅ Arquivo salvo como: $output_path")
end

# Exemplo de uso:
# pdbRen("entrada.pdb", [1, 1008], "saida_filtrada.pdb")