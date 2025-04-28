using Pkg ; Pkg.activate("/home/m3g/Documentos/1_analysis/1_analises/p100_900ns_analysis/p100_900ns_analysis/environments")
using PDBTools, ComplexMixtures, Plots, LaTeXStrings, EasyFit

function main()
    # Leitura do sistema e carregamento dos resultados
    system = read_pdb("./1.1_100_s1_05/p100_box.pdb")
    acr = select(system, "resname AAMD")
    #results = load("./json_median/water_extended_05_average.json")
    results = load("./json_median/glu_full_05_average.json")

    # Cálculo das contribuições por resíduos
    cn1 = ResidueContributions(results, select(acr, at -> name(at) in ["O1", "C3"]); dmin=1.5, dmax=6)
    cn2 = ResidueContributions(results, select(acr, at -> name(at) in ["N1", "HN1", "HN2"]); dmin=1.5, dmax=6)
    cn3 = ResidueContributions(results, select(acr, at -> name(at) in ["C1", "H11", "H13", "C2", "H21"]); dmin=1.5, dmax=6)

    # Ajustes de eixos e escala
    idmin, idmax = 20, 40
    dmin, dmax = 1.5, 6
    clims=(0, 0.02)
    font_size_label = 19
    font_size_stick = 20
    cor_plots = :tempo
    plt_type = contourf

    # Ajustar as coordenadas da caixa para centralizá-la no gráfico
    x_min, x_max = idmin + 0.50, idmin + 2.5  # Coordenadas horizontais ajustadas
    y_min, y_max = dmax - 1.0, dmax - 0.5    # Coordenadas verticais ajustadas

    # Plot de cada grupo
    p1 = plt_type(cn1[idmin:idmax]; 
            clims=clims, 
            color=cor_plots, 
            ylabel=L"r/\AA", 
            xlabel="",
            xguidefontsize=font_size_stick,
            yguidefontsize=font_size_stick,
            xtickfont=font(font_size_stick, "Computer Modern"),
            ytickfont=font(font_size_stick, "Computer Modern"),
            linewidth=1,
            xticks=nothing,
            xlim=(idmin, idmax),          # Definir limites do eixo x
            ylim=(dmin, dmax))            # Definir limites do eixo y

    # Adicionar o quadrado preenchido no gráfico menor
    plot!([x_min, x_max, x_max, x_min, x_min],  # Coordenadas x
            [y_min, y_min, y_max, y_max, y_min],  # Coordenadas y
            seriestype = :shape,                  # Tipo de série: forma preenchida
            color = :white,                       # Cor do preenchimento
            alpha = 0.8,                          # Transparência ajustada
            label = nothing)                      # Sem legenda

    # Adicionar anotação dentro do quadrado no gráfico menor
    annotate!((x_min + (x_max - x_min) / 2,     # Coordenada x centralizada na caixa
            y_min + (y_max - y_min) / 2,     # Coordenada y centralizada na caixa
            text(L"\textrm{CO}", "Computer Modern", font_size_label, :black)))  # Anotação no centro

    p2 = plt_type(cn2[idmin:idmax];  # contourf
            clims=clims, 
            color=cor_plots, 
            ylabel=L"r/\AA", 
            xlabel="",
            xguidefontsize=font_size_stick,
            yguidefontsize=font_size_stick,
            xtickfont=font(font_size_stick, "Computer Modern"),
            ytickfont=font(font_size_stick, "Computer Modern"),
            linewidth=1,
            xticks=nothing,
            xlim=(idmin, idmax),          # Definir limites do eixo x
            ylim=(dmin, dmax))            # Definir limites do eixo y)

    # Adicionar o quadrado preenchido no gráfico menor
    plot!([x_min, x_max, x_max, x_min, x_min],  # Coordenadas x
            [y_min, y_min, y_max, y_max, y_min],  # Coordenadas y
            seriestype = :shape,                  # Tipo de série: forma preenchida
            color = :white,                       # Cor do preenchimento
            alpha = 0.8,                          # Transparência ajustada
            label = nothing)                      # Sem legenda

    # Adicionar anotação dentro do quadrado no gráfico menor
    annotate!((x_min + (x_max - x_min) / 2,     # Coordenada x centralizada na caixa
            y_min + (y_max - y_min) / 2,     # Coordenada y centralizada na caixa
            text(L"\textrm{NH_2}", "Computer Modern", font_size_label, :black)))  # Anotação no centro

    p3 = plt_type(cn3[idmin:idmax]; 
            clims=clims, 
            color=cor_plots, 
            ylabel=L"r/\AA", 
            xguidefontsize=font_size_stick,
            yguidefontsize=font_size_stick,
            xtickfont=font(font_size_stick, "Computer Modern"),
            ytickfont=font(font_size_stick, "Computer Modern"),
            linewidth=1,
            xlim=(idmin, idmax),          # Definir limites do eixo x
            ylim=(dmin, dmax),            # Definir limites do eixo y)
            xticks=(idmin:idmax,idmin:idmax) #20:40,20:40)
            )

    # Adicionar o quadrado preenchido no gráfico menor
    plot!([x_min, x_max, x_max, x_min, x_min],  # Coordenadas x
            [y_min, y_min, y_max, y_max, y_min],  # Coordenadas y
            seriestype = :shape,                  # Tipo de série: forma preenchida
            color = :white,                       # Cor do preenchimento
            alpha = 0.8,                          # Transparência ajustada
            label = nothing)                      # Sem legenda

    # Adicionar anotação dentro do quadrado no gráfico menor
    annotate!((x_min + (x_max - x_min) / 2,     # Coordenada x centralizada na caixa
            y_min + (y_max - y_min) / 2,     # Coordenada y centralizada na caixa
            text(L"\textrm{CHCH_2}", "Computer Modern", font_size_label, :black)))  # Anotação no centro

    # Layout final
    p_final = plot([p1, p2, p3]...; 
            layout=(3, 1), 
            size=(1200, 1200), 
            leftmargin=5Plots.Measures.mm,    # Definindo as margem esquerda do plot final
            rightmargin=5Plots.Measures.mm,   # Definição margem direita
            topmargin=0Plots.Measures.mm)

    # Salvar o gráfico
    savefig("residue_contributions_glucose.png")
    println("Gráfico salvo como 'residue_contributions.png'")

    return p_final
end

main()