
export get_contribution_matrix, plot_diff_contributions

using Plots, LaTeXStrings, ComplexMixtures, PDBTools, EasyFit
const CM, PDB = ComplexMixtures, PDBTools
import EasyFit: movavg

# Grupos químicos dos monômeros
const GROUPS = [
  (["CF","HF1","HF2","HF3"], L"\\textrm{CH_3-t}"),
  (["C1","H11","H13","C2","H21"], L"\\textrm{CHCH_2}"),
  (["O1","C3"], "CO"),
  (["N1","HN1","HN2"], L"\\textrm{NH_2}"),
  (["CZ","HZ1","HZ2","HZ3"], L"\\textrm{CH_3}")
]

function configure_plots()
    default(
        fontfamily  = "Computer Modern",
        linewidth   = 2.5,
        framestyle  = :box,
        label       = nothing,
        grid        = false,
    )
    scalefontsizes(1.3)
end

function extract_group_contributions(pdb_file::String, json_file::String, groups::Vector=GROUPS)
    system = readPDB(pdb_file)
    acr = select(system, "resname AAMD")
    results = CM.load(json_file)

    group_contribs = Vector{Float64}[]
    labels = LaTeXString[]

    for (imer, mer) in enumerate(eachresidue(acr))
        for (group_atoms, group_label) in groups
            if occursin("CH_3-t", String(group_label)) && imer ≠ 1 && imer ≠ length(eachresidue(acr))
                continue
            end

            mer_group_atoms = filter(at -> name(at) in group_atoms, mer)
            if isempty(mer_group_atoms)
                continue
            end

            contrib = contributions(results, SoluteGroup(mer_group_atoms))
            contrib = movavg(contrib; n=10).x

            push!(group_contribs, contrib)
            push!(labels, group_label)
        end
    end

    return results.d, hcat(group_contribs...), labels
end

function get_contribution_matrix(pdb_file::String, json_file::String)
    return extract_group_contributions(pdb_file, json_file)
end

function annotate_residues!(plt, res_b, res_e, x_pos)
    y_pos = 3.65
    for res in res_b:res_e
        txt = L"$" * string(res) * L"$"
        annotate!(plt, x_pos + 3 * ((res - res_b) % 25), y_pos, text(txt, 30))
    end
end

function plot_diff_contributions(d::Vector, diff_contribs::Matrix, labels::Vector{LaTeXString}; filename="diferenca_grupos.png")
    configure_plots()

    idmin = findfirst(>(1.9), d)
    idmax = findfirst(>(3.5), d)
    plotvec = []

    for i in 1:4
        b, e, res_b, res_e, x_respos, title = if i == 1
            (1, 76, 1, 25, 4, L"\\mathrm{Residue}")
        elseif i == 2
            (77, 151, 26, 50, 3, "")
        elseif i == 3
            (152, 226, 51, 75, 3, "")
        else
            (227, 302, 76, 100, 3, "")
        end

        b = min(b, length(labels))
        e = min(e, length(labels))
        new_labels = labels[b:e]

        p = contourf(
            1:length(new_labels),
            d[idmin:idmax],
            diff_contribs[idmin:idmax, b:e],
            title           = title,
            xlabel          = (i == 4 ? L"\\mathrm{Group}" : ""),
            ylabel          = L"r/\\AA",
            xticks          = (1:length(new_labels), new_labels),
            titlefontsize   = 40,
            xguidefontsize  = 40,
            yguidefontsize  = 30,
            xrotation       = 60,
            xtickfont       = font(16, "Computer Modern"),
            ytickfont       = font(25, "Computer Modern"),
            color           = cgrad(:balance),
            linewidth       = 1,
            linecolor       = :black,
            colorbar        = :none,
            levels          = 10,
            size            = (1800, 600),
            dpi             = 500,
        )
        annotate_residues!(p, res_b, res_e, x_respos)
        push!(plotvec, p)
    end

    final_plot = plot(plotvec..., layout = (4, 1), size = (2800, 2800), margin = 1.75Plots.Measures.cm)
    savefig(final_plot, filename)
end



