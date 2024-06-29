module TargetAquisitionReport

using Statistics

import ArgParse
import CSV
import DataFrames
import ProgressMeter: @showprogress
import UniMZ
import UniMZUtil: Proteomics, TMS, pFind

include("../util.jl")

prepare(args) = begin
    path_ms = args["ms"]
    @info "file path of selected data:"
    println("\t$(path_ms)")
    paths_ms_old = reduce(vcat, UniMZ.match_path.(args["ms_old"], ".umz")) |> unique |> sort
    @info "file paths of selected original data:"
    foreach(x -> println("$(x[1]):\t$(x[2])"), enumerate(paths_ms_old))
    path_psm = args["psm"]
    out = mkpath(args["out"])
    path_ft = args["ft"]
    fmt = args["fmt"] |> Symbol
    ε = parse(Float64, args["error"]) * 1.0e-6
    fdr = parse(Float64, args["fdr"]) / 100
    decoy = args["decoy"]::Bool
    τ_ms_sim = parse(Float64, args["ms_sim_thres"])
    tab_ele, tab_aa, tab_mod = pFind.read_mass_table(args["cfg"])
    return (; path_ms, paths_ms_old, path_psm, out, path_ft, fmt, ε, fdr, decoy, τ_ms_sim, tab_ele, tab_aa, tab_mod)
end

process(path; path_ms, paths_ms_old, path_psm, out, path_ft, fmt, ε, fdr, decoy, τ_ms_sim, tab_ele, tab_aa, tab_mod) = begin
    M = UniMZ.read_ms(path_ms)
    df_m2 = map(m -> (; m.id, mz=m.activation_center, rt=m.retention_time, m.peaks), M.MS2) |> DataFrames.DataFrame
    M2I = map(x -> x[2] => x[1], enumerate(df_m2.id)) |> Dict

    M_old = map(p -> splitext(basename(p))[1] => UniMZ.dict_by_id(UniMZ.read_ms(p).MS2), paths_ms_old) |> Dict

    df = pFind.read_psm(path_psm)
    df.engine .= :pFind
    filter!(r -> r.fdr .≤ fdr, df)
    !decoy && filter!(r -> r.td == :T, df)

    ns = [
        "Order", "Peptide", "Peptide_Type", "mh_calc", "Modifications", "Evalue", "Precursor_Mass_Error(Da)",
        "Proteins", "prot_type", "FileID", "LabelID", "Alpha_Matched", "Beta_Matched", "Alpha_Evalue", "Beta_Evalue",
        "Alpha_Seq_Coverage", "Beta_Seq_Coverage",
    ]
    DataFrames.select!(df, DataFrames.Not(filter(x -> x ∈ names(df), ns)))
    ns = ["engine", "mh", "mz", "z", "pep", "mod", "prot", "title", "file", "scan", "idx_pre"]
    DataFrames.select!(df, ns, DataFrames.Not(ns))

    ion_syms = ["b", "y"]
    ion_types = map(i -> getfield(UniMZ, Symbol("ion_$(i)")), ion_syms)
    M_ = [splitext(basename(path_ms))[1] => UniMZ.dict_by_id(M.MS2)] |> Dict
    Proteomics.calc_cov!(df, M_, ε, ion_syms, ion_types, tab_ele, tab_aa, tab_mod)

    df.credible = map(r -> r.cov_ion_y ≥ 0.6 && r.cov_ion_b ≥ 0.4, eachrow(df))

    ns = [
        "Scan_No", "Sequence", "mh_calc", "Mass_Shift(Exp.-Calc.)", "score_raw", "Modification",
        "Specificity", "Positions", "Label", "Miss.Clv.Sites", "Avg.Frag.Mass.Shift", "Others", "mz_calc"
    ]
    DataFrames.select!(df, DataFrames.Not(filter(x -> x ∈ names(df), ns)))

    df.id = Vector(1:size(df, 1))
    DataFrames.select!(df, :id, DataFrames.Not([:id]))
    df.rt = [df_m2[M2I[r.scan], :rt] for r in eachrow(df)]

    df_m2.psm = [df[df.scan .== r.id, :id] for r in eachrow(df_m2)]

    !isempty(path_ft) && @info "Feature loading from "* path_ft
    df_ft = (isempty(path_ft) ? [] : CSV.File(path_ft)) |> DataFrames.DataFrame
    df_ft.id = Vector(1:size(df_ft, 1))
    DataFrames.select!(df_ft, :id, DataFrames.Not([:id]))

    @info "Target loading from " * path
    df_tg = CSV.File(path) |> DataFrames.DataFrame
    df_tg.id = Vector(1:size(df_tg, 1))
    TMS.parse_target_list!(df_tg, fmt)
    DataFrames.select!(df_tg, [:id, :mz, :z, :start, :stop], DataFrames.Not([:id, :mz, :z, :start, :stop]))
    "mod_a" ∈ names(df_tg) && (df_tg.mod_a = parse.(Array{UniMZ.Mod}, unify_mods_str.(df_tg.mod_a)))
    "mod_b" ∈ names(df_tg) && (df_tg.mod_b = parse.(Array{UniMZ.Mod}, unify_mods_str.(df_tg.mod_b)))

    @info "Feature mapping"
    tmp = sort!([(x.mz::Float64, x.id::Int) for x in eachrow(df_ft)])
    mzs = map(x -> x[1], tmp)
    ids = map(x -> x[2], tmp)
    df_tg.ft_ = [sort(filter(x -> df_ft[x, :z] == r.z, ids[UniMZ.argquery_ε(mzs, r.mz, ε)])) for r in eachrow(df_tg)]
    df_tg.n_ft = length.(df_tg.ft_)

    Ks = ["", "_allsim", "_all"]

    calc_sim(dda, tda) = map(p -> !isempty(UniMZ.argquery_ε(tda.peaks, p.mz, ε)), M_old[dda.file][dda.scan].peaks) |> mean

    @info "MS2 mapping"
    tmp = sort!([(x.mz::Float64, x.id::Int) for x in eachrow(df_m2)])
    mzs = map(x -> x[1], tmp)
    ids = map(x -> x[2], tmp)
    df_tg.m2_all_ = [map(x -> M2I[x], sort(ids[UniMZ.argquery_ε(mzs, r.mz, ε)])) for r in eachrow(df_tg)]
    df_tg.m2_allsim_ = [filter(i -> calc_sim(r, df_m2[i, :]) ≥ τ_ms_sim, r.m2_all_) for r in eachrow(df_tg)]
    df_tg.m2_ = [filter(i -> r.start ≤ df_m2.rt[i] ≤ r.stop, r.m2_all_) for r in eachrow(df_tg)]
    for K in Ks
        df_tg[!, "n_m2$(K)"] = length.(df_tg[!, "m2$(K)_"])
    end
    for K in Ks
        df_tg[!, "m2$(K)_id_"] = [map(x -> df_m2.id[x], r["m2$(K)_"]) for r in eachrow(df_tg)]
    end

    @info "MS2 similarity calculating"
    for K in Ks
        df_tg[!, "m2$(K)_sim_"] = @showprogress map(eachrow(df_tg)) do r
            map(s -> "$(s.id):$(round(calc_sim(r, s); digits=2))", eachrow(df_m2[r["m2$(K)_"], :])) |> xs -> join(xs, ";")
        end
    end

    @info "PSM mapping"
    for K in Ks
        df_tg[!, "psm$(K)_"] = [vcat(df_m2[r["m2$(K)_"], :psm]...) for r in eachrow(df_tg)]
        df_tg[!, "n_psm$(K)"] = length.(df_tg[!, "psm$(K)_"])
    end

    ns = filter(n -> !endswith(n, '_'), names(df_tg))
    DataFrames.select!(df_tg, ns, DataFrames.Not(ns))

    for K in Ks
        df_tg[!, "iden$(K)"] = map(df_tg[!, "psm$(K)_"]) do psms
            map(psmstr, eachrow(df[psms, :])) |> xs -> join(xs, ";")
        end
        df_tg[!, "have_iden$(K)"] = .!isempty.(df_tg[!, "iden$(K)"])
        df_tg[!, "iden$(K)_credible"] = map(df_tg[!, "psm$(K)_"]) do psms
            map(psmstr, eachrow(df[filter(i -> df.credible[i], psms), :])) |> xs -> join(xs, ";")
        end
        df_tg[!, "have_iden$(K)_credible"] = .!isempty.(df_tg[!, "iden$(K)_credible"])
    end

    for K in Ks
        df_tg[!, "same_iden$(K)"] = map(eachrow(df_tg)) do r
            filter(eachrow(df[r["psm$(K)_"], :])) do s
                (s.engine != :pFind) && return false
                return is_same_pepmod(r, s)
            end .|> psmstr |> xs -> join(xs, ";")
        end
        df_tg[!, "have_same_iden$(K)"] = .!isempty.(df_tg[!, "same_iden$(K)"])

        df_tg[!, "same_iden$(K)_credible"] = map(eachrow(df_tg)) do r
            filter(eachrow(df[r["psm$(K)_"], :])) do s
                (s.engine != :pFind) && return false
                (!s.credible) && return false
                return is_same_pepmod(r, s)
            end .|> psmstr |> xs -> join(xs, ";")
        end
        df_tg[!, "have_same_iden$(K)_credible"] = .!isempty.(df_tg[!, "same_iden$(K)_credible"])

        vs = map(eachrow(df_tg)) do r
            b = nothing
            for s in eachrow(df[r["psm$(K)_"], :])
                (s.engine != :pFind) && continue
                isnothing(b) && (b = s)
                (s.cov ≤ b.cov) && continue
                b = s
            end
            return isnothing(b) ? ("", false) : (psmstr(b), is_same_pepmod(r, b))
        end
        df_tg[!, "best_iden$(K)"] = first.(vs)
        df_tg[!, "best_iden$(K)_is_same_iden"] = last.(vs)
    end

    UniMZ.safe_save(p -> CSV.write(p, df_tg), joinpath(out, "$(basename(splitext(path_ms)[1])).target.TargetAquisitionReport.csv"))
    UniMZ.safe_save(p -> CSV.write(p, df), joinpath(out, "$(basename(splitext(path_ms)[1])).psm.TargetAquisitionReport.csv"))
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="TargetAquisitionReport")
    ArgParse.@add_arg_table! settings begin
        "target"
            help = "target list"
            required = true
        "--ms"
            help = ".umz or .ms1/2 file; .ms2/1 file should be in the same directory for .ms1/2"
            required = true
        "--ms_old"
            help = "origianl .umz or .ms1/2 files; .ms2/1 files should be in the same directory for .ms1/2"
            nargs = '+'
            required = true
        "--psm"
            help = "pFind PSM path"
            required = true
        "--out"
            help = "output directory"
            default = "./out/"
            metavar = "./out/"
        "--ft"
            help = "feature list"
            default = ""
        "--fmt", "-f"
            help = "target list format: auto, TW, TmQE, TmFu"
            metavar = "auto|TW|TmQE|TmFu"
            default = "auto"
        "--error"
            help = "m/z error"
            metavar = "ppm"
            default = "20.0"
        "--fdr"
            help = "FDR threshold (%)"
            default = "Inf"
        "--decoy"
            help = "preserve decoy identifications"
            action = :store_true
        "--ms_sim_thres"
            help = "threshold of MS similarity"
            default = "0.5"
        "--cfg"
            help = "pFind config directory"
            default = ""
    end
    args = ArgParse.parse_args(settings)
    process(args["target"]; prepare(args)...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

julia_main()::Cint = begin
    main()
    return 0
end

end
