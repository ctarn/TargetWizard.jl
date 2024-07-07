"""
Target Acquisition Report for Crosslink
"""
module TargetAcquisitionXLReport

"""
- target list (multiple formats supported), e.g., `TargetWizard.all.TW.target.csv`.
- traditional and targeted mass spectrometry data, e.g., `DDA.raw` and `TMS.raw`.
  - The raw data should be converted into an open-source format such as MS1/MS2. [ThermoRawRead](http://thermorawread.ctarn.io) is recommended.
- (filtered) crosslink identification results of targeted mass spectrometry data, e.g., `TMS.plink.csv`.
- optional: (filtered) linear peptide identification results of targeted mass spectrometry data, e.g., `TMS_fdr.pfind.csv`.
- optional: precursor list detected by [`PepPre`](http://peppre.ctarn.io)
- optional: candidate crosslink list
"""
require = true

"""
Once finished, TargetWizard will save two reports to `Output Directroy`.
- `csv` report of all targets, e.g., `TMS.target.TargetAcquisitionXLReport.csv`.
- `csv` report of all crosslink PSMs, e.g., `TMS.crosslink.TargetAcquisitionXLReport.csv`.
- `csv` report of all linear peptide PSMs, e.g., `TMS.linear.TargetAcquisitionXLReport.csv`.
- `csv` report of all monolink PSMs, e.g., `TMS.monolink.TargetAcquisitionXLReport.csv`.
- `csv` report of all looplink PSMs, e.g., `TMS.looplink.TargetAcquisitionXLReport.csv`.
"""
output = true

"""
# Max. MS1 Mass Error
mass error used to match targets, PSMs, and MS scans.

# FDR Threshold
used to filter PSM list.

# MS Sim. Thres.
used to match traditional and targeted MS scans.

![Target Acquisition Report for Crosslink](../assets/report/TargetAcquisitionXLReport.png)
"""
usage = true

"""
"""
example = true

# TODO remove pfind support

using Statistics

import ArgParse
import CSV
import DataFrames
import ProgressMeter: @showprogress
import UniMZ
import UniMZUtil: Proteomics, Crosslink, TMS, pFind, pLink

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
    path_xl = args["xl"]
    path_ft = args["ft"]
    path_psm_pf = args["psm_pf"]
    fmt = args["fmt"] |> Symbol
    linker = Symbol(args["linker"])
    ε = parse(Float64, args["error"]) * 1.0e-6
    fdr = parse(Float64, args["fdr"]) / 100
    decoy = args["decoy"]::Bool
    τ_ms_sim = parse(Float64, args["ms_sim_thres"])
    tab_ele_pl, tab_aa_pl, tab_mod_pl, tab_xl_pl = pLink.read_mass_table(args["cfg"])
    tab_ele_pf, tab_aa_pf, tab_mod_pf = pFind.read_mass_table(args["cfg_pf"])
    return (; path_ms, paths_ms_old, path_psm, out, path_xl, path_ft, path_psm_pf, fmt, linker, ε, fdr, decoy, τ_ms_sim, tab_ele_pl, tab_aa_pl, tab_mod_pl, tab_xl_pl, tab_ele_pf, tab_aa_pf, tab_mod_pf)
end

process(path; path_ms, paths_ms_old, path_psm, out, path_xl, path_ft, path_psm_pf, fmt, linker, ε, fdr, decoy, τ_ms_sim, tab_ele_pl, tab_aa_pl, tab_mod_pl, tab_xl_pl, tab_ele_pf, tab_aa_pf, tab_mod_pf) = begin
    M = UniMZ.read_ms(path_ms)
    df_m2 = map(m -> (; m.id, mz=m.activation_center, rt=m.retention_time, m.peaks), M.MS2) |> DataFrames.DataFrame
    M2I = map(x -> x[2] => x[1], enumerate(df_m2.id)) |> Dict

    M_old = map(p -> splitext(basename(p))[1] => UniMZ.dict_by_id(UniMZ.read_ms(p).MS2), paths_ms_old) |> Dict

    dfs = pLink.read_psm_full(path_psm; linker)
    df_psm = dfs.xl
    df_linear = dfs.linear
    df_mono = dfs.mono
    df_loop = dfs.loop

    for df in [df_psm, df_linear, df_mono, df_loop]
        df.engine .= :pLink
        filter!(r -> r.fdr .≤ fdr, df)
        !decoy && filter!(r -> r.td == :TT || r.td == :T, df)
    end

    ns = [
        "Order", "Peptide", "Peptide_Type", "mh_calc", "Modifications", "Evalue", "Precursor_Mass_Error(Da)",
        "Proteins", "prot_type", "FileID", "LabelID", "Alpha_Matched", "Beta_Matched", "Alpha_Evalue", "Beta_Evalue",
        "Alpha_Seq_Coverage", "Beta_Seq_Coverage",
    ]
    DataFrames.select!(df_psm, DataFrames.Not(filter(x -> x ∈ names(df_psm), ns)))
    ns = [
        "engine", "mh", "mz", "z", "pep_a", "pep_b", "mod_a", "mod_b", "site_a", "site_b",
        "prot_a", "prot_b", "error", "title", "file", "scan", "idx_pre",
    ]
    DataFrames.select!(df_psm, ns, DataFrames.Not(ns))

    ion_syms = ["b", "y"]
    ion_types = map(i -> getfield(UniMZ, Symbol("ion_$(i)")), ion_syms)
    M_ = [splitext(basename(path_ms))[1] => UniMZ.dict_by_id(M.MS2)] |> Dict
    Crosslink.calc_cov_crosslink!(df_psm, M_, ε, ion_syms, ion_types, tab_ele_pl, tab_aa_pl, tab_mod_pl, tab_xl_pl)
    Crosslink.calc_cov_linear!(df_linear, M_, ε, ion_syms, ion_types, tab_ele_pl, tab_aa_pl, tab_mod_pl)
    Crosslink.calc_cov_monolink!(df_mono, M_, ε, ion_syms, ion_types, tab_ele_pl, tab_aa_pl, tab_mod_pl, tab_xl_pl)
    Crosslink.calc_cov_looplink!(df_loop, M_, ε, ion_syms, ion_types, tab_ele_pl, tab_aa_pl, tab_mod_pl, tab_xl_pl)

    df_psm.cov_min = min.(df_psm.cov_a, df_psm.cov_b)
    df_psm.credible = map(eachrow(df_psm)) do r
        r.cov_a_ion_y ≥ 0.6 && r.cov_b_ion_y ≥ 0.6 && r.cov_a_ion_b ≥ 0.4 && r.cov_b_ion_b ≥ 0.4
    end
    df_linear.credible = map(r -> r.cov_ion_y ≥ 0.6 && r.cov_ion_b ≥ 0.4, eachrow(df_linear))
    df_mono.credible = map(r -> r.cov_ion_y ≥ 0.6 && r.cov_ion_b ≥ 0.4, eachrow(df_mono))
    df_loop.credible = map(r -> r.cov_ion_y ≥ 0.6 && r.cov_ion_b ≥ 0.4, eachrow(df_loop))

    if !isempty(path_psm_pf)
        df_psm_pf = pFind.read_psm(path_psm_pf)
        df_psm_pf.engine .= :pFind
        ns = [
            "Scan_No", "Sequence", "mh_calc", "Mass_Shift(Exp.-Calc.)", "score_raw", "Modification",
            "Specificity", "Positions", "Label", "Miss.Clv.Sites", "Avg.Frag.Mass.Shift", "Others", "mz_calc"
        ]
        DataFrames.select!(df_psm_pf, DataFrames.Not(filter(x -> x ∈ names(df_psm_pf), ns)))
        Proteomics.calc_cov!(df_psm_pf, M_, ε, ion_syms, ion_types, tab_ele_pf, tab_aa_pf, tab_mod_pf)
        DataFrames.rename!(df_psm_pf, :pep => :pep_a, :mod => :mod_a, :prot => :prot_a)
        df_psm_pf.credible = map(r -> r.cov_ion_y ≥ 0.6 && r.cov_ion_b ≥ 0.4, eachrow(df_psm_pf))
        df_psm = vcat(df_psm, df_psm_pf; cols=:union)
    end

    for df in [df_psm, df_linear, df_mono, df_loop]
        df.id = Vector(1:size(df, 1))
        DataFrames.select!(df, :id, DataFrames.Not([:id]))
        df.rt = [df_m2[M2I[r.scan], :rt] for r in eachrow(df)]
    end

    df_m2.psm = [df_psm[df_psm.scan .== r.id, :id] for r in eachrow(df_m2)]
    df_m2.psm_linear = [df_linear[df_linear.scan .== r.id, :id] for r in eachrow(df_m2)]
    df_m2.psm_mono = [df_mono[df_mono.scan .== r.id, :id] for r in eachrow(df_m2)]
    df_m2.psm_loop = [df_loop[df_loop.scan .== r.id, :id] for r in eachrow(df_m2)]

    !isempty(path_xl) && @info "XL Candidtes loading from " * path_xl
    df_xl = (isempty(path_xl) ? [] : CSV.File(path_xl)) |> DataFrames.DataFrame
    df_xl.id = Vector(1:size(df_xl, 1))
    DataFrames.select!(df_xl, :id, DataFrames.Not([:id]))

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

    @info "XL Candidtes mapping"
    tmp = sort!([(x.mz::Float64, x.id::Int) for x in eachrow(df_xl)])
    mzs = map(x -> x[1], tmp)
    ids = map(x -> x[2], tmp)
    df_tg.xl_ = [sort(filter(x -> df_xl[x, :z] == r.z, ids[UniMZ.argquery_ε(mzs, r.mz, ε)])) for r in eachrow(df_tg)]

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
    filter_plink(x) = filter(i -> df_psm.engine[i] == :pLink, x)
    filter_pfind(x) = filter(i -> df_psm.engine[i] == :pFind, x)
    for K in Ks
        df_tg[!, "psm$(K)_"] = [vcat(df_m2[r["m2$(K)_"], :psm]...) for r in eachrow(df_tg)]
        df_tg[!, "n_psm$(K)"] = length.(df_tg[!, "psm$(K)_"])

        df_tg[!, "psm_linear$(K)_"] = [vcat(df_m2[r["m2$(K)_"], :psm_linear]...) for r in eachrow(df_tg)]
        df_tg[!, "psm_mono$(K)_"] = [vcat(df_m2[r["m2$(K)_"], :psm_mono]...) for r in eachrow(df_tg)]
        df_tg[!, "psm_loop$(K)_"] = [vcat(df_m2[r["m2$(K)_"], :psm_loop]...) for r in eachrow(df_tg)]

        df_tg[!, "n_psm_plink$(K)"] = length.(map(filter_plink, df_tg[!, "psm$(K)_"]))
        df_tg[!, "n_psm_pfind$(K)"] = length.(map(filter_pfind, df_tg[!, "psm$(K)_"]))
    end

    ns = filter(n -> !endswith(n, '_'), names(df_tg))
    DataFrames.select!(df_tg, ns, DataFrames.Not(ns))

    for (k, s, d) in zip(["", "_linear", "_mono", "_loop"], [psmstr_link, psmstr, psmstr_mono, psmstr_loop], [df_psm, df_linear, df_mono, df_loop])
        for K in Ks
            df_tg[!, "iden$(k)$(K)"] = map(df_tg[!, "psm$(k)$(K)_"]) do psms
                map(s, eachrow(d[psms, :])) |> xs -> join(xs, ";")
            end
            df_tg[!, "have_iden$(k)$(K)"] = .!isempty.(df_tg[!, "iden$(k)$(K)"])
            df_tg[!, "iden$(k)$(K)_credible"] = map(df_tg[!, "psm$(k)$(K)_"]) do psms
                map(s, eachrow(d[filter(i -> d.credible[i], psms), :])) |> xs -> join(xs, ";")
            end
            df_tg[!, "have_iden$(k)$(K)_credible"] = .!isempty.(df_tg[!, "iden$(k)$(K)_credible"])
        end
    end

    for K in Ks
        for (f, n) in zip([is_same_xl, is_same_xl_pepmod], ["iden", "pepmod"])
            df_tg[!, "same_$(n)$(K)"] = map(eachrow(df_tg)) do r
                filter(eachrow(df_psm[r["psm$(K)_"], :])) do s
                    (s.engine != :pLink) && return false
                    return f(r, s)
                end .|> psmstr_link |> xs -> join(xs, ";")
            end
            df_tg[!, "have_same_$(n)$(K)"] = .!isempty.(df_tg[!, "same_$(n)$(K)"])

            df_tg[!, "same_$(n)$(K)_credible"] = map(eachrow(df_tg)) do r
                filter(eachrow(df_psm[r["psm$(K)_"], :])) do s
                    (s.engine != :pLink) && return false
                    (!s.credible) && return false
                    return f(r, s)
                end .|> psmstr_link |> xs -> join(xs, ";")
            end
            df_tg[!, "have_same_$(n)$(K)_credible"] = .!isempty.(df_tg[!, "same_$(n)$(K)_credible"])

            vs = map(eachrow(df_tg)) do r
                b = nothing
                for s in eachrow(df_psm[r["psm$(K)_"], :])
                    (s.engine != :pLink) && continue
                    isnothing(b) && (b = s)
                    (s.cov_min ≤ b.cov_min) && continue
                    b = s
                end
                return isnothing(b) ? ("", false) : (psmstr_link(b), f(r, b))
            end
            df_tg[!, "best_iden$(K)"] = first.(vs)
            df_tg[!, "best_iden$(K)_is_same_$(n)"] = last.(vs)
        end
    end

    UniMZ.safe_save(p -> CSV.write(p, df_tg), joinpath(out, "$(basename(splitext(path_ms)[1])).target.TargetAcquisitionXLReport.csv"))
    UniMZ.safe_save(p -> CSV.write(p, df_psm), joinpath(out, "$(basename(splitext(path_ms)[1])).crosslink.TargetAcquisitionXLReport.csv"))
    UniMZ.safe_save(p -> CSV.write(p, df_linear), joinpath(out, "$(basename(splitext(path_ms)[1])).linear.TargetAcquisitionXLReport.csv"))
    UniMZ.safe_save(p -> CSV.write(p, df_mono), joinpath(out, "$(basename(splitext(path_ms)[1])).monolink.TargetAcquisitionXLReport.csv"))
    UniMZ.safe_save(p -> CSV.write(p, df_loop), joinpath(out, "$(basename(splitext(path_ms)[1])).looplink.TargetAcquisitionXLReport.csv"))
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="TargetAcquisitionXLReport")
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
            help = "pLink PSM path"
            required = true
        "--out"
            help = "output directory"
            default = "./out/"
            metavar = "./out/"
        "--xl"
            help = "candidate xl list"
            default = ""
        "--ft"
            help = "feature list"
            default = ""
        "--psm_pf"
            help = "pFind PSM path"
            default = ""
        "--fmt", "-f"
            help = "target list format: auto, TW, TmQE, TmFu"
            metavar = "auto|TW|TmQE|TmFu"
            default = "auto"
        "--linker"
            help = "default linker"
            default = "DSSO"
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
            help = "pLink config directory"
            default = ""
        "--cfg_pf"
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
