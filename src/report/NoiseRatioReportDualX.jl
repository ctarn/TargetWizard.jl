module NoiseRatioReportDualX

using Printf
using Statistics

import ArgParse
import CSV
import DataFrames
import ProgressMeter: @showprogress
import RelocatableFolders: @path
import UniMS
import UniMSUtil: pLink

include("../util.jl")

prepare(args) = begin
    path_ms = args["ms"]
    path_psm = args["psm"]
    out = mkpath(args["out"])
    fmt = args["fmt"] |> Symbol
    linker = Symbol(args["linker"])
    ε = parse(Float64, args["error"]) * 1.0e-6
    fdr = parse(Float64, args["fdr"]) / 100
    ion_syms = split(args["ion"], ",") .|> strip
    cfg = args["cfg"]
    return (; path_ms, path_psm, out, fmt, linker, ε, fdr, ion_syms, cfg)
end

process(path; path_ms, path_psm, out, fmt, linker, ε, fdr, ion_syms, cfg) = begin
    ion_types = map(i -> getfield(UniMS, Symbol("ion_$(i)")), ion_syms)

    M = UniMS.read_all(p -> UniMS.dict_by_id(UniMS.read_ums(p, split=false)), path_ms, "ums")

    if isempty(cfg)
        tab_ele = pLink.read_element() |> NamedTuple
        tab_aa = map(x -> UniMS.mass(x, tab_ele), pLink.read_amino_acid() |> NamedTuple)
        tab_mod = UniMS.mapvalue(x -> x.mass, pLink.read_modification())
        tab_xl = pLink.read_linker() |> NamedTuple
    else
        tab_ele = pLink.read_element(joinpath(cfg, "element.ini")) |> NamedTuple
        tab_aa = map(x -> UniMS.mass(x, tab_ele), pLink.read_amino_acid(joinpath(cfg, "aa.ini")) |> NamedTuple)
        tab_mod = UniMS.mapvalue(x -> x.mass, pLink.read_modification(joinpath(cfg, "modification.ini")))
        tab_xl = pLink.read_linker(joinpath(cfg, "xlink.ini")) |> NamedTuple
    end

    dfs_psm = map(path_psm) do p
        df = pLink.read_psm_full(p).xl
        df = df[df.fdr .≤ fdr, :]
        ("linker" ∉ names(df)) && (df.linker .= linker)
        df.id = Vector(1:size(df, 1))
        DataFrames.select!(df, :id, DataFrames.Not([:id]))
        df.rt = [M[r.file][r.scan].retention_time for r in eachrow(df)]
        calc_cov_crosslink!(df, M, ε, ion_syms, ion_types, tab_ele, tab_aa, tab_mod, tab_xl)
        df.ion_ = map(eachrow(df)) do r
            m = M[r.file][r.scan]
            map(xs -> map(x -> (; x..., snr=m.peaks[x.peak].inten / m.noises[x.peak]), xs), r.ion_)
        end
        df.snr_a = map(xs -> join(map(x -> @sprintf("%s:%.2f", x.text_abbr, x.snr), xs), ','), first.(df.ion_))
        df.snr_b = map(xs -> join(map(x -> @sprintf("%s:%.2f", x.text_abbr, x.snr), xs), ','), last.(df.ion_))
        df.db_a = map(xs -> join(map(x -> @sprintf("%s:%.2f", x.text_abbr, 10*log10(x.snr)), xs), ','), first.(df.ion_))
        df.db_b = map(xs -> join(map(x -> @sprintf("%s:%.2f", x.text_abbr, 10*log10(x.snr)), xs), ','), last.(df.ion_))
        return df
    end

    @info "Target loading from " * path
    df_tg = DataFrames.DataFrame(CSV.File(path))
    df_tg.id = Vector(1:size(df_tg, 1))
    parse_target_list!(df_tg, fmt)
    DataFrames.select!(df_tg, [:id, :mz, :z, :start, :stop], DataFrames.Not([:id, :mz, :z, :start, :stop]))
    "mod_a" ∈ names(df_tg) && (df_tg.mod_a = parse.(Array{UniMS.Mod}, unify_mods_str.(df_tg.mod_a)))
    "mod_b" ∈ names(df_tg) && (df_tg.mod_b = parse.(Array{UniMS.Mod}, unify_mods_str.(df_tg.mod_b)))

    @info "PSM mapping"
    df_tg.psm_A, df_tg.psm_B = map(dfs_psm) do df_psm
        tmp = sort!([(x.mz::Float64, x.id::Int) for x in eachrow(df_psm)])
        mzs = map(x -> x[1], tmp)
        ids = map(x -> x[2], tmp)
        return map(eachrow(df_tg)) do r
            psm = sort(filter(x -> df_psm[x, :z] == r.z, ids[UniMS.argquery_δ(mzs, r.mz, 10)]))
            psm = filter(i -> r.start ≤ df_psm.rt[i] ≤ r.stop, psm)
            psm = filter(i -> is_same_xl(df_psm[i, :], r), psm)
            return psm
        end
    end

    df_tg.psm_A_best .= 0
    df_tg.psm_B_best .= 0
    df_tg.snr_num_a .= 0
    df_tg.snr_num_b .= 0
    df_tg.snr_mean_a .= ""
    df_tg.snr_mean_b .= ""
    df_tg.snr_sum_a .= ""
    df_tg.snr_sum_b .= ""
    df_tg.snr_min_a .= ""
    df_tg.snr_min_b .= ""
    df_tg.snr_ion_a .= ""
    df_tg.snr_ion_b .= ""

    for r in eachrow(df_tg)
        r.psm_A_best = isempty(r.psm_A) ? 0 : first(r.psm_A)
        r.psm_B_best = isempty(r.psm_B) ? 0 : first(r.psm_B)
        if r.psm_A_best == 0 || r.psm_B_best == 0
            continue
        end
        ion_A_a, ion_A_b = dfs_psm[1].ion_[r.psm_A_best]
        ion_B_a, ion_B_b = dfs_psm[2].ion_[r.psm_B_best]
        s = getfield.(ion_A_a, :text) ∩ getfield.(ion_B_a, :text)
        ion_A_a = filter(i -> i.text ∈ s, ion_A_a)
        ion_B_a = filter(i -> i.text ∈ s, ion_B_a)
        s = getfield.(ion_A_b, :text) ∩ getfield.(ion_B_b, :text)
        ion_A_b = filter(i -> i.text ∈ s, ion_A_b)
        ion_B_b = filter(i -> i.text ∈ s, ion_B_b)
        sort!.([ion_A_a, ion_B_a, ion_A_b, ion_B_b]; by=x -> x.text)
        r.snr_ion_a = join([@sprintf("%s:%.2f|%.2f", x.text_abbr, x.snr, y.snr) for (x, y) in zip(ion_A_a, ion_B_a)], ",")
        r.snr_ion_b = join([@sprintf("%s:%.2f|%.2f", x.text_abbr, x.snr, y.snr) for (x, y) in zip(ion_A_b, ion_B_b)], ",")
        r.snr_num_a = length(ion_A_a)
        r.snr_num_b = length(ion_A_b)
        if r.snr_num_a != 0
            r.snr_mean_a = @sprintf("%.2f|%.2f", mean(getfield.(ion_A_a, :snr)), mean(getfield.(ion_B_a, :snr)))
            r.snr_sum_a = @sprintf("%.2f|%.2f", sum(getfield.(ion_A_a, :snr)), sum(getfield.(ion_B_a, :snr)))
            r.snr_min_a = @sprintf("%.2f|%.2f", minimum(getfield.(ion_A_a, :snr)), minimum(getfield.(ion_B_a, :snr)))
        end
        if r.snr_num_b != 0
            r.snr_mean_b = @sprintf("%.2f|%.2f", mean(getfield.(ion_A_b, :snr)), mean(getfield.(ion_B_b, :snr)))
            r.snr_sum_b = @sprintf("%.2f|%.2f", sum(getfield.(ion_A_b, :snr)), sum(getfield.(ion_B_b, :snr)))
            r.snr_min_b = @sprintf("%.2f|%.2f", maximum(getfield.(ion_A_b, :snr)), maximum(getfield.(ion_B_b, :snr)))
        end
    end

    for df in [df_tg, dfs_psm...]
        DataFrames.select!(df, filter(n -> !endswith(n, '_'), names(df)))
    end

    UniMS.safe_save(p -> CSV.write(p, df_tg), joinpath(out, "$(basename(splitext(path)[1])).NoiseRatioReportDualX.csv"))
    UniMS.safe_save(p -> CSV.write(p, dfs_psm[1]), joinpath(out, "$(basename(splitext(path_psm[1])[1])).NoiseRatioReportDualX.csv"))
    UniMS.safe_save(p -> CSV.write(p, dfs_psm[2]), joinpath(out, "$(basename(splitext(path_psm[2])[1])).NoiseRatioReportDualX.csv"))
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="NoiseRatioReportDualX")
    ArgParse.@add_arg_table! settings begin
        "target"
            help = "target list"
            required = true
        "--ms"
            help = "two .ums files"
            required = true
            nargs = 2
        "--psm"
            help = "two pLink PSM paths"
            required = true
            nargs = 2
        "--out"
            help = "output directory"
            default = "./out/"
            metavar = "./out/"
        "--fmt", "-f"
            help = "target list format: auto, TW, TmQE, TmFu"
            metavar = "auto|TW|TmQE|TmFu"
            default = "auto"
        "--linker"
            help = "default linker"
            metavar = "DSSO"
            default = "DSSO"
        "--error"
            help = "m/z error"
            metavar = "ppm"
            default = "20.0"
        "--fdr"
            help = "FDR threshold (%)"
            default = "Inf"
        "--ion"
            help = "fragment ion type: a, b, c, x, y, z, a_NH3, b_NH3, y_NH3, a_H2O, b_H2O, y_H2O"
            metavar = "b,y,b_NH3,b_H2O,y_NH3,y_H2O"
            default = "b,y,b_NH3,b_H2O,y_NH3,y_H2O"
        "--cfg"
            help = "pLink config directory"
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
