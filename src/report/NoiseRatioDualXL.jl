module NoiseRatioDualXLReport

using Printf
using Statistics

import ArgParse
import CSV
import DataFrames
import ProgressMeter: @showprogress
import RelocatableFolders: @path
import UniMZ
import UniMZUtil: Crosslink, TMS, pLink

include("../util.jl")

prepare(args) = begin
    path_ms = args["ms"]
    path_psm = args["psm"]
    out = args["out"]
    mkpath(out)
    fmt = args["fmt"] |> Symbol
    linker = Symbol(args["linker"])
    fdr = parse(Float64, args["fdr"]) / 100
    ion_syms = split(args["ion"], ",") .|> strip
    ε = parse(Float64, args["error"]) * 1.0e-6
    tab_ele, tab_aa, tab_mod, tab_xl = pLink.read_mass_table(args["cfg"])
    return (; path_ms, path_psm, out, fmt, linker, fdr, ion_syms, ε, tab_ele, tab_aa, tab_mod, tab_xl)
end

process(path; path_ms, path_psm, out, fmt, linker, fdr, ion_syms, ε, tab_ele, tab_aa, tab_mod, tab_xl) = begin
    ion_types = map(i -> getfield(UniMZ, Symbol("ion_$(i)")), ion_syms)

    M = UniMZ.read_all(p -> UniMZ.dict_by_id(UniMZ.read_umz(p, split=false)), path_ms, "umz")

    dfs_psm = map(path_psm) do p
        df = pLink.read_psm_full(p; linker).xl
        df = df[df.fdr .≤ fdr, :]
        UniMZ.assign_id!(df; force=true)
        df.rt = [M[r.file][r.scan].retention_time for r in eachrow(df)]
        df.ion_ = @showprogress map(r -> Crosslink.match_crosslink_frag_ion(r, M, ε, ion_types, tab_ele, tab_aa, tab_mod, tab_xl), eachrow(df))
        Crosslink.calc_cov_crosslink!(df, df.ion_, ion_syms, ion_types)
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
    UniMZ.assign_id!(df_tg)
    TMS.parse_target_list!(df_tg, fmt)
    DataFrames.select!(df_tg, [:id, :mz, :z, :start, :stop], DataFrames.Not([:id, :mz, :z, :start, :stop]))
    "mod_a" ∈ names(df_tg) && (df_tg.mod_a = parse.(Array{UniMZ.Mod}, unify_mods_str.(df_tg.mod_a)))
    "mod_b" ∈ names(df_tg) && (df_tg.mod_b = parse.(Array{UniMZ.Mod}, unify_mods_str.(df_tg.mod_b)))

    @info "PSM mapping"
    df_tg.psm_A, df_tg.psm_B = map(dfs_psm) do df_psm
        tmp = sort!([(x.mz::Float64, x.id::Int) for x in eachrow(df_psm)])
        mzs, ids = first.(tmp), last.(tmp)
        return map(eachrow(df_tg)) do r
            # ok to use δ = 1.0 since will be further filtered
            psm = ids[UniMZ.argquery_δ(mzs, r.mz, 1.0)]
            psm = filter(i -> df_psm[i, :z] == r.z, psm)
            psm = filter(i -> r.start ≤ df_psm.rt[i] ≤ r.stop, psm)
            psm = filter(i -> is_same_xl(df_psm[i, :], r), psm)
            return psm |> sort
        end
    end

    df_tg.psm_A_best .= 0
    df_tg.psm_B_best .= 0
    df_tg.snr_num_a .= 0
    df_tg.snr_num_b .= 0
    df_tg.snr_median_a .= ""
    df_tg.snr_median_b .= ""
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
            xs, ys = getfield.(ion_A_a, :snr), getfield.(ion_B_a, :snr)
            r.snr_median_a = @sprintf("%.2f|%.2f", median(xs), median(ys))
            r.snr_mean_a = @sprintf("%.2f|%.2f", mean(xs), mean(ys))
            r.snr_sum_a = @sprintf("%.2f|%.2f", sum(xs), sum(ys))
            r.snr_min_a = @sprintf("%.2f|%.2f", minimum(xs), minimum(ys))
        end
        if r.snr_num_b != 0
            xs, ys = getfield.(ion_A_b, :snr), getfield.(ion_B_b, :snr)
            r.snr_median_b = @sprintf("%.2f|%.2f", median(xs), median(ys))
            r.snr_mean_b = @sprintf("%.2f|%.2f", mean(xs), mean(ys))
            r.snr_sum_b = @sprintf("%.2f|%.2f", sum(xs), sum(ys))
            r.snr_min_b = @sprintf("%.2f|%.2f", minimum(xs), minimum(ys))
        end
    end

    for df in [df_tg, dfs_psm...]
        DataFrames.select!(df, filter(n -> !endswith(n, '_'), names(df)))
    end

    UniMZ.safe_save(p -> CSV.write(p, df_tg), joinpath(out, "$(basename(splitext(path)[1])).NoiseRatioDualXLReport.csv"))
    UniMZ.safe_save(p -> CSV.write(p, dfs_psm[1]), joinpath(out, "$(basename(splitext(path_psm[1])[1])).NoiseRatioDualXLReport.csv"))
    UniMZ.safe_save(p -> CSV.write(p, dfs_psm[2]), joinpath(out, "$(basename(splitext(path_psm[2])[1])).NoiseRatioDualXLReport.csv"))
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="NoiseRatioDualXLReport")
    ArgParse.@add_arg_table! settings begin
        "target"
            help = "target list"
            required = true
        "--ms"
            help = "two .umz files"
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
        "--fdr"
            help = "FDR threshold (%)"
            default = "Inf"
        "--ion"
            help = "fragment ion type: a, b, c, x, y, z, a_NH3, b_NH3, y_NH3, a_H2O, b_H2O, y_H2O"
            metavar = "b,y,b_NH3,b_H2O,y_NH3,y_H2O"
            default = "b,y,b_NH3,b_H2O,y_NH3,y_H2O"
        "--error"
            help = "m/z error"
            metavar = "ppm"
            default = "20.0"
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
