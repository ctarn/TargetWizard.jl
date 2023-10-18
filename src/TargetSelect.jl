module TargetSelect

import ArgParse
import CSV
import DataFrames
import MesMS
import MesUtil: TMS, pFind, pLink

prepare(args) = begin
    df = pLink.read_psm_full(args["psm"]).xl
    out = mkpath(args["out"])
    name = args["name"]
    ε = parse(Float64, args["error"]) * 1.0e-6
    fdr_min = parse(Float64, args["fdr_min"]) / 100
    fdr_max = parse(Float64, args["fdr_max"]) / 100
    fdr_ge = args["fdr_ge"]
    fdr_le = args["fdr_le"]
    td = split(args["td"], ",") .|> strip .|> Symbol
    pt = split(args["pt"], ",") .|> strip .|> Symbol
    batch_size = parse(Float64, args["batch"])
    rt = parse(Float64, args["rtime"])
    lc = parse(Float64, args["lc"])
    fmt = split(args["fmt"], ",") .|> strip .|> Symbol
    return (; df, out, name, ε, fdr_min, fdr_max, fdr_ge, fdr_le, td, pt, batch_size, rt, lc, fmt)
end

_nearbymax(M, i, mz, ε, τ, δ) = begin
    skip = 0
    v = -Inf
    i_max = i
    while skip ≤ τ
        v_ = MesMS.max_inten_ε(M[i].peaks, mz, ε)
        if v_ > v
            skip = 0
            v = v_
            i_max = i
        else
            skip += 1
        end
        i += δ
    end
    return i_max, v
end

nearbymax(M, i, mz, ε, τ=2) = begin
    li, lv = _nearbymax(M, i, mz, ε, τ, -1)
    ri, rv = _nearbymax(M, i, mz, ε, τ, 1)
    return lv ≥ rv ? (li, lv) : (ri, rv)
end

process(paths; df, out, name, ε, fdr_min, fdr_max, fdr_ge, fdr_le, td, pt, batch_size, rt, lc, fmt) = begin
    Ms = map(p -> MesMS.read_ms(p), paths)
    M1 = map((p, M) -> splitext(basename(p))[1] => MesMS.dict_by_id(M.MS1), paths, Ms) |> Dict
    M2 = map((p, M) -> splitext(basename(p))[1] => MesMS.dict_by_id(M.MS2), paths, Ms) |> Dict
    M1V = map((p, M) -> splitext(basename(p))[1] => M.MS1, paths, Ms) |> Dict
    M1I = map((p, M) -> splitext(basename(p))[1] => [m.id => i for (i, m) in enumerate(M.MS1)] |> Dict, paths, Ms) |> Dict

    s = trues(size(df, 1))
    s .&= fdr_ge ? (df.fdr .≥ fdr_min) : (df.fdr .> fdr_min)
    s .&= fdr_le ? (df.fdr .≤ fdr_max) : (df.fdr .< fdr_max)
    s .&= reduce(.|, [df.td .== t for t in td])
    s .&= reduce(.|, [df.prot_type .== t for t in pt])

    df = df[s, :]

    gd = DataFrames.groupby(df, [:pep_a, :pep_b, :site_a, :site_b, :mod_a, :mod_b, :z])
    df = DataFrames.combine(gd,
        [:td, :prot_type, :fdr, :score, :mz, :mz_calc, :scan, :file] .=> first,
        renamecols=false,
    )

    df.rt = [M2[r.file][r.scan].retention_time for r in eachrow(df)]

    df.inten = map(eachrow(df)) do r
        m2 = M2[r.file][r.scan]
        m1 = M1[r.file][m2.pre]
        return MesMS.max_inten_ε(m1.peaks, r.mz, ε)
    end
    vs = map(eachrow(df)) do r
        m2 = M2[r.file][r.scan]
        i = M1I[r.file][m2.pre]
        i, v = nearbymax(M1V[r.file], i, r.mz, ε, 2)
        return v, (M1V[r.file][i].retention_time - m2.retention_time)
    end
    df.inten_max = first.(vs)
    df.inten_max_delta_rt = last.(vs)

    df.start = min.(lc * 60, max.(0, df.rt .- (rt / 2)))
    df.stop = min.(lc * 60, max.(0, df.rt .+ (rt / 2)))

    n_batch = isinf(batch_size) ? 1 : ceil(Int, size(df, 1) / batch_size)
    df = sort(df, :rt)
    df.id = 1:size(df, 1)
    df.batch = (df.id .- 1) .% n_batch .+ 1

    @info "$(size(df, 1)) features splitting into $(n_batch) batches"

    tw = :TW ∈ fmt
    tmqe = :TmQE ∈ fmt
    tmfu = :TmFu ∈ fmt
    p = joinpath(out, name)
    tw && MesMS.safe_save(p -> CSV.write(p, df), "$(p).all.TW.target.csv", "list")

    for i in 1:n_batch
        df_ = df[df.batch .== i, :]
        @info "batch $(i): $(size(df_, 1))"
        tw && MesMS.safe_save(p -> CSV.write(p, df_), "$(p).batch$(i).TW.target.csv", "list")
        tmqe && MesMS.safe_save(p -> CSV.write(p, TMS.build_target_TmQE(df_)), "$(p).batch$(i).TmQE.target.csv", "list (Thermo Q Exactive)")
        tmfu && MesMS.safe_save(p -> CSV.write(p, TMS.build_target_TmFu(df_)), "$(p).batch$(i).TmFu.target.csv", "list (Thermo Fusion)")
    end
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="TargetSelect")
    ArgParse.@add_arg_table! settings begin
        "data"
            help = "list of .mes or .ms1/2 files; .ms2/1 files should be in the same directory for .ms1/2"
            nargs = '+'
            required = true
        "--psm"
            help = "pLink PSM file (full list)"
            metavar = "PSM"
            required = true
        "--out", "-o"
            help = "output directory"
            metavar = "./out/"
            default = "./out/"
        "--name"
            help = "task name"
            metavar = "name"
            default = "TargetWizard"
        "--error"
            help = "m/z error"
            metavar = "ppm"
            default = "20.0"
        "--fdr_min"
            help = "min. FDR (%)"
            metavar = "min"
            default = "-Inf"
        "--fdr_max"
            help = "max. FDR (%)"
            metavar = "max"
            default = "Inf"
        "--fdr_ge"
            help = "include min. FDR (compared with `≥ min` instead of `> min`)"
            action = :store_true
        "--fdr_le"
            help = "include max. FDR (compared with `≤ max` instead of `< max`)"
            action = :store_true
        "--td"
            help = "target/decoy types (split by `,`)"
            metavar = "TT,TD,DD"
            default = "TT,TD,DD"
        "--pt"
            help = "inter/intra-protein types (split by `,`)"
            metavar = "Inter,Intra"
            default = "Inter,Intra"
        "--batch"
            help = "batch size"
            metavar = "num"
            default = "Inf"
        "--rtime"
            help = "retention time window (sec)"
            metavar = "sec"
            default = "240"
        "--lc"
            help = "LC gradient length (min)"
            metavar = "min"
            default = "Inf"
        "--fmt"
            help = "format(s) of target list (split by `,`)        TW: TargetWizard, TmQE: Thermo Q Exactive, TmFu: Thermo Fusion"
            metavar = "TW,TmQE,TmFu"
            default = "TW,TmQE,TmFu"
    end
    args = ArgParse.parse_args(settings)
    paths = reduce(vcat, MesMS.match_path.(args["data"], ".mes")) |> unique |> sort
    @info "file paths of selected data:"
    foreach(x -> println("$(x[1]):\t$(x[2])"), enumerate(paths))
    process(paths; prepare(args)...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

julia_main()::Cint = begin
    main()
    return 0
end

end
