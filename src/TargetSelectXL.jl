module TargetSelectXL

import ArgParse
import CSV
import DataFrames
import UniMZ
import UniMZUtil: TMS, pLink

include("util.jl")

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

process(paths; df, out, name, ε, fdr_min, fdr_max, fdr_ge, fdr_le, td, pt, batch_size, rt, lc, fmt) = begin
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
    Ms = map(p -> UniMZ.read_ms(p), paths)
    build_target(df, Ms, paths, out, name, ε, batch_size, rt, lc, fmt)
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="TargetSelectXL")
    ArgParse.@add_arg_table! settings begin
        "data"
            help = "list of .umz or .ms1/2 files; .ms2/1 files should be in the same directory for .ms1/2"
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
    paths = reduce(vcat, UniMZ.match_path.(args["data"], ".umz")) |> unique |> sort
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
