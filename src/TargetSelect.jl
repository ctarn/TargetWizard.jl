module TargetSelect

import ArgParse
import CSV
import DataFrames
import MesMS
import MesUtil: pFind, pLink

build_target_hf(df) = begin
    return DataFrames.DataFrame(
        Symbol("Mass [m/z]") => df.mz,
        Symbol("Formula [M]") => "",
        Symbol("Formula type") => "",
        Symbol("Species") => "",
        Symbol("CS [z]") => df.z,
        Symbol("Polarity") => "Positive",
        Symbol("Start [min]") => df.start ./ 60,
        Symbol("End [min]") => df.stop ./ 60,
        Symbol("(N)CE") => "",
        Symbol("(N)CE type") => "",
        Symbol("MSX ID") => "",
        Symbol("Comment") => "",
    )
end

build_target_lumos(df) = begin
    return DataFrames.DataFrame(
        Symbol("Compound") => "",
        Symbol("Formula") => "",
        Symbol("Adduct") => "",
        Symbol("m/z") => df.mz,
        Symbol("z") => df.z,
        Symbol("t start (min)") => df.start ./ 60,
        Symbol("t stop (min)") => df.stop ./ 60,
    )
end

prepare(args) = begin
    name = args["name"]
    df = pLink.read_psm_full(args["psm"]).xl
    df.mod_a = pFind.modstr.(df.mod_a)
    df.mod_b = pFind.modstr.(df.mod_b)
    fdr_min = parse(Float64, args["fdr_min"]) / 100
    fdr_max = parse(Float64, args["fdr_max"]) / 100
    fdr_min_close = args["fdr_min_close"]
    fdr_max_close = args["fdr_max_close"]
    td = split(args["td_type"], ",") .|> strip .|> Symbol
    pt = split(args["prot_type"], ",") .|> strip .|> Symbol
    batch_size = parse(Float64, args["batch"])
    rt = parse(Float64, args["rtime"])
    rt_max = parse(Float64, args["rtime_max"])
    out = mkpath(args["out"])
    return (; name, df, fdr_min, fdr_max, fdr_min_close, fdr_max_close, td, pt, batch_size, rt, rt_max, out)
end

target_select(paths; name, df, fdr_min, fdr_max, fdr_min_close, fdr_max_close, td, pt, batch_size, rt, rt_max, out) = begin
    Ms = map(paths) do path
        return MesMS.mapvalue(MesMS.dict_by_id, MesMS.read_all(MesMS.read_ms2, path))
    end
    M = reduce(merge, Ms)

    s = trues(size(df, 1))
    s .&= fdr_min_close ? (df.fdr .≥ fdr_min) : (df.fdr .> fdr_min)
    s .&= fdr_max_close ? (df.fdr .≤ fdr_max) : (df.fdr .< fdr_max)
    s .&= reduce(.|, [df.td .== t for t in td])
    s .&= reduce(.|, [df.prot_type .== t for t in pt])

    df = df[s, :]

    gd = DataFrames.groupby(df, [:pep_a, :pep_b, :site_a, :site_b, :mod_a, :mod_b, :z])
    df = DataFrames.combine(gd,
        [:td, :prot_type, :score, :mz, :mz_calc, :scan, :raw] .=> first,
        renamecols=false,
    )

    df.rt = [M[r.raw][r.scan].retention_time for r in eachrow(df)]

    df.start = min.(rt_max * 60, max.(0, df.rt .- (rt / 2)))
    df.stop = min.(rt_max * 60, max.(0, df.rt .+ (rt / 2)))

    n_batch = isinf(batch_size) ? 1 : ceil(Int, size(df, 1) / batch_size)
    df = sort(df, :rt)
    df.id = 1:size(df, 1)
    df.batch = (df.id .- 1) .% n_batch .+ 1

    @info "$(size(df, 1)) features splitting into $(n_batch) batches"

    path_pre = joinpath(out, name)
    MesMS.safe_save(p -> CSV.write(p, df), "$(path_pre).generic.all.aims.tw.csv", "AIMS list")

    for i in 1:n_batch
        df_ = df[df.batch .== i, :]
        @info "batch $(i): $(size(df_, 1))"
        MesMS.safe_save(p -> CSV.write(p, df_), "$(path_pre).generic.batch$(i).aims.tw.csv", "AIMS list")
        MesMS.safe_save(p -> CSV.write(p, build_target_lumos(df_)), "$(path_pre).lumos.batch$(i).aims.tw.csv", "AIMS (Lumos) list")
        MesMS.safe_save(p -> CSV.write(p, build_target_hf(df_)), "$(path_pre).hf.batch$(i).aims.tw.csv", "AIMS (HF) list")
    end
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="XeekAIMS")
    ArgParse.@add_arg_table! settings begin
        "--name"
            help = "task name"
            metavar = "name"
            default = "TargetWizard"
        "--psm"
            help = "pLink PSM file (full list)"
            metavar = "PSM"
            required = true
        "--fdr_min"
            help = "min. FDR (%)"
            metavar = "min"
            default = "-Inf"
        "--fdr_max"
            help = "max. FDR (%)"
            metavar = "max"
            default = "Inf"
        "--fdr_min_close"
            help = "include min. FDR (compared with `≥ min` instead of `> min`)"
            action = :store_true
        "--fdr_max_close"
            help = "include max. FDR (compared with `≤ max` instead of `< max`)"
            action = :store_true
        "--td_type"
            help = "target/decoy types (split by `,`)"
            metavar = "td"
            default = "TT,TD,DD"
        "--prot_type"
            help = "inter/intra types (split by `,`)"
            metavar = "pt"
            default = "Inter,Intra"
        "--batch", "-b"
            help = "batch size"
            metavar = "num"
            default = "Inf"
        "--rtime", "-t"
            help = "retention time window"
            metavar = "sec"
            default = "240"
        "--rtime_max"
            help = "max. retention time"
            metavar = "min"
            default = "Inf"
        "--out", "-o"
            help = "output directory"
            metavar = "output"
            default = "./out/"
        "data"
            help = "list of .MS2 files"
            nargs = '+'
            required = true
    end
    args = ArgParse.parse_args(settings)
    paths = (sort∘unique∘reduce)(vcat, MesMS.match_path.(args["data"], ".ms2"); init=String[])
    @info "file paths of selected data:"
    foreach(x -> println("$(x[1]):\t$(x[2])"), enumerate(paths))
    sess = prepare(args)
    target_select(paths; sess...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

julia_main()::Cint = begin
    main()
    return 0
end

end
