module TargetSelectionReport

import ArgParse
import CSV
import DataFrames
import MesMS
import RelocatableFolders: @path

const DIR_DATA = @path joinpath(@__DIR__, "../../data")

parse_target_list!(df, fmt) = begin
    if fmt == :auto
        cols_TW = ["mz", "z", "start", "stop"]
        cols_TmQE = ["Mass [m/z]", "CS [z]", "Start [min]", "End [min]"]
        cols_TmFu = ["m/z", "z", "t start (min)", "t stop (min)"]
        if all(n -> n ∈ names(df), cols_TW)
            @info "list treated as `TargetWizard` format based on following columns: $(cols_TW)"
            fmt = :TW
        elseif all(n -> n ∈ names(df), cols_TmQE)
            @info "list treated as `Thermo Q Exactive` format based on following columns: $(cols_TmQE)"
            fmt = :TmQE
        elseif all(n -> n ∈ names(df), cols_TmFu)
            @info "list treated as `Thermo Fusion` format based on following columns: $(cols_TmFu)"
            fmt = :TmFu
        else
            @error "failed to detect list format"
            fmt = :unknown
        end
    end
    if fmt == :TmQE
        DataFrames.rename!(df, "Mass [m/z]" => "mz", "CS [z]" => "z", "Start [min]" => "start", "End [min]" => "stop")
        df.start = df.start .* 60
        df.stop = df.stop .* 60
    elseif fmt == :TmFu
        DataFrames.rename!(df, "m/z" => "mz", "t start (min)" => "start", "t stop (min)" => "stop")
        df.start = df.start .* 60
        df.stop = df.stop .* 60
    end
    df.m = MesMS.mz_to_m.(df.mz, df.z)
    df.rt = (df.start .+ df.stop) ./ 2
    return df
end

prepare(args) = begin
    fmt = Symbol(args["fmt"])
    @info "specified format: $(fmt)"
    out = mkpath(args["out"])
    return (; fmt, out)
end

plot(path; fmt, out) = begin
    @info "reading " * path
    df = DataFrames.DataFrame(CSV.File(path))
    parse_target_list!(df, fmt)
    name = basename(path)

    data = """
const M = [$(join(string.(df.m), ","))]
const M_MIN = $(minimum(df.m))
const M_MAX = $(maximum(df.m))
const MZ = [$(join(string.(df.mz), ","))]
const MZ_MIN = $(minimum(df.mz))
const MZ_MAX = $(maximum(df.mz))
const Z = [$(join(string.(df.z), ","))]
const Z_MIN = $(minimum(df.z))
const Z_MAX = $(maximum(df.z))
const RT = [$(join(string.(df.rt), ","))]
const RT_MIN = $(minimum(df.rt))
const RT_MAX = $(maximum(df.rt))
const RT_START = [$(join(string.(df.start), ","))]
const RT_START_MIN = $(minimum(df.start))
const RT_START_MAX = $(maximum(df.start))
const RT_STOP = [$(join(string.(df.stop), ","))]
const RT_STOP_MIN = $(minimum(df.stop))
const RT_STOP_MAX = $(maximum(df.stop))
"""

    html = replace(read(joinpath(DIR_DATA, "base.html"), String),
        "{{ title }}" => "TargetWizard Target Report",
        "{{ subtitle }}" => basename(path),
        "{{ main }}" => read(joinpath(DIR_DATA, "TargetSelectionReport.html"), String),
        "{{ lib }}" => read(joinpath(DIR_DATA, "lib", "chartjs-4.2.1.js"), String),
        "{{ data }}" => data,
        "{{ script }}" => read(joinpath(DIR_DATA, "TargetSelectionReport.js"), String),
    )
    path_out = joinpath(out, basename(path) * ".TargetSelectionReport.html")
    MesMS.safe_save(io -> write(io, html), path_out)
    MesMS.open_url(path_out)
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="TargetSelectionReport")
    ArgParse.@add_arg_table! settings begin
        "--fmt", "-f"
            help = "target list format: auto, TW, TmQE, TmFu"
            metavar = "auto|TW|TmQE|TmFu"
            default = "auto"
        "--out", "-o"
            help = "output directory"
            metavar = "./out/"
            default = "./out/"
        "data"
            help = "list of target files"
            nargs = '+'
            required = true
    end
    args = ArgParse.parse_args(settings)
    paths = (sort∘unique∘reduce)(vcat, MesMS.match_path.(args["data"], ".target.csv"); init=String[])
    @info "file paths of selected data:"
    foreach(x -> println("$(x[1]):\t$(x[2])"), enumerate(paths))
    sess = prepare(args)
    plot.(paths; sess...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

julia_main()::Cint = begin
    main()
    return 0
end

end
