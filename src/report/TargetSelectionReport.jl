module TargetSelectionReport

import ArgParse
import CSV
import DataFrames
import RelocatableFolders: @path
import UniMZ
import UniMZUtil: TMS

const DIR_DATA = @path joinpath(@__DIR__, "../../data")

include("../util.jl")

prepare(args) = begin
    out = mkpath(args["out"])
    fmt = Symbol(args["fmt"])
    @info "specified format: $(fmt)"
    return (; out, fmt)
end

process(path; out, fmt) = begin
    @info "reading " * path
    df = DataFrames.DataFrame(CSV.File(path))
    TMS.parse_target_list!(df, fmt)

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
    UniMZ.safe_save(io -> write(io, html), path_out)
    UniMZ.open_url(path_out)
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="TargetSelectionReport")
    ArgParse.@add_arg_table! settings begin
        "data"
            help = "list of target files"
            nargs = '+'
            required = true
        "--out", "-o"
            help = "output directory"
            metavar = "./out/"
            default = "./out/"
        "--fmt", "-f"
            help = "target list format: auto, TW, TmQE, TmFu"
            metavar = "auto|TW|TmQE|TmFu"
            default = "auto"
    end
    args = ArgParse.parse_args(settings)
    paths = reduce(vcat, UniMZ.match_path.(args["data"], ".target.csv")) |> unique |> sort
    @info "file paths of selected data:"
    foreach(x -> println("$(x[1]):\t$(x[2])"), enumerate(paths))
    process.(paths; prepare(args)...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

julia_main()::Cint = begin
    main()
    return 0
end

end
