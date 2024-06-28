module TargetSelectionReport

import ArgParse
import CSV
import DataFrames
import UniMZ
import UniMZUtil: TMS

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
    TMS.Report.target_selection(df, joinpath(out, basename(path)); name=path)
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
