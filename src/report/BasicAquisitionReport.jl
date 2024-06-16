module BasicAquisitionReport

import ArgParse
import UniMZ
import UniMZUtil: Proteomics

prepare(args) = begin
    out = args["out"]
    mkpath(out)
    return (; out)
end

process(path; out) = begin
    Proteomics.Report.basic_aquisition(path, out)
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="BasicAquisitionReport")
    ArgParse.@add_arg_table! settings begin
        "data"
            help = "list of .umz or .ms1/2 files; .ms2/1 files should be in the same directory for .ms1/2"
            nargs = '+'
            required = true
        "--out", "-o"
            help = "output directory"
            metavar = "./out/"
            default = "./out/"
    end
    args = ArgParse.parse_args(settings)
    paths = reduce(vcat, UniMZ.match_path.(args["data"], ".umz")) |> unique |> sort
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
