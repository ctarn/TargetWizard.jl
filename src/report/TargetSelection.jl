"""
Target Selection Report
"""
module TargetSelectionReport

"""
- target list (multiple formats supported), e.g., `TargetWizard.all.TW.target.csv`.
"""
require = true

"""
Once finished, TargetWizard will save a report to `Output Directroy`, and open the report automatically in a web browser.
- `html` report, e.g., `TargetSelectionReport.html` ([example](../assets/manual/TargetSelectionReport.html)).
"""
output = true

"""
![Target Selection Report](../assets/manual/TargetSelectionReport.png)
"""
usage = true

"""
The report will show statistics plots of the MS data, including:
- Acquisition Load
- Charge State Distribution
- Mass Distribution
- …

For many plots, you can change the settings including RT range, barplot bin size, etc.
The plot will be updated automatically.

You can also click the legend of a plot to hide or display some items.

![Acquisition Load](../assets/manual/TargetSelectionReport_AL.png)
![Mass Distribution](../assets/manual/TargetSelectionReport_MD.png)
"""
example = true

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
