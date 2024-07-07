"""
Basic Data Acquisition Report
"""
module BasicAcquisitionReport

"""
- traditional or targeted mass spectrometry data, e.g., `TMS.raw`.
  - The raw data should be converted into an open-source format such as MS1/MS2. [ThermoRawRead](http://thermorawread.ctarn.io) is recommended.
"""
require = true

"""
Once finished, TargetWizard will save a report to `Output Directroy`, and open the report automatically in a web browser.
- `html` report, e.g., `BasicAcquisitionReport.html` ([example](../assets/manual/BasicAcquisitionReport.html)).
"""
output = true

"""
![Basic Acquisition Report](../assets/manual/BasicAcquisitionReport.png)
"""
usage = true

"""
The report will show statistics plots of the MS data, including:
- Acquisition Speed
- Activation Center
- Ion Injection Time
- Total Ion Current (TIC), Base Peak Intensity (BPI), Base Peak Mass (BPM)
- …

For many plots, you can change the settings including RT range, barplot bin size, etc.
The plot will be updated automatically.

You can also click the legend of a plot to hide or display some items.

![Acquisition Speed](../assets/manual/BasicAcquisitionReport_AS.png)
![Ion Injection Time](../assets/manual/BasicAcquisitionReport_IT.png)
"""
example = true

import ArgParse
import UniMZ
import UniMZUtil: Proteomics

prepare(args) = begin
    out = args["out"]
    mkpath(out)
    return (; out)
end

process(path; out) = begin
    Proteomics.Report.basic_acquisition(path, out)
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="BasicAcquisitionReport")
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
