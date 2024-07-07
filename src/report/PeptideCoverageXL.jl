"""
Peptide Fragment Ion Coverage Report for Crosslink
"""
module PeptideCoverageXLReport

"""
- (filtered) identification results of targeted mass spectrometry data, e.g., `TMS_fdr.pfind.csv`.
- traditional or targeted mass spectrometry data, e.g., `TMS.raw`.
  - The raw data should be converted into an open-source format such as MS1/MS2. [ThermoRawRead](http://thermorawread.ctarn.io) is recommended.
"""
require = true

"""
Once finished, TargetWizard will save reports to `Output Directroy`, and open a `html` report automatically in a web browser.
- `html` report for linear peptides, e.g., `PeptideCoverageReport.html` ([example](../assets/report/PeptideCoverageReport.html)).
- `html` report for crosslinks, e.g., `PeptideCoverageXLReport.html` ([example](../assets/report/PeptideCoverageXLReport.html)).
- `csv` report including detailed fragment information for linear, monolink, looplink, and crosslink peptides.
"""
output = true

"""
![Peptide Coverage Report for Crosslink](../assets/report/PeptideCoverageReport.png)
"""
usage = true

"""
The report will show statistics plots of fragment ion coverage of peptide pairs or α/β peptides, including:
- Overall Peptide Coverage Distribution
- Peptide Coverage Distribution of Specific Ion Types
- …

For many plots, you can change the settings including barplot bin size, FDR, etc.
The plot will be updated automatically.

You can also click the legend of a plot to hide or display some items.

![Cross-linked Peptide Fragment Ion Coverage of b Ions](../assets/report/PeptideCoverageXLReport_b.png)
"""
example = true

import ArgParse
import UniMZ
import UniMZUtil: Crosslink, pLink

prepare(args) = begin
    out = args["out"]
    mkpath(out)
    linker = Symbol(args["linker"])
    ε = parse(Float64, args["error"]) * 1.0e-6
    ion_syms = split(args["ion"], ",") .|> strip
    @info "selected fragment ion type: $(join(ion_syms, ", "))"
    tab_ele, tab_aa, tab_mod, tab_xl = pLink.read_mass_table(args["cfg"])
    return (; out, linker, ε, ion_syms, tab_ele, tab_aa, tab_mod, tab_xl)
end

process(path, paths_ms; out, linker, ε, ion_syms, tab_ele, tab_aa, tab_mod, tab_xl) = begin
    M = map(p -> splitext(basename(p))[1] => UniMZ.dict_by_id(UniMZ.read_ms(p).MS2), paths_ms) |> Dict
    dfs = pLink.read_psm_full(path; linker)
    Crosslink.Report.peptide_coverage_linear(dfs.linear, M, joinpath(out, basename(path)); name=path, ε, ion_syms, tab_ele, tab_aa, tab_mod)
    Crosslink.Report.peptide_coverage_monolink(dfs.mono, M, joinpath(out, basename(path)); name=path, ε, ion_syms, tab_ele, tab_aa, tab_mod, tab_xl)
    Crosslink.Report.peptide_coverage_looplink(dfs.loop, M, joinpath(out, basename(path)); name=path, ε, ion_syms, tab_ele, tab_aa, tab_mod, tab_xl)
    Crosslink.Report.peptide_coverage_crosslink(dfs.xl, M, joinpath(out, basename(path)); name=path, ε, ion_syms, tab_ele, tab_aa, tab_mod, tab_xl)
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="PeptideCoverageXLReport")
    ArgParse.@add_arg_table! settings begin
        "data"
            help = "pLink PSM path"
            required = true
        "--ms"
            help = "list of .umz or .ms2 files"
            nargs = '+'
            required = true
        "--out", "-o"
            help = "output directory"
            metavar = "./out/"
            default = "./out/"
        "--linker"
            help = "default linker"
            default = "DSSO"
        "--error"
            help = "m/z error"
            metavar = "ppm"
            default = "20.0"
        "--ion"
            help = "fragment ion type: a, b, c, x, y, z, a_NH3, b_NH3, y_NH3, a_H2O, b_H2O, y_H2O"
            metavar = "b,y,b_NH3,b_H2O,y_NH3,y_H2O"
            default = "b,y,b_NH3,b_H2O,y_NH3,y_H2O"
        "--cfg"
            help = "pLink config directory"
            default = ""
    end
    args = ArgParse.parse_args(settings)
    paths = reduce(vcat, UniMZ.match_path.(args["ms"], ".umz")) |> unique |> sort
    @info "file paths of selected data:"
    foreach(x -> println("$(x[1]):\t$(x[2])"), enumerate(paths))
    process(args["data"], paths; prepare(args)...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

julia_main()::Cint = begin
    main()
    return 0
end

end
