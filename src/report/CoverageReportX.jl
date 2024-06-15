module CoverageReportX

using Statistics

import ArgParse
import CSV
import DataFrames
import ProgressMeter: @showprogress
import RelocatableFolders: @path
import UniMZ: UniMZ, Plot
import UniMZUtil: Crosslink, pFind, pLink

const DIR_DATA = @path joinpath(@__DIR__, "../../data")

include("../util.jl")

prepare(args) = begin
    out = mkpath(args["out"])
    linker = Symbol(args["linker"])
    ε = parse(Float64, args["error"]) * 1.0e-6
    ion_syms = split(args["ion"], ",") .|> strip
    @info "selected fragment ion type: $(join(ion_syms, ", "))"
    cfg = args["cfg"]
    return (; out, linker, ε, ion_syms, cfg)
end

process(path, paths_ms; out, linker, ε, ion_syms, cfg) = begin
    ion_types = map(i -> getfield(UniMZ, Symbol("ion_$(i)")), ion_syms)

    M = map(p -> splitext(basename(p))[1] => UniMZ.dict_by_id(UniMZ.read_ms(p).MS2), paths_ms) |> Dict

    if isempty(cfg)
        tab_ele = pLink.read_element() |> NamedTuple
        tab_aa = map(x -> UniMZ.mass(x, tab_ele), pLink.read_amino_acid() |> NamedTuple)
        tab_mod = UniMZ.mapvalue(x -> x.mass, pLink.read_modification())
        tab_xl = pLink.read_linker() |> NamedTuple
    else
        tab_ele = pLink.read_element(joinpath(cfg, "element.ini")) |> NamedTuple
        tab_aa = map(x -> UniMZ.mass(x, tab_ele), pLink.read_amino_acid(joinpath(cfg, "aa.ini")) |> NamedTuple)
        tab_mod = UniMZ.mapvalue(x -> x.mass, pLink.read_modification(joinpath(cfg, "modification.ini")))
        tab_xl = pLink.read_linker(joinpath(cfg, "xlink.ini")) |> NamedTuple
    end

    psm = pLink.read_psm_full(path)
    df_xl = psm.xl
    df_linear = psm.linear
    df_mono = psm.mono
    df_loop = psm.loop
    df_xl.linker .= linker
    df_mono.linker .= linker
    df_loop.linker .= linker

    Crosslink.calc_cov_crosslink!(df_xl, M, ε, ion_syms, ion_types, tab_ele, tab_aa, tab_mod, tab_xl)
    Crosslink.calc_cov_linear!(df_linear, M, ε, ion_syms, ion_types, tab_ele, tab_aa, tab_mod)
    Crosslink.calc_cov_monolink!(df_mono, M, ε, ion_syms, ion_types, tab_ele, tab_aa, tab_mod, tab_xl)
    Crosslink.calc_cov_looplink!(df_loop, M, ε, ion_syms, ion_types, tab_ele, tab_aa, tab_mod, tab_xl)

    UniMZ.safe_save(p -> CSV.write(p, df_xl), joinpath(out, basename(path) * ".crosslink.CoverageReportX.csv"))
    UniMZ.safe_save(p -> CSV.write(p, df_linear), joinpath(out, basename(path) * ".linear.CoverageReportX.csv"))
    UniMZ.safe_save(p -> CSV.write(p, df_mono), joinpath(out, basename(path) * ".monolink.CoverageReportX.csv"))
    UniMZ.safe_save(p -> CSV.write(p, df_loop), joinpath(out, basename(path) * ".looplink.CoverageReportX.csv"))

    data = """
const FDR = [$(join(string.(df_xl.fdr .* 100), ","))]
const COV_all = [$(join(string.(df_xl.cov .* 100), ","))]
const COV_A_all = [$(join(string.(df_xl.cov_a .* 100), ","))]
const COV_B_all = [$(join(string.(df_xl.cov_b .* 100), ","))]
"""
    data *= map(ion_syms) do sym
"""
const COV_$(sym) = [$(join(string.(df_xl[!, "cov_ion_$(sym)"] .* 100), ","))]
const COV_A_$(sym) = [$(join(string.(df_xl[!, "cov_a_ion_$(sym)"] .* 100), ","))]
const COV_B_$(sym) = [$(join(string.(df_xl[!, "cov_b_ion_$(sym)"] .* 100), ","))]
"""
    end |> join

    main = replace(read(joinpath(DIR_DATA, "CoverageReportX.html"), String),
        "{{ section }}" => "Overall Peptide Fragment Ion Coverage",
        "{{ id }}" => "all",
    )
    main *= map(ion_syms) do sym
        replace(read(joinpath(DIR_DATA, "CoverageReportX.html"), String),
            "{{ section }}" => "Peptide Fragment Ion Coverage of $(sym) Ions",
            "{{ id }}" => "$(sym)",
        )
    end |> join

    script = replace(read(joinpath(DIR_DATA, "CoverageReportX.js"), String), "{{ id }}" => "all")
    script *= map(ion_syms) do sym
        replace(read(joinpath(DIR_DATA, "CoverageReportX.js"), String), "{{ id }}" => "$(sym)")
    end |> join

    html = replace(read(joinpath(DIR_DATA, "base.html"), String),
        "{{ title }}" => "TargetWizard Cross-link Fragment Ion Coverage Report",
        "{{ subtitle }}" => basename(path),
        "{{ main }}" => main,
        "{{ lib }}" => read(joinpath(DIR_DATA, "lib", "chartjs-4.2.1.js"), String),
        "{{ data }}" => data,
        "{{ script }}" => script,
    )
    path_out = joinpath(out, basename(path) * ".CoverageReportX.html")
    UniMZ.safe_save(io -> write(io, html), path_out)
    UniMZ.open_url(path_out)
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="CoverageReportX")
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
