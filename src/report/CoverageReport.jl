module CoverageReport

using Statistics

import ArgParse
import CSV
import DataFrames
import ProgressMeter: @showprogress
import RelocatableFolders: @path
import UniMZ: UniMZ, Plot
import UniMZUtil: UniMZUtil, pFind, pLink

const DIR_DATA = @path joinpath(@__DIR__, "../../data")

include("../util.jl")

prepare(args) = begin
    out = mkpath(args["out"])
    ε = parse(Float64, args["error"]) * 1.0e-6
    ion_syms = split(args["ion"], ",") .|> strip
    @info "selected fragment ion type: $(join(ion_syms, ", "))"
    cfg = args["cfg"]
    return (; out, ε, ion_syms, cfg)
end

process(path, paths_ms; out, ε, ion_syms, cfg) = begin
    ion_types = map(i -> getfield(UniMZ, Symbol("ion_$(i)")), ion_syms)

    M = map(p -> splitext(basename(p))[1] => UniMZ.dict_by_id(UniMZ.read_ms(p).MS2), paths_ms) |> Dict

    if isempty(cfg)
        tab_ele = pFind.read_element() |> NamedTuple
        tab_aa = map(x -> UniMZ.mass(x, tab_ele), pFind.read_amino_acid() |> NamedTuple)
        tab_mod = UniMZ.mapvalue(x -> x.mass, pFind.read_modification())
    else
        tab_ele = pFind.read_element(joinpath(cfg, "element.ini")) |> NamedTuple
        tab_aa = map(x -> UniMZ.mass(x, tab_ele), pFind.read_amino_acid(joinpath(cfg, "aa.ini")) |> NamedTuple)
        tab_mod = UniMZ.mapvalue(x -> x.mass, pFind.read_modification(joinpath(cfg, "modification.ini")))
    end

    df = pFind.read_psm(path)

    UniMZUtil.calc_cov_linear!(df, M, ε, ion_syms, ion_types, tab_ele, tab_aa, tab_mod)

    UniMZ.safe_save(p -> CSV.write(p, df), joinpath(out, basename(path) * ".CoverageReport.csv"))

    data = """
const FDR = [$(join(string.(df.fdr .* 100), ","))]
const COV_all = [$(join(string.(df.cov .* 100), ","))]
"""
    data *= map(ion_syms) do sym
"""
const COV_$(sym) = [$(join(string.(df[!, "cov_ion_$(sym)"] .* 100), ","))]
"""
    end |> join

    main = replace(read(joinpath(DIR_DATA, "CoverageReport.html"), String),
        "{{ section }}" => "Overall Peptide Fragment Ion Coverage",
        "{{ id }}" => "all",
    )
    main *= map(ion_syms) do sym
        replace(read(joinpath(DIR_DATA, "CoverageReport.html"), String),
            "{{ section }}" => "Peptide Fragment Ion Coverage of $(sym) Ions",
            "{{ id }}" => "$(sym)",
        )
    end |> join

    script = replace(read(joinpath(DIR_DATA, "CoverageReport.js"), String), "{{ id }}" => "all")
    script *= map(ion_syms) do sym
        replace(read(joinpath(DIR_DATA, "CoverageReport.js"), String), "{{ id }}" => "$(sym)")
    end |> join

    html = replace(read(joinpath(DIR_DATA, "base.html"), String),
        "{{ title }}" => "TargetWizard Fragment Ion Coverage Report",
        "{{ subtitle }}" => basename(path),
        "{{ main }}" => main,
        "{{ lib }}" => read(joinpath(DIR_DATA, "lib", "chartjs-4.2.1.js"), String),
        "{{ data }}" => data,
        "{{ script }}" => script,
    )
    path_out = joinpath(out, basename(path) * ".CoverageReport.html")
    UniMZ.safe_save(io -> write(io, html), path_out)
    UniMZ.open_url(path_out)
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="CoverageReport")
    ArgParse.@add_arg_table! settings begin
        "data"
            help = "pFind PSM path"
            required = true
        "--ms"
            help = "list of .umz or .ms2 files"
            nargs = '+'
            required = true
        "--out", "-o"
            help = "output directory"
            metavar = "./out/"
            default = "./out/"
        "--error"
            help = "m/z error"
            metavar = "ppm"
            default = "20.0"
        "--ion"
            help = "fragment ion type: a, b, c, x, y, z, a_NH3, b_NH3, y_NH3, a_H2O, b_H2O, y_H2O"
            metavar = "b,y,b_NH3,b_H2O,y_NH3,y_H2O"
            default = "b,y,b_NH3,b_H2O,y_NH3,y_H2O"
        "--cfg"
            help = "pFind config directory"
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
