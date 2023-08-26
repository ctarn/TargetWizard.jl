module XLCoverageReport

using Statistics

import ArgParse
import CSV
import DataFrames
import MesMS: MesMS, Plot
import MesUtil: pFind, pLink
import ProgressMeter: @showprogress
import RelocatableFolders: @path

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

process_crosslink!(df, M, ε, ion_syms, ion_types, tab_ele, tab_aa, tab_mod, tab_xl) = begin
    @info "fragment ion calculating"
    df.ion = @showprogress map(eachrow(df)) do r
        peaks = M[r.file][r.scan].peaks
        seqs = (r.pep_a, r.pep_b)
        modss = (r.mod_a, r.mod_b)
        sites = (r.site_a, r.site_b)
        types = [(i, 1:(r.z-1)) for i in ion_types]
        ionss = Plot.build_ions_crosslink(peaks, seqs, modss, tab_xl[r.linker], sites, ε, tab_ele, tab_aa, tab_mod; types)
        return filter.(i -> i.peak > 0 && i.loc > 0, ionss)
    end

    @info "coverage calculating"
    df.ion_a = first.(df.ion)
    df.ion_b = last.(df.ion)
    foreach(r -> filter!(i -> i.loc < length(r.pep_a), r.ion_a), eachrow(df))
    foreach(r -> filter!(i -> i.loc < length(r.pep_b), r.ion_b), eachrow(df))

    match_a = [falses(length(r.pep_a)-1) for r in eachrow(df)]
    match_b = [falses(length(r.pep_b)-1) for r in eachrow(df)]
    for idx in 1:size(df, 1)
        foreach(i -> match_a[idx][i.loc] = true, df.ion_a[idx])
        foreach(i -> match_b[idx][i.loc] = true, df.ion_b[idx])
    end

    df.cov = round.(mean.(vcat.(match_a, match_b)); digits=4)
    df.cov_a = round.(mean.(match_a); digits=4)
    df.cov_b = round.(mean.(match_b); digits=4)

    @info "coverage of each type of fragment ion calculating"
    @showprogress for (sym, ion_type) in zip(ion_syms, ion_types)
        ion_a = map(r -> filter(i -> i.type == ion_type.type, r.ion_a), eachrow(df))
        ion_b = map(r -> filter(i -> i.type == ion_type.type, r.ion_b), eachrow(df))
        match_a = [falses(length(r.pep_a)-1) for r in eachrow(df)]
        match_b = [falses(length(r.pep_b)-1) for r in eachrow(df)]
        for idx in 1:size(df, 1)
            foreach(i -> match_a[idx][i.loc] = true, ion_a[idx])
            foreach(i -> match_b[idx][i.loc] = true, ion_b[idx])
        end
        df[!, "cov_ion_$(sym)"] = round.(mean.(vcat.(match_a, match_b)); digits=4)
        df[!, "cov_a_ion_$(sym)"] = round.(mean.(match_a); digits=4)
        df[!, "cov_b_ion_$(sym)"] = round.(mean.(match_b); digits=4)
    end

    df.prot_a = pFind.protstr.(df.prot_a)
    df.prot_b = pFind.protstr.(df.prot_b)

    DataFrames.select!(df, DataFrames.Not([:ion, :ion_a, :ion_b]), :ion_a, :ion_b)
    df.ion_a = map(ions -> join(getfield.(ions, :text_abbr), ','), df.ion_a)
    df.ion_b = map(ions -> join(getfield.(ions, :text_abbr), ','), df.ion_b)

    return df
end

process_linear!(df, M, ε, ion_syms, ion_types, tab_ele, tab_aa, tab_mod) = begin
    @info "fragment ion calculating"
    df.ion = @showprogress map(eachrow(df)) do r
        peaks = M[r.file][r.scan].peaks
        types = [(i, 1:(r.z-1)) for i in ion_types]
        ions = Plot.build_ions(peaks, r.pep, r.mod, ε, tab_ele, tab_aa, tab_mod; types)
        return filter(i -> i.peak > 0 && 0 < i.loc < length(r.pep), ions)
    end

    @info "coverage calculating"
    match = [falses(length(r.pep)-1) for r in eachrow(df)]
    for idx in 1:size(df, 1)
        foreach(i -> match[idx][i.loc] = true, df.ion[idx])
    end

    df.cov = round.(mean.(match); digits=4)

    @info "coverage of each type of fragment ion calculating"
    @showprogress for (sym, ion_type) in zip(ion_syms, ion_types)
        ion = map(r -> filter(i -> i.type == ion_type.type, r.ion), eachrow(df))
        match = [falses(length(r.pep)-1) for r in eachrow(df)]
        for idx in 1:size(df, 1)
            foreach(i -> match[idx][i.loc] = true, ion[idx])
        end
        df[!, "cov_ion_$(sym)"] = round.(mean.(match); digits=4)
    end

    DataFrames.select!(df, DataFrames.Not([:ion]), :ion)
    df.ion = map(ions -> join(getfield.(ions, :text_abbr), ','), df.ion)

    return df
end

process(path, paths_ms; out, linker, ε, ion_syms, cfg) = begin
    ion_types = map(i -> getfield(MesMS, Symbol("ion_$(i)")), ion_syms)

    M = map(p -> splitext(basename(p))[1] => MesMS.dict_by_id(MesMS.read_ms(p).MS2), paths_ms) |> Dict

    if isempty(cfg)
        tab_ele = pLink.read_element() |> NamedTuple
        tab_aa = map(x -> MesMS.mass(x, tab_ele), pLink.read_amino_acid() |> NamedTuple)
        tab_mod = MesMS.mapvalue(x -> x.mass, pLink.read_modification())
        tab_xl = pLink.read_linker() |> NamedTuple
    else
        tab_ele = pLink.read_element(joinpath(cfg, "element.ini")) |> NamedTuple
        tab_aa = map(x -> MesMS.mass(x, tab_ele), pLink.read_amino_acid(joinpath(cfg, "aa.ini")) |> NamedTuple)
        tab_mod = MesMS.mapvalue(x -> x.mass, pLink.read_modification(joinpath(cfg, "modification.ini")))
        tab_xl = pLink.read_linker(joinpath(cfg, "xlink.ini")) |> NamedTuple
    end

    psm = pLink.read_psm_full(path)
    df_xl = psm.xl
    df_linear = psm.linear
    df_mono = psm.mono
    df_loop = psm.loop
    df_xl.linker .= linker

    process_crosslink!(df_xl, M, ε, ion_syms, ion_types, tab_ele, tab_aa, tab_mod, tab_xl)
    process_linear!(df_linear, M, ε, ion_syms, ion_types, tab_ele, tab_aa, tab_mod)

    MesMS.safe_save(p -> CSV.write(p, df_xl), joinpath(out, basename(path) * ".crosslink.XLCoverageReport.csv"))
    MesMS.safe_save(p -> CSV.write(p, df_linear), joinpath(out, basename(path) * ".linear.XLCoverageReport.csv"))

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

    main = replace(read(joinpath(DIR_DATA, "XLCoverageReport.html"), String),
        "{{ section }}" => "Overall Peptide Fragment Ion Coverage",
        "{{ id }}" => "all",
    )
    main *= map(ion_syms) do sym
        replace(read(joinpath(DIR_DATA, "XLCoverageReport.html"), String),
            "{{ section }}" => "Peptide Fragment Ion Coverage of $(sym) Ions",
            "{{ id }}" => "$(sym)",
        )
    end |> join

    script = replace(read(joinpath(DIR_DATA, "XLCoverageReport.js"), String), "{{ id }}" => "all")
    script *= map(ion_syms) do sym
        replace(read(joinpath(DIR_DATA, "XLCoverageReport.js"), String), "{{ id }}" => "$(sym)")
    end |> join

    html = replace(read(joinpath(DIR_DATA, "base.html"), String),
        "{{ title }}" => "TargetWizard Cross-link Fragment Ion Coverage Report",
        "{{ subtitle }}" => basename(path),
        "{{ main }}" => main,
        "{{ lib }}" => read(joinpath(DIR_DATA, "lib", "chartjs-4.2.1.js"), String),
        "{{ data }}" => data,
        "{{ script }}" => script,
    )
    path_out = joinpath(out, basename(path) * ".XLCoverageReport.html")
    MesMS.safe_save(io -> write(io, html), path_out)
    MesMS.open_url(path_out)
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="XLCoverageReport")
    ArgParse.@add_arg_table! settings begin
        "data"
            help = "pLink PSM path"
            required = true
        "--ms"
            help = "list of .mes or .ms2 files"
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
    paths = reduce(vcat, MesMS.match_path.(args["ms"], ".mes")) |> unique |> sort
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
