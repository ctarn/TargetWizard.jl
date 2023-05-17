module AquisitionReport

import ArgParse
import CSV
import DataFrames
import MesMS
import MesUtil: pFind, pLink
import RelocatableFolders: @path

const DIR_DATA = @path joinpath(@__DIR__, "../data")

prepare(args) = begin
    out = mkpath(args["out"])
    return (; out)
end

plot(path; out) = begin
    fname_m1 = splitext(path)[1] * ".ms1"
    @info "MS1 loading from " * fname_m1
    M1 = MesMS.read_ms1(fname_m1)
    fname_m2 = splitext(path)[1] * ".ms2"
    @info "MS2 loading from " * fname_m2
    M2 = MesMS.read_ms2(fname_m2)
    name = splitext(basename(path))[1]

    rt1 = map(m -> m.retention_time, M1)
    rt2 = map(m -> m.retention_time, M2)
    rt_max = ceil(Int, max(maximum(rt1; init=0.0), maximum(rt2; init=0.0)))
    mz2 = map(m -> m.activation_center, M2)
    it1 = map(m -> m.injection_time, M1)
    it2 = map(m -> m.injection_time, M2)
    tic1 = log10.(map(m -> m.total_ion_current, M1) .+ 1)
    tic2 = log10.(map(m -> m.total_ion_current, M2) .+ 1)
    bpi1 = log10.(map(m -> m.base_peak_intensity, M1) .+ 1)
    bpi2 = log10.(map(m -> m.base_peak_intensity, M2) .+ 1)
    bpm1 = map(m -> m.base_peak_mass, M1)
    bpm2 = map(m -> m.base_peak_mass, M2)

    data = """
rt1 = [$(join(string.(rt1), ","))]
rt2 = [$(join(string.(rt2), ","))]
rt_max = $(rt_max)
f1_mean = $(length(M1) / rt_max)
f2_mean = $(length(M2) / rt_max)
mz2 = [$(join(string.(mz2), ","))]
mz2_min = $(minimum(mz2))
mz2_max = $(maximum(mz2))
it1 = [$(join(string.(it1), ","))]
it2 = [$(join(string.(it2), ","))]
it_max = $(max(maximum(it1), maximum(it2)))
tic1 = [$(join(string.(tic1), ","))]
tic2 = [$(join(string.(tic2), ","))]
tic_min = $(min(minimum(tic1), minimum(tic2)))
tic_max = $(max(maximum(tic1), maximum(tic2)))
bpi1 = [$(join(string.(bpi1), ","))]
bpi2 = [$(join(string.(bpi2), ","))]
bpi_min = $(min(minimum(bpi1), minimum(bpi2)))
bpi_max = $(max(maximum(bpi1), maximum(bpi2)))
bpm1 = [$(join(string.(bpm1), ","))]
bpm2 = [$(join(string.(bpm2), ","))]
bpm_min = $(min(minimum(bpm1), minimum(bpm2)))
bpm_max = $(max(maximum(bpm1), maximum(bpm2)))
"""

    html = replace(read(joinpath(DIR_DATA, "base.html"), String),
        "{{ title }}" => "TargetWizard Data Aquisition Report",
        "{{ subtitle }}" => name,
        "{{ main }}" => read(joinpath(DIR_DATA, "AquisitionReport.html"), String),
        "{{ lib }}" => read(joinpath(DIR_DATA, "lib", "chartjs-4.2.1.js"), String),
        "{{ data }}" => data,
        "{{ script }}" => read(joinpath(DIR_DATA, "AquisitionReport.js"), String),
    )
    path_out = joinpath(out, name * ".AquisitionReport.html")
    MesMS.safe_save(io -> write(io, html), path_out)
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="AquisitionReport")
    ArgParse.@add_arg_table! settings begin
        "--out", "-o"
            help = "output directory"
            metavar = "output"
            default = "./out/"
        "data"
            help = "list of .MS1 files"
            nargs = '+'
            required = true
    end
    args = ArgParse.parse_args(settings)
    paths = (sort∘unique∘reduce)(vcat, MesMS.match_path.(args["data"], ".ms1"); init=String[])
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
