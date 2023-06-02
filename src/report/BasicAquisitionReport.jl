module BasicAquisitionReport

import ArgParse
import CSV
import DataFrames
import MesMS
import MesUtil: pFind, pLink
import RelocatableFolders: @path

const DIR_DATA = @path joinpath(@__DIR__, "../../data")

prepare(args) = begin
    out = mkpath(args["out"])
    return (; out)
end

plot(path; out) = begin
    name = splitext(basename(path))[1]
    M = MesMS.read_ms(path)
    rt1 = map(m -> m.retention_time, M.MS1)
    rt2 = map(m -> m.retention_time, M.MS2)
    rt_max = ceil(Int, max(maximum(rt1; init=0.0), maximum(rt2; init=0.0)))
    mz2 = map(m -> m.activation_center, M.MS2)
    it1 = map(m -> m.injection_time, M.MS1)
    it2 = map(m -> m.injection_time, M.MS2)
    tic1 = log10.(map(m -> m.total_ion_current, M.MS1) .+ 1)
    tic2 = log10.(map(m -> m.total_ion_current, M.MS2) .+ 1)
    bpi1 = log10.(map(m -> m.base_peak_intensity, M.MS1) .+ 1)
    bpi2 = log10.(map(m -> m.base_peak_intensity, M.MS2) .+ 1)
    bpm1 = map(m -> m.base_peak_mass, M.MS1)
    bpm2 = map(m -> m.base_peak_mass, M.MS2)

    data = """
const RT1 = [$(join(string.(rt1), ","))]
const RT2 = [$(join(string.(rt2), ","))]
const RT_MAX = $(rt_max)
const F1_MEAN = $(length(M.MS1) / rt_max)
const F2_MEAN = $(length(M.MS2) / rt_max)
const MZ2 = [$(join(string.(mz2), ","))]
const MZ2_MIN = $(minimum(mz2))
const MZ2_MAX = $(maximum(mz2))
const IT1 = [$(join(string.(it1), ","))]
const IT2 = [$(join(string.(it2), ","))]
const IT_MAX = $(max(maximum(it1), maximum(it2)))
const TIC1 = [$(join(string.(tic1), ","))]
const TIC2 = [$(join(string.(tic2), ","))]
const TIC_MIN = $(min(minimum(tic1), minimum(tic2)))
const TIC_MAX = $(max(maximum(tic1), maximum(tic2)))
const BPI1 = [$(join(string.(bpi1), ","))]
const BPI2 = [$(join(string.(bpi2), ","))]
const BPI_MIN = $(min(minimum(bpi1), minimum(bpi2)))
const BPI_MAX = $(max(maximum(bpi1), maximum(bpi2)))
const BPM1 = [$(join(string.(bpm1), ","))]
const BPM2 = [$(join(string.(bpm2), ","))]
const BPM_MIN = $(min(minimum(bpm1), minimum(bpm2)))
const BPM_MAX = $(max(maximum(bpm1), maximum(bpm2)))
"""

    html = replace(read(joinpath(DIR_DATA, "base.html"), String),
        "{{ title }}" => "TargetWizard Data Aquisition Report",
        "{{ subtitle }}" => name,
        "{{ main }}" => read(joinpath(DIR_DATA, "BasicAquisitionReport.html"), String),
        "{{ lib }}" => read(joinpath(DIR_DATA, "lib", "chartjs-4.2.1.js"), String),
        "{{ data }}" => data,
        "{{ script }}" => read(joinpath(DIR_DATA, "BasicAquisitionReport.js"), String),
    )
    path_out = joinpath(out, name * ".BasicAquisitionReport.html")
    MesMS.safe_save(io -> write(io, html), path_out)
    MesMS.open_url(path_out)
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="BasicAquisitionReport")
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
    paths = (sort∘unique∘reduce)(vcat, MesMS.match_path.(args["data"], ".mes"); init=String[])
    @info "file paths of selected MS data:"
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
