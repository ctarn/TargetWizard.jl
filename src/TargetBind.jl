module TargetBind

import ArgParse
import CSV
import DataFrames
import MesMS: MesMS, XL
import ProgressMeter: @showprogress

prepare(args) = begin
    @info "reading from " * args["target"]
    df = args["target"] |> CSV.File |> DataFrames.DataFrame
    out = mkpath(args["out"])
    mode = args["mode"] |> strip |> Symbol
    (mode ∉ [:center, :window, :extended_window]) && error("unknown mode: $(mode)")
    εt = parse(Float64, args["error_rt"])
    εm = parse(Float64, args["error_mz"]) * 1.0e-6
    fmts = split(args["fmt"], ",") .|> strip .|> Symbol
    return (; df, out, mode, εt, εm, fmts)
end

process(path; df, out, mode, εt, εm, fmts) = begin
    M = MesMS.read_ms(path; MS1=false).MS2
    name = basename(path) |> splitext |> first

    S = @showprogress map(M) do ms
        rt = ((df.start .- εt) .≤ ms.retention_time .≤ (df.stop .+ εt))
        if mode == :extended_window
            return ((ms.activation_center - ms.isolation_width / 2 - 1) .< df.mz .< (ms.activation_center + ms.isolation_width / 2)) .& rt
        elseif mode == :window
            return ((ms.activation_center - ms.isolation_width / 2) .< df.mz .< (ms.activation_center + ms.isolation_width / 2)) .& rt
        else # :center
            return MesMS.in_moe.(df.mz, ms.activation_center, εm) .& rt
        end
    end

    I = @showprogress map(S) do s
        map(r -> MesMS.Ion(r.mz, r.z), eachrow(df[s, :]))
    end
    for fmt in fmts
        ext = fmt ∈ [:csv, :tsv] ? "scan_precursor.$(fmt)" : fmt
        MesMS.safe_save(p -> MesMS.write_ms_with_precursor(p, M, I; fmt, name), joinpath(out, "$(name).$(ext)"))
    end
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="TargetBind")
    ArgParse.@add_arg_table! settings begin
        "data"
            help = "list of .mes or .ms2 files"
            nargs = '+'
            required = true
        "--target", "--tg"
            help = "target list file"
            metavar = "target"
            required = true
        "--out", "-o"
            help = "output directory"
            metavar = "output"
            default = "./out/"
        "--mode"
            help = "bind by isolation `center`, isolation `window`, or `extended_window`"
            metavar = "center|window|extended_window"
            default = "extended_window"
        "--error_rt", "--et"
            help = "retention time error"
            metavar = "sec"
            default = "16"
        "--error_mz", "--em"
            help = "MS1 mass error"
            metavar = "ppm"
            default = "10.0"
        "--fmt", "-f"
            help = "output format"
            metavar = "csv,tsv,ms2,mgf,pf2"
            default = "csv"
    end
    args = ArgParse.parse_args(settings)
    paths = reduce(vcat, MesMS.match_path.(args["data"], ".mes")) |> unique |> sort
    @info "file paths of selected MS data:"
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
