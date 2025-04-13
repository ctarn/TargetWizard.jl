module TargetBind

import ArgParse
import CSV
import DataFrames
import ProgressMeter: @showprogress
import UniMZ
import UniMZUtil: TMS

include("util.jl")

Δ = 1.0033

cosine_sim(xs, ys) = sum(abs.(xs .* ys)) / (sqrt(sum(xs .^ 2) * sum(ys .^2)))

getms2s(M, ion, i, Δrt, ε) = begin
    ms2s = [i]
    id = i
    n = 2
    while id > 1 && n > 0
        id = id - 1
        m = M[id]
        abs(M[i].retention_time - m.retention_time) > Δrt && break
        !UniMZ.in_moe(m.activation_center, ion.mz, ε) && continue
        push!(ms2s, id)
        n -= 1
    end
    id = i
    n = 2
    while id < length(M) && n > 0
        id = id + 1
        m = M[id]
        abs(M[i].retention_time - m.retention_time) > Δrt && break
        !UniMZ.in_moe(m.activation_center, ion.mz, ε) && continue
        push!(ms2s, id)
        n -= 1
    end
    return M[sort!(ms2s)]
end

prepare(args) = begin
    @info "reading from " * args["target"]
    df = args["target"] |> CSV.File |> DataFrames.DataFrame
    out = mkpath(args["out"])
    mode = args["mode"] |> strip |> Symbol
    (mode ∉ [:center, :window, :extended_window]) && error("unknown mode: $(mode)")
    εt = parse(Float64, args["error_rt"])
    εm = parse(Float64, args["error_mz"]) * 1.0e-6
    fmt_target = Symbol(args["fmt_target"])
    @info "specified target format: $(fmt_target)"
    fmts = split(args["fmt"], ",") .|> strip .|> Symbol
    return (; df, out, mode, εt, εm, fmt_target, fmts)
end

process(path; df, out, mode, εt, εm, fmt_target, fmts) = begin
    M = UniMZ.read_ms(path)
    M1 = M.MS1 |> UniMZ.dict_by_id
    M2 = M.MS2
    name = basename(path) |> splitext |> first
    df = TMS.parse_target_list!(copy(df), fmt_target)
    @info "MS2 matching..."
    I = @showprogress map(M2) do ms
        if mode == :extended_window
            s = ((ms.activation_center - ms.isolation_width / 2 - 1) .< df.mz .< (ms.activation_center + ms.isolation_width / 2))
        elseif mode == :window
            s = ((ms.activation_center - ms.isolation_width / 2) .< df.mz .< (ms.activation_center + ms.isolation_width / 2))
        else # :center
            s = UniMZ.in_moe.(df.mz, ms.activation_center, εm)
        end
        s = s .& ((df.start .- εt) .≤ ms.retention_time .≤ (df.stop .+ εt))
        return map(r -> UniMZ.Ion(r.mz, r.z), eachrow(df[s, :]))
    end
    if false
        @info "MS2 preprocessing using isotopic pattern..."
        M2 = @showprogress map(M2) do ms
            ps = map(ms.peaks) do p
                if any(z -> UniMZ.query_ε(ms.peaks, p.mz - Δ / z, εm) |> !isempty, 1:3)
                    return UniMZ.Peak[]
                end
                zs = filter(z -> UniMZ.query_ε(ms.peaks, p.mz + Δ / z, εm) |> !isempty, 1:3)
                zs = isempty(zs) ? [1] : zs
                return map(z -> UniMZ.Peak(UniMZ.mz_to_mh(p.mz, z), p.inten), zs)
            end
            return UniMZ.fork(ms; peaks=reduce(vcat, ps)|>sort)
        end
    end
    @info "MS2 preprocessing using xic pattern..."
    S = @showprogress map(eachindex(M2), M2, I) do idx_ms, ms, ions
        ss = map(ions) do ion
            ms2s = getms2s(M2, ion, idx_ms, 30, εm)
            ms1s = map(ms2 ->M1[ms2.pre], ms2s)
            xic1 = map(ms1 -> UniMZ.max_inten_ε(ms1.peaks, ion.mz, εm), ms1s)
            return map(ms.peaks) do p
                xic2 = map(ms2 -> UniMZ.max_inten_ε(ms2.peaks, p.mz, εm), ms2s)
                return (1 - acos(cosine_sim(xic1, xic2)) / π * 2) ≥ 0.5
            end
        end
        return isempty(ions) ? trues(length(ms.peaks)) : reduce(.|, ss)
    end
    M2 = @showprogress map(M2, S) do ms, s
        UniMZ.fork(ms; peaks=ms.peaks[s])
    end
    for fmt in fmts
        ext = fmt ∈ [:csv, :tsv] ? "scan_precursor.$(fmt)" : fmt
        UniMZ.safe_save(p -> UniMZ.write_ms_with_precursor(p, M2, I; fmt, name), joinpath(out, "$(name).$(ext)"))
    end
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="TargetBind")
    ArgParse.@add_arg_table! settings begin
        "data"
            help = "list of .umz or .ms2 files"
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
        "--fmt_target", "--ft"
            help = "target list format: auto, TW, TmQE, TmFu"
            metavar = "auto|TW|TmQE|TmFu"
            default = "auto"
        "--fmt", "-f"
            help = "output format"
            metavar = "csv,tsv,ms2,mgf,pf2"
            default = "csv"
    end
    args = ArgParse.parse_args(settings)
    paths = reduce(vcat, UniMZ.match_path.(args["data"], ".umz")) |> unique |> sort
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
