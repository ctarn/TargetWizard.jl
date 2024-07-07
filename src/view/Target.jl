"""
Target View
"""
module TargetView

"""
- target list (multiple formats supported), e.g., `TargetWizard.all.TW.target.csv`.
- traditional and targeted mass spectrometry data, e.g., `DDA.raw` and `TMS.raw`.
  - The raw data should be converted into an open-source format such as MS1/MS2. [ThermoRawRead](http://thermorawread.ctarn.io) is recommended.
- (filtered) identification results of targeted mass spectrometry data, e.g., `TMS_fdr.pfind.csv`.
- optional: precursor list detected by [`PepPre`](http://peppre.ctarn.io)
"""
require = true

"""
Once finished, TargetWizard will start a web server for user to review each target.
([example](../assets/manual/TargetView.pdf))
"""
output = true

"""
# Max. MS1 Mass Error
mass error used to match targets, PSMs, and MS scans.

# FDR Threshold
used to filter PSM list.

# MS Sim. Thres.
used to match traditional and targeted MS scans.

![Target View](../assets/manual/TargetView.png)
"""
usage = true

"""
# Select Target
![Target List](../assets/manual/TargetView_target.png)

# View MS Acquisition
![MS Acquisition](../assets/manual/TargetView_ms.png)

# Select PSM
![PSM List](../assets/manual/TargetView_psm_list.png)

# View PSM
![PSM](../assets/manual/TargetView_psm.png)
"""
example = true

using Sockets
using Statistics

import ArgParse
import CSV
import DataFrames
import ProgressMeter: @showprogress
import RelocatableFolders: @path
import UniMZ
import UniMZUtil: Proteomics, TMS, pFind

using Dash
using PlotlyBase

include("../util.jl")

const DIR_DATA = @path joinpath(@__DIR__, "../../data/dash")

Δ = 1.00335

plot_hit(df_m2) = begin
    ls = [scatter(x=df_m2.rt, y=df_m2.mz, customdata=1:size(df_m2, 1), mode="markers", name="All MS2", line_color="grey")]
    return Plot(ls, Layout(; xaxis_title="retention time (s)", yaxis_title="m/z"))
end

get_rect(ft, y) = (x=[ft.rtime_start, ft.rtime_start, ft.rtime_stop, ft.rtime_stop, ft.rtime_start], y=[y, -y, -y, y, y])

plot_lc(tg, df_ft, df_m1, df_m2, p_hit, ε) = begin
    x = df_m1.rt
    ys = map(n -> map(p -> UniMZ.max_inten_ε(p, tg.mz + n * Δ / tg.z, ε), df_m1.peaks), -1:2)
    ls = [scatter(x=x, y=ys[2], mode="lines", name="M")]
    push!(ls, scatter(x=x, y=ys[3], mode="lines", name="M + 1 Da"))
    push!(ls, scatter(x=x, y=ys[4], mode="lines", name="M + 2 Da"))
    push!(ls, scatter(x=x, y=ys[1], mode="lines", name="M - 1 Da"))
    p1 = Plot(ls, Layout(; xaxis_title="retention time (s)", yaxis_title="abundance"))
    ls = [scatter(; get_rect(df_ft[i, :], ε/2)..., mode="lines", line_dash="dash", line_color="green", name="FT#$(i)") for i in tg.ft_]
    push!(ls, scatter(x=ones(2) * tg.start, y=[-ε, ε], mode="lines", line_dash="dash", line_color="red", name="RT start"))
    push!(ls, scatter(x=ones(2) * tg.stop, y=[-ε, ε], mode="lines", line_dash="dash", line_color="red", name="RT stop"))
    append!(ls, [scatter(x=df_m2.rt[i:i], y=[UniMZ.error_rel(tg.mz, df_m2.mz[i])], name="MS2#$(df_m2.id[i])", customdata=[i], mode="markers", marker_size=4 * (1 + length(df_m2.psm[i]))) for i in tg.m2_all_])
    p2 = Plot(ls, Layout(; yaxis_title="m/z error"))
    p = [p1; p2; p_hit]
    relayout!(p, Layout(Subplots(rows=3, cols=1, row_heights=[2, 1, 2], vertical_spacing=0.02, shared_xaxes=true), clickmode="event+select"), height=600)
    return p
end

build_app(df_tg, df_ft, df_m1, df_m2, df_psm, M2I, ele_pfind, aa_pfind, mod_pfind, ε) = begin
    df_tg_tb = DataFrames.select(df_tg, filter(n -> !endswith(n, "_"), names(df_tg)))
    app = dash(; assets_folder=DIR_DATA)
    app.layout = html_div() do
        html_h1("TargetView", style=Dict("text-align"=>"center")),
        dash_datatable(
            id="tg_table",
            style_table=Dict("min-width"=>"100%", "overflow-x"=>"auto"),
            style_cell=Dict("overflow"=>"hidden", "text-overflow"=>"ellipsis", "min-width"=>"64px", "max-width"=>"256px"),
            columns=[(; name=i, id=i) for i in names(df_tg_tb)],
            data=Dict.(pairs.(eachrow(string.(df_tg_tb)))),
            filter_action="native",
            sort_action="native",
            sort_mode="multi",
            row_selectable="single",
            page_action="native",
            page_size=10,
            export_format="csv",
            export_headers="display",
        ),
        dcc_graph(id="lc_graph", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"TargetView_LC"))),
        dash_datatable(
            id="psm_table",
            style_table=Dict("min-width"=>"100%", "overflow-x"=>"auto"),
            style_cell=Dict("overflow"=>"hidden", "min-width"=>"64px",),
            columns=[(; name=i, id=i) for i in names(df_psm)],
            filter_action="native",
            sort_action="native",
            sort_mode="multi",
            row_selectable="single",
            page_action="native",
            page_size=10,
            export_format="csv",
            export_headers="display",
        ),
        dcc_graph(id="seq_graph", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"TargetView_SEQ"))),
        dcc_graph(id="psm_graph", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"TargetView_PSM")))
    end

    p_hit = plot_hit(df_m2)

    callback!(app,
        Output("lc_graph", "figure"),
        Input("tg_table", "derived_virtual_data"),
        Input("tg_table", "derived_virtual_selected_rows"),
    ) do v1, v2
        id = parse(Int, v1[v2[begin] + 1].id)
        tg = df_tg[id, :]
        return plot_lc(tg, df_ft, df_m1, df_m2, p_hit, ε)
    end

    callback!(app,
        Output("psm_table", "data"),
        Input("lc_graph", "selectedData"),
    ) do v
        ids = map(p -> p.customdata, v.points)
        psms = df_psm[vcat(df_m2[ids, :psm]...), :]
        return Dict.(pairs.(eachrow(string.(psms))))
    end

    callback!(app,
        Output("seq_graph", "figure"),
        Output("psm_graph", "figure"),
        Input("psm_table", "derived_virtual_data"),
        Input("psm_table", "derived_virtual_selected_rows"),
    ) do v1, v2
        id = parse(Int, v1[v2[begin] + 1].id)
        r = df_psm[id, :]
        m2 = df_m2[M2I[r.scan], :]
        ions = UniMZ.build_ions(m2.peaks, r.pep, r.mod, ε, ele_pfind, aa_pfind, mod_pfind)
        p_seq = UniMZ.Plotly.seq(r.pep, r.mod, ions)
        p_psm = UniMZ.Plotly.spec(m2.peaks, filter(i -> i.peak > 0, ions))
        return p_seq, p_psm
    end
    return app
end

prepare(args) = begin
    path_ms = args["ms"]
    @info "file path of selected data:"
    println("\t$(path_ms)")
    paths_ms_old = reduce(vcat, UniMZ.match_path.(args["ms_old"], ".umz")) |> unique |> sort
    @info "file paths of selected original data:"
    foreach(x -> println("$(x[1]):\t$(x[2])"), enumerate(paths_ms_old))
    path_psm = args["psm"]
    out = mkpath(args["out"])
    path_ft = args["ft"]
    fmt = args["fmt"] |> Symbol
    ε = parse(Float64, args["error"]) * 1.0e-6
    fdr = parse(Float64, args["fdr"]) / 100
    decoy = args["decoy"]::Bool
    τ_ms_sim = parse(Float64, args["ms_sim_thres"])
    tab_ele, tab_aa, tab_mod = pFind.read_mass_table(args["cfg"])
    host = parse(IPAddr, args["host"])
    port = parse(Int, args["port"])
    return (; path_ms, paths_ms_old, path_psm, out, path_ft, fmt, ε, fdr, decoy, τ_ms_sim, tab_ele, tab_aa, tab_mod, host, port)
end

process(path; path_ms, paths_ms_old, path_psm, out, path_ft, fmt, ε, fdr, decoy, τ_ms_sim, tab_ele, tab_aa, tab_mod, host, port) = begin
    M = UniMZ.read_ms(path_ms)
    df_m1 = map(m -> (; m.id, rt=m.retention_time, m.peaks), M.MS1) |> DataFrames.DataFrame
    df_m2 = map(m -> (; m.id, mz=m.activation_center, rt=m.retention_time, m.peaks), M.MS2) |> DataFrames.DataFrame
    M2I = map(x -> x[2] => x[1], enumerate(df_m2.id)) |> Dict

    M_old = map(p -> splitext(basename(p))[1] => UniMZ.dict_by_id(UniMZ.read_ms(p).MS2), paths_ms_old) |> Dict

    df = pFind.read_psm(path_psm)
    df.engine .= :pFind
    filter!(r -> r.fdr .≤ fdr, df)
    !decoy && filter!(r -> r.td == :T, df)

    ns = [
        "Order", "Peptide", "Peptide_Type", "mh_calc", "Modifications", "Evalue", "Precursor_Mass_Error(Da)",
        "Proteins", "prot_type", "FileID", "LabelID", "Alpha_Matched", "Beta_Matched", "Alpha_Evalue", "Beta_Evalue",
        "Alpha_Seq_Coverage", "Beta_Seq_Coverage",
    ]
    DataFrames.select!(df, DataFrames.Not(filter(x -> x ∈ names(df), ns)))
    ns = ["engine", "mh", "mz", "z", "pep", "mod", "prot", "title", "file", "scan", "idx_pre"]
    DataFrames.select!(df, ns, DataFrames.Not(ns))

    ion_syms = ["b", "y"]
    ion_types = map(i -> getfield(UniMZ, Symbol("ion_$(i)")), ion_syms)
    M_ = [splitext(basename(path_ms))[1] => UniMZ.dict_by_id(M.MS2)] |> Dict
    Proteomics.calc_cov!(df, M_, ε, ion_syms, ion_types, tab_ele, tab_aa, tab_mod)

    df.credible = map(r -> r.cov_ion_y ≥ 0.6 && r.cov_ion_b ≥ 0.4, eachrow(df))

    ns = [
        "Scan_No", "Sequence", "mh_calc", "Mass_Shift(Exp.-Calc.)", "score_raw", "Modification",
        "Specificity", "Positions", "Label", "Miss.Clv.Sites", "Avg.Frag.Mass.Shift", "Others", "mz_calc"
    ]
    DataFrames.select!(df, DataFrames.Not(filter(x -> x ∈ names(df), ns)))

    df.id = Vector(1:size(df, 1))
    DataFrames.select!(df, :id, DataFrames.Not([:id]))
    df.rt = [df_m2[M2I[r.scan], :rt] for r in eachrow(df)]

    df_m2.psm = [df[df.scan .== r.id, :id] for r in eachrow(df_m2)]

    !isempty(path_ft) && @info "Feature loading from "* path_ft
    df_ft = (isempty(path_ft) ? [] : CSV.File(path_ft)) |> DataFrames.DataFrame
    df_ft.id = Vector(1:size(df_ft, 1))
    DataFrames.select!(df_ft, :id, DataFrames.Not([:id]))

    @info "Target loading from " * path
    df_tg = CSV.File(path) |> DataFrames.DataFrame
    df_tg.id = Vector(1:size(df_tg, 1))
    TMS.parse_target_list!(df_tg, fmt)
    DataFrames.select!(df_tg, [:id, :mz, :z, :start, :stop], DataFrames.Not([:id, :mz, :z, :start, :stop]))
    "mod_a" ∈ names(df_tg) && (df_tg.mod_a = parse.(Array{UniMZ.Mod}, unify_mods_str.(df_tg.mod_a)))
    "mod_b" ∈ names(df_tg) && (df_tg.mod_b = parse.(Array{UniMZ.Mod}, unify_mods_str.(df_tg.mod_b)))

    @info "Feature mapping"
    tmp = sort!([(x.mz::Float64, x.id::Int) for x in eachrow(df_ft)])
    mzs = map(x -> x[1], tmp)
    ids = map(x -> x[2], tmp)
    df_tg.ft_ = [sort(filter(x -> df_ft[x, :z] == r.z, ids[UniMZ.argquery_ε(mzs, r.mz, ε)])) for r in eachrow(df_tg)]
    df_tg.n_ft = length.(df_tg.ft_)

    Ks = ["", "_allsim", "_all"]

    calc_sim(dda, tda) = map(p -> !isempty(UniMZ.argquery_ε(tda.peaks, p.mz, ε)), M_old[dda.file][dda.scan].peaks) |> mean

    @info "MS2 mapping"
    tmp = sort!([(x.mz::Float64, x.id::Int) for x in eachrow(df_m2)])
    mzs = map(x -> x[1], tmp)
    ids = map(x -> x[2], tmp)
    df_tg.m2_all_ = [map(x -> M2I[x], sort(ids[UniMZ.argquery_ε(mzs, r.mz, ε)])) for r in eachrow(df_tg)]
    df_tg.m2_allsim_ = [filter(i -> calc_sim(r, df_m2[i, :]) ≥ τ_ms_sim, r.m2_all_) for r in eachrow(df_tg)]
    df_tg.m2_ = [filter(i -> r.start ≤ df_m2.rt[i] ≤ r.stop, r.m2_all_) for r in eachrow(df_tg)]
    for K in Ks
        df_tg[!, "n_m2$(K)"] = length.(df_tg[!, "m2$(K)_"])
    end
    for K in Ks
        df_tg[!, "m2$(K)_id_"] = [map(x -> df_m2.id[x], r["m2$(K)_"]) for r in eachrow(df_tg)]
    end

    @info "MS2 similarity calculating"
    for K in Ks
        df_tg[!, "m2$(K)_sim_"] = @showprogress map(eachrow(df_tg)) do r
            map(s -> "$(s.id):$(round(calc_sim(r, s); digits=2))", eachrow(df_m2[r["m2$(K)_"], :])) |> xs -> join(xs, ";")
        end
    end

    @info "PSM mapping"
    for K in Ks
        df_tg[!, "psm$(K)_"] = [vcat(df_m2[r["m2$(K)_"], :psm]...) for r in eachrow(df_tg)]
        df_tg[!, "n_psm$(K)"] = length.(df_tg[!, "psm$(K)_"])
    end

    ns = filter(n -> !endswith(n, '_'), names(df_tg))
    DataFrames.select!(df_tg, ns, DataFrames.Not(ns))

    @async begin
        sleep(4)
        UniMZ.open_url("http://$(host):$(port)")
    end
    app = build_app(df_tg, df_ft, df_m1, df_m2, df, M2I, tab_ele, tab_aa, tab_mod, ε)
    run_server(app, host, port)
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="TargetView")
    ArgParse.@add_arg_table! settings begin
        "target"
            help = "target list"
            required = true
        "--ms"
            help = ".umz or .ms1/2 file; .ms2/1 file should be in the same directory for .ms1/2"
            required = true
        "--ms_old"
            help = "origianl .umz or .ms1/2 files; .ms2/1 files should be in the same directory for .ms1/2"
            nargs = '+'
            required = true
        "--psm"
            help = "pFind PSM path"
            required = true
        "--out"
            help = "output directory"
            default = "./out/"
            metavar = "./out/"
        "--ft"
            help = "feature list"
            default = ""
        "--fmt", "-f"
            help = "target list format: auto, TW, TmQE, TmFu"
            metavar = "auto|TW|TmQE|TmFu"
            default = "auto"
        "--error"
            help = "m/z error"
            metavar = "ppm"
            default = "20.0"
        "--fdr"
            help = "FDR threshold (%)"
            default = "Inf"
        "--decoy"
            help = "preserve decoy identifications"
            action = :store_true
        "--ms_sim_thres"
            help = "threshold of MS similarity"
            default = "0.5"
        "--cfg"
            help = "pFind config directory"
            default = ""
        "--host"
            help = "hostname"
            metavar = "hostname"
            default = "127.0.0.1"
        "--port"
            help = "port"
            metavar = "port"
            default = "30030"
    end
    args = ArgParse.parse_args(settings)
    process(args["target"]; prepare(args)...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

julia_main()::Cint = begin
    main()
    return 0
end

end
