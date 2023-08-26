module TargetXView

using Sockets

import ArgParse
import CSV
import DataFrames
import MesMS
import MesUtil: pFind, pLink
import ProgressMeter: @showprogress
import RelocatableFolders: @path

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
    ys = map(n -> map(p -> MesMS.max_inten_ε(p, tg.mz + n * Δ / tg.z, ε), df_m1.peaks), -1:2)
    ls = [scatter(x=x, y=ys[2], mode="lines", name="M")]
    push!(ls, scatter(x=x, y=ys[3], mode="lines", name="M + 1 Da"))
    push!(ls, scatter(x=x, y=ys[4], mode="lines", name="M + 2 Da"))
    push!(ls, scatter(x=x, y=ys[1], mode="lines", name="M - 1 Da"))
    p1 = Plot(ls, Layout(; xaxis_title="retention time (s)", yaxis_title="abundance"))
    ls = [scatter(; get_rect(df_ft[i, :], ε/2)..., mode="lines", line_dash="dash", line_color="green", name="FT#$(i)") for i in tg.ft_]
    push!(ls, scatter(x=ones(2) * tg.start, y=[-ε, ε], mode="lines", line_dash="dash", line_color="red", name="RT start"))
    push!(ls, scatter(x=ones(2) * tg.stop, y=[-ε, ε], mode="lines", line_dash="dash", line_color="red", name="RT stop"))
    append!(ls, [scatter(x=df_m2.rt[i:i], y=[MesMS.error_rel(tg.mz, df_m2.mz[i])], name="MS2#$(df_m2.id[i])", customdata=[i], mode="markers", marker_size=4 * (1 + length(df_m2.psm[i]))) for i in tg.m2_all_])
    p2 = Plot(ls, Layout(; yaxis_title="m/z error"))
    p = [p1; p2; p_hit]
    relayout!(p, Layout(Subplots(rows=3, cols=1, row_heights=[2, 1, 2], vertical_spacing=0.02, shared_xaxes=true), clickmode="event+select"), height=600)
    return p
end

build_app(df_tg, df_xl, df_ft, df_m1, df_m2, df_psm, M2I, ele_plink, aa_plink, mod_plink, xl_plink, ele_pfind, aa_pfind, mod_pfind, ε) = begin
    df_tg_tb = DataFrames.select(df_tg, filter(n -> !endswith(n, "_"), names(df_tg)))
    app = dash(; assets_folder=DIR_DATA)
    app.layout = html_div() do
        html_h1("TargetXView", style=Dict("text-align"=>"center")),
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
        dash_datatable(
            id="xl_table",
            style_table=Dict("min-width"=>"100%", "overflow-x"=>"auto"),
            style_cell=Dict("overflow"=>"hidden", "min-width"=>"64px"),
            columns=[(; name=i, id=i) for i in names(df_xl)],
            filter_action="native",
            sort_action="native",
            sort_mode="multi",
            page_action="native",
            page_size=10,
            export_format="csv",
            export_headers="display",
        ),
        dcc_graph(id="lc_graph"),
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
        dcc_graph(id="seq_graph"),
        dcc_graph(id="psm_graph")
    end

    p_hit = plot_hit(df_m2)

    callback!(app,
        Output("xl_table", "data"),
        Output("lc_graph", "figure"),
        Input("tg_table", "derived_virtual_data"),
        Input("tg_table", "derived_virtual_selected_rows"),
    ) do v1, v2
        id = parse(Int, v1[v2[begin] + 1].id)
        tg = df_tg[id, :]
        table_data = Dict.(pairs.(eachrow(string.(df_xl[tg.xl_, :]))))
        fig = plot_lc(tg, df_ft, df_m1, df_m2, p_hit, ε)
        return table_data, fig
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
        if r.engine == :pLink
            seqs = (r.pep_a, r.pep_b)
            modss = (r.mod_a, r.mod_b)
            linker = xl_plink[Symbol(r.linker)]
            sites = (r.site_a, r.site_b)
            ionss = MesMS.Plot.build_ions_xl(m2.peaks, seqs, modss, linker, sites, ε, ele_plink, aa_plink, mod_plink)
            p_seq = MesMS.Plotly.seq_xl(seqs, modss, sites, ionss)
            p_psm = MesMS.Plotly.spec(m2.peaks, filter(i -> i.peak > 0, vcat(ionss...)))
        elseif r.engine == :pFind
            ions = MesMS.Plot.build_ions(m2.peaks, r.pep_a, r.mod_a, ε, ele_pfind, aa_pfind, mod_pfind)
            p_seq = MesMS.Plotly.seq(r.pep_a, r.mod_a, ions)
            p_psm = MesMS.Plotly.spec(m2.peaks, filter(i -> i.peak > 0, ions))
        end
        return p_seq, p_psm
    end
    return app
end

prepare(args) = begin
    path_ms = args["ms"]
    path_psm = args["psm"]
    out = mkpath(args["out"])
    path_xl = args["xl"]
    path_ft = args["ft"]
    path_psm_pf = args["psm_pf"]
    fmt = args["fmt"] |> Symbol
    linker = Symbol(args["linker"])
    ε = parse(Float64, args["error"]) * 1.0e-6
    fdr = parse(Float64, args["fdr"]) / 100
    cfg = args["cfg"]
    cfg_pf = args["cfg_pf"]
    host = parse(IPAddr, args["host"])
    port = parse(Int, args["port"])
    return (; path_ms, path_psm, out, path_xl, path_ft, path_psm_pf, fmt, linker, ε, fdr, cfg, cfg_pf, host, port)
end

process(path; path_ms, path_psm, out, path_xl, path_ft, path_psm_pf, fmt, linker, ε, fdr, cfg, cfg_pf, host, port) = begin
    M = MesMS.read_ms(path_ms)
    df_m1 = map(m -> (; m.id, rt=m.retention_time, m.peaks), M.MS1) |> DataFrames.DataFrame
    df_m2 = map(m -> (; m.id, mz=m.activation_center, rt=m.retention_time, m.peaks), M.MS2) |> DataFrames.DataFrame
    M2I = map(x -> x[2] => x[1], enumerate(df_m2.id)) |> Dict

    df_psm = pLink.read_psm_full(path_psm).xl
    df_psm = df_psm[df_psm.fdr .≤ fdr, :]
    df_psm.engine .= :pLink
    ns = [
        "Order", "Peptide", "Peptide_Type", "mh_calc", "Modifications", "Evalue", "Precursor_Mass_Error(Da)",
        "Proteins", "prot_type", "FileID", "LabelID", "Alpha_Matched", "Beta_Matched", "Alpha_Evalue", "Beta_Evalue",
        "Alpha_Seq_Coverage", "Beta_Seq_Coverage",
    ]
    DataFrames.select!(df_psm, DataFrames.Not(filter(x -> x ∈ names(df_psm), ns)))
    ns = [
        "engine", "mh", "mz", "z", "pep_a", "pep_b", "mod_a", "mod_b", "site_a", "site_b",
        "prot_a", "prot_b", "error", "title", "file", "scan", "idx_pre",
    ]
    DataFrames.select!(df_psm, ns, DataFrames.Not(ns))
    ("linker" ∉ names(df_psm)) && (df_psm.linker .= linker)

    if !isempty(path_psm_pf)
        df_psm_pf = pFind.read_psm(path_psm_pf)
        df_psm_pf.engine .= :pFind
        ns = [
            "Scan_No", "Sequence", "mh_calc", "Mass_Shift(Exp.-Calc.)", "score_raw", "Modification",
            "Specificity", "Positions", "Label", "Miss.Clv.Sites", "Avg.Frag.Mass.Shift", "Others", "mz_calc"
        ]
        DataFrames.select!(df_psm_pf, DataFrames.Not(filter(x -> x ∈ names(df_psm_pf), ns)))
        DataFrames.rename!(df_psm_pf, :pep => :pep_a, :mod => :mod_a, :proteins => :prot_a)
        df_psm = vcat(df_psm, df_psm_pf; cols=:union)
    end

    df_psm.id = Vector(1:size(df_psm, 1))
    DataFrames.select!(df_psm, :id, DataFrames.Not([:id]))
    df_psm.rt = [df_m2[M2I[r.scan], :rt] for r in eachrow(df_psm)]

    if isempty(cfg)
        ele_plink = pLink.read_element() |> NamedTuple
        aa_plink = map(x -> MesMS.mass(x, ele_plink), pLink.read_amino_acid() |> NamedTuple)
        mod_plink = MesMS.mapvalue(x -> x.mass, pLink.read_modification())
        xl_plink = pLink.read_linker() |> NamedTuple
    else
        ele_plink = pLink.read_element(joinpath(cfg, "element.ini")) |> NamedTuple
        aa_plink = map(x -> MesMS.mass(x, ele_plink), pLink.read_amino_acid(joinpath(cfg, "aa.ini")) |> NamedTuple)
        mod_plink = MesMS.mapvalue(x -> x.mass, pLink.read_modification(joinpath(cfg, "modification.ini")))
        xl_plink = pLink.read_linker(joinpath(cfg, "xlink.ini")) |> NamedTuple
    end

    if isempty(cfg_pf)
        ele_pfind = pFind.read_element() |> NamedTuple
        aa_pfind = map(x -> MesMS.mass(x, ele_pfind), pFind.read_amino_acid() |> NamedTuple)
        mod_pfind = MesMS.mapvalue(x -> x.mass, pFind.read_modification())
    else
        ele_pfind = pFind.read_element(joinpath(cfg_pf, "element.ini")) |> NamedTuple
        aa_pfind = map(x -> MesMS.mass(x, ele_pfind), pFind.read_amino_acid(joinpath(cfg_pf, "aa.ini")) |> NamedTuple)
        mod_pfind = MesMS.mapvalue(x -> x.mass, pFind.read_modification(joinpath(cfg_pf, "modification.ini")))
    end

    df_m2.psm = [df_psm[df_psm.scan .== r.id, :id] for r in eachrow(df_m2)]

    !isempty(path_xl) && @info "XL Candidtes loading from " * path_xl
    df_xl = (isempty(path_xl) ? [] : CSV.File(path_xl)) |> DataFrames.DataFrame
    df_xl.id = Vector(1:size(df_xl, 1))
    DataFrames.select!(df_xl, :id, DataFrames.Not([:id]))

    !isempty(path_ft) && @info "Feature loading from "* path_ft
    df_ft = (isempty(path_ft) ? [] : CSV.File(path_ft)) |> DataFrames.DataFrame
    df_ft.id = Vector(1:size(df_ft, 1))
    DataFrames.select!(df_ft, :id, DataFrames.Not([:id]))

    @info "Target loading from " * path
    df_tg = CSV.File(path) |> DataFrames.DataFrame
    df_tg.id = Vector(1:size(df_tg, 1))
    parse_target_list!(df_tg, fmt)
    DataFrames.select!(df_tg, [:id, :mz, :z, :start, :stop], DataFrames.Not([:id, :mz, :z, :start, :stop]))

    @info "XL Candidtes mapping"
    tmp = sort!([(x.mz::Float64, x.id::Int) for x in eachrow(df_xl)])
    mzs = map(x -> x[1], tmp)
    ids = map(x -> x[2], tmp)
    df_tg.xl_ = [sort(filter(x -> df_xl[x, :z] == r.z, ids[MesMS.argquery_ε(mzs, r.mz, ε)])) for r in eachrow(df_tg)]

    @info "Feature mapping"
    tmp = sort!([(x.mz::Float64, x.id::Int) for x in eachrow(df_ft)])
    mzs = map(x -> x[1], tmp)
    ids = map(x -> x[2], tmp)
    df_tg.ft_ = [sort(filter(x -> df_ft[x, :z] == r.z, ids[MesMS.argquery_ε(mzs, r.mz, ε)])) for r in eachrow(df_tg)]
    df_tg.n_ft = length.(df_tg.ft_)

    @info "MS2 mapping"
    tmp = sort!([(x.mz::Float64, x.id::Int) for x in eachrow(df_m2)])
    mzs = map(x -> x[1], tmp)
    ids = map(x -> x[2], tmp)
    df_tg.m2_all_ = [map(x -> M2I[x], sort(ids[MesMS.argquery_ε(mzs, r.mz, ε)])) for r in eachrow(df_tg)]
    df_tg.m2_ = [filter(i -> r.start ≤ df_m2.rt[i] ≤ r.stop, r.m2_all_) for r in eachrow(df_tg)]
    df_tg.n_m2_all = length.(df_tg.m2_all_)
    df_tg.n_m2 = length.(df_tg.m2_)
    df_tg.m2_all_id_ = [map(x -> df_m2.id[x], r.m2_all_) for r in eachrow(df_tg)]
    df_tg.m2_id_ = [map(x -> df_m2.id[x], r.m2_) for r in eachrow(df_tg)]

    @info "PSM mapping"
    df_tg.psm_all_ = [vcat(df_m2[r.m2_all_, :psm]...) for r in eachrow(df_tg)]
    df_tg.psm_ = [vcat(df_m2[r.m2_, :psm]...) for r in eachrow(df_tg)]
    df_tg.n_psm_all = length.(df_tg.psm_all_)
    df_tg.n_psm = length.(df_tg.psm_)

    filter_plink(x) = filter(i -> df_psm.engine[i] == :pLink, x)
    filter_pfind(x) = filter(i -> df_psm.engine[i] == :pFind, x)

    df_tg.n_psm_plink_all = length.(map(filter_plink, df_tg.psm_all_))
    df_tg.n_psm_plink = length.(map(filter_plink, df_tg.psm_))
    df_tg.n_psm_pfind_all = length.(map(filter_pfind, df_tg.psm_all_))
    df_tg.n_psm_pfind = length.(map(filter_pfind, df_tg.psm_))

    ns = filter(n -> !endswith(n, '_'), names(df_tg))
    DataFrames.select!(df_tg, ns, DataFrames.Not(ns))

    MesMS.safe_save(p -> CSV.write(p, df_tg), joinpath(out, "$(basename(splitext(path_ms)[1])).TargetXView.csv"))
    
    @async begin
        sleep(4)
        MesMS.open_url("http://$(host):$(port)")
    end
    app = build_app(df_tg, df_xl, df_ft, df_m1, df_m2, df_psm, M2I, ele_plink, aa_plink, mod_plink, xl_plink, ele_pfind, aa_pfind, mod_pfind, ε)
    run_server(app, host, port)
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="TargetXView")
    ArgParse.@add_arg_table! settings begin
        "target"
            help = "target list"
            required = true
        "--ms"
            help = ".mes or .ms1/2 file; .ms2/1 file should be in the same directory for .ms1/2"
            required = true
        "--psm"
            help = "pLink PSM path"
            required = true
        "--out"
            help = "output directory"
            default = "./out/"
            metavar = "./out/"
        "--xl"
            help = "candidate xl list"
            default = ""
        "--ft"
            help = "feature list"
            default = ""
        "--psm_pf"
            help = "pFind PSM path"
            default = ""
        "--fmt", "-f"
            help = "target list format: auto, TW, TmQE, TmFu"
            metavar = "auto|TW|TmQE|TmFu"
            default = "auto"
        "--linker"
            help = "default linker"
            default = "DSSO"
        "--error"
            help = "m/z error"
            metavar = "ppm"
            default = "20.0"
        "--fdr"
            help = "FDR threshold (%)"
            default = "Inf"
        "--cfg"
            help = "pLink config directory"
            default = ""
        "--cfg_pf"
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
