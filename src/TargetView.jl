module TargetView

using Sockets

import ArgParse
import CSV
import DataFrames
import MesMS
import MesUtil: pFind, pLink
import ProgressMeter: @showprogress

using Dash
using PlotlyBase

Δ = 1.00335

plot_hit(df_m2) = begin
    ls = [scatter(x=df_m2.rt, y=df_m2.mz, customdata=1:DataFrames.nrow(df_m2), mode="markers", name="All MS2", line_color="grey")]
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
    app = dash()
    app.layout = html_div() do
        html_h1("TargetView", style=Dict("text-align"=>"center")),
        dash_datatable(
            id="tg_table",
            style_table=Dict("min-width"=>"100%", "overflow-x"=>"auto"),
            style_cell=Dict("overflow"=>"hidden", "text-overflow"=>"ellipsis", "min-width"=>"64px", "max-width"=>"256px"),
            columns=[(; name=i, id=i) for i in names(df_tg_tb)],
            data=Dict.(pairs.(eachrow(df_tg_tb))),
            filter_action="native",
            sort_action="native",
            sort_mode="multi",
            row_selectable="single",
            page_action="native",
            page_size=10,
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
        id = v1[v2[begin] + 1].id
        tg = df_tg[id, :]
        table_data = Dict.(pairs.(eachrow(df_xl[tg.xl_, :])))
        fig = plot_lc(tg, df_ft, df_m1, df_m2, p_hit, ε)
        return table_data, fig
    end

    callback!(app,
        Output("psm_table", "data"),
        Input("lc_graph", "selectedData"),
    ) do v
        ids = map(p -> p.customdata, v.points)
        psms = df_psm[vcat(df_m2[ids, :psm]...), :]
        return Dict.(pairs.(eachrow(psms)))
    end

    callback!(app,
        Output("seq_graph", "figure"),
        Output("psm_graph", "figure"),
        Input("psm_table", "derived_virtual_data"),
        Input("psm_table", "derived_virtual_selected_rows"),
    ) do v1, v2
        id = v1[v2[begin] + 1].id
        r = df_psm[id, :]
        m2 = df_m2[M2I[r.scan], :]
        if r.engine == :pLink
            seqs = (r.pep_a, r.pep_b)
            modss = (r.mod_a, r.mod_b)
            linker = xl_plink[Symbol(r.linker)]
            sites = (r.site_a, r.site_b)
            ionss = MesMS.Plot.build_ions_xl(m2.peaks, seqs, modss, linker, sites, ε, ele_plink, aa_plink, mod_plink)
            ionss = map(ionss) do ions
                map(ions) do ion
                    a, b, c = match(r"\$(.+)_\{(.+)\}\^\{(.+)\}\$", ion.text).captures
                    return (; ion..., text="$(a)($(b))$(c)")
                end
            end
            p_seq = MesMS.Plotly.seq_xl(seqs, modss, sites, ionss)
            p_psm = MesMS.Plotly.spec(m2.peaks, filter(i -> i.peak > 0, vcat(ionss...)))
        elseif r.engine == :pFind
            ions = MesMS.Plot.build_ions(m2.peaks, r.pep_a, r.mod_a, ε, ele_pfind, aa_pfind, mod_pfind)
            ions = map(ions) do ion
                a, b, c = match(r"\$(.+)_\{(.+)\}\^\{(.+)\}\$", ion.text).captures
                return (; ion..., text="$(a)($(b))$(c)")
            end
            p_seq = MesMS.Plotly.seq(r.pep_a, r.mod_a, ions)
            p_psm = MesMS.Plotly.spec(m2.peaks, filter(i -> i.peak > 0, ions))
        end
        return p_seq, p_psm
    end
    return app
end

prepare(args) = begin
    linker = Symbol(args["linker"])
    host = parse(IPAddr, args["host"])
    port = parse(Int, args["port"])
    path_xl = args["xl"]
    path_ft = args["ft"]
    path_tg = args["tg"]
    path_psm = args["psm"]
    path_psm_pf = args["psm_pf"]
    fdr = parse(Float64, args["fdr"]) / 100
    cfg = args["cfg"]
    cfg_pf = args["cfg_pf"]
    ε = parse(Float64, args["error"]) * 1.0e-6
    return (; linker, host, port, path_xl, path_ft, path_tg, path_psm, fdr, path_psm_pf, cfg, cfg_pf, ε)
end

report(path; linker, host, port, path_xl, path_ft, path_tg, path_psm, fdr, path_psm_pf, cfg, cfg_pf, ε) = begin
    path_ms1 = splitext(path)[1] * ".ms1"
    path_ms2 = splitext(path)[1] * ".ms2"
    df_m1 = map(MesMS.read_ms1(path_ms1)) do m
        (; m.id, rt=m.retention_time, m.peaks)
    end |> DataFrames.DataFrame
    df_m2 = map(MesMS.read_ms2(path_ms2)) do m
        (; m.id, mz=m.activation_center, rt=m.retention_time, m.peaks)
    end |> DataFrames.DataFrame

    M1I = map(m -> m[2].id => m[1], enumerate(eachrow(df_m1))) |> Dict
    M2I = map(m -> m[2].id => m[1], enumerate(eachrow(df_m2))) |> Dict

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
        "prot_a", "prot_b", "error", "title", "raw", "scan", "idx_pre",
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

    df_psm.id = Vector(1:DataFrames.nrow(df_psm))
    DataFrames.select!(df_psm, :id, DataFrames.Not([:id]))
    df_psm.rt = [df_m2[M2I[r.scan], :rt] for r in eachrow(df_psm)]

    if isempty(cfg)
        ele_plink = pLink.read_element() |> NamedTuple
        aa_plink = map(x -> MesMS.calc_mass(x, ele_plink), pLink.read_amino_acid() |> NamedTuple)
        mod_plink = MesMS.mapvalue(x -> x.mass, pLink.read_mod())
        xl_plink = pLink.read_linker() |> NamedTuple
    else
        ele_plink = pLink.read_element(joinpath(cfg, "element.ini")) |> NamedTuple
        aa_plink = map(x -> MesMS.calc_mass(x, ele_plink), pLink.read_amino_acid(joinpath(cfg, "aa.ini")) |> NamedTuple)
        mod_plink = MesMS.mapvalue(x -> x.mass, pLink.read_mod(joinpath(cfg, "modification.ini")))
        xl_plink = pLink.read_linker(joinpath(cfg, "xlink.ini")) |> NamedTuple
    end

    if isempty(cfg_pf)
        ele_pfind = pFind.read_element() |> NamedTuple
        aa_pfind = map(x -> MesMS.calc_mass(x, ele_pfind), pFind.read_amino_acid() |> NamedTuple)
        mod_pfind = MesMS.mapvalue(x -> x.mass, pFind.read_mod())
    else
        ele_pfind = pFind.read_element(joinpath(cfg_pf, "element.ini")) |> NamedTuple
        aa_pfind = map(x -> MesMS.calc_mass(x, ele_pfind), pFind.read_amino_acid(joinpath(cfg_pf, "aa.ini")) |> NamedTuple)
        mod_pfind = MesMS.mapvalue(x -> x.mass, pFind.read_mod(joinpath(cfg_pf, "modification.ini")))
    end

    df_m2.psm = [df_psm[df_psm.scan .== r.id, :id] for r in eachrow(df_m2)]

    !isempty(path_xl) && @info "XL Candidtes loading from " * path_xl
    df_xl = (isempty(path_xl) ? [] : CSV.File(path_xl)) |> DataFrames.DataFrame
    df_xl.id = Vector(1:DataFrames.nrow(df_xl))
    DataFrames.select!(df_xl, :id, DataFrames.Not([:id]))

    !isempty(path_ft) && @info "Feature loading from "* path_ft
    df_ft = (isempty(path_ft) ? [] : CSV.File(path_ft)) |> DataFrames.DataFrame
    df_ft.id = Vector(1:DataFrames.nrow(df_ft))
    DataFrames.select!(df_ft, :id, DataFrames.Not([:id]))

    @info "Target loading from " * path_tg
    df_tg = CSV.File(path_tg) |> DataFrames.DataFrame
    df_tg.id = Vector(1:DataFrames.nrow(df_tg))
    DataFrames.rename!(df_tg, "m/z" => "mz", "t start (min)" => "start", "t stop (min)" => "stop")
    df_tg.start .*= 60
    df_tg.stop .*= 60
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
    
    @async begin
        sleep(4)
        MesMS.open_url("http://$(host):$(port)")
    end
    app = build_app(df_tg, df_xl, df_ft, df_m1, df_m2, df_psm, M2I, ele_plink, aa_plink, mod_plink, xl_plink, ele_pfind, aa_pfind, mod_pfind, ε)
    run_server(app, host, port)
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="TargetView")
    ArgParse.@add_arg_table! settings begin
        "--host"
            help = "hostname"
            metavar = "hostname"
            default = "127.0.0.1"
        "--port"
            help = "port"
            metavar = "port"
            default = "30030"
        "--error"
            help = "m/z error"
            metavar = "ppm"
            default = "20.0"
        "--cfg"
            help = "pLink config directory"
            default = ""
        "--cfg_pf"
            help = "pFind config directory"
            default = ""
        "--linker"
            help = "default linker"
            default = "DSSO"
        "--xl"
            help = "candidate xl list"
            default = ""
        "--ft"
            help = "feature list"
            default = ""
        "--tg"
            help = "target list"
            required = true
        "--psm"
            help = "pLink PSM path"
            required = true
        "--psm_pf"
            help = "pFind PSM path"
            default = ""
        "--fdr"
            help = "FDR threshold (%)"
            default = "Inf"
        "data"
            help = ".ms2 file; .ms1 files should be in the same directory"
            required = true
    end
    args = ArgParse.parse_args(settings)
    sess = prepare(args)
    report(args["data"]; sess...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

julia_main()::Cint = begin
    main()
    return 0
end

end
