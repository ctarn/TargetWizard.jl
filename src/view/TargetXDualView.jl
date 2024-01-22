module TargetXDualView

using Sockets

import ArgParse
import CSV
import DataFrames
import MesMS
import MesUtil: pLink
import ProgressMeter: @showprogress
import RelocatableFolders: @path

using Dash
using PlotlyBase

include("../util.jl")

const DIR_DATA = @path joinpath(@__DIR__, "../../data/dash")

Δ = 1.00335

plot_lc(tg, dfs_ft, dfs_m1, dfs_m2, ε) = begin
    ps = map(["A", "B"], dfs_m1) do s, df_m1
        x = df_m1.rt
        ys = map(n -> map(p -> MesMS.max_inten_ε(p, tg.mz + n * Δ / tg.z, ε), df_m1.peaks), -1:2)
        ls = [scatter(x=x, y=ys[2], mode="lines", name="$(s) M")]
        push!(ls, scatter(x=x, y=ys[3], mode="lines", name="$(s) M + 1 Da"))
        push!(ls, scatter(x=x, y=ys[4], mode="lines", name="$(s) M + 2 Da"))
        push!(ls, scatter(x=x, y=ys[1], mode="lines", name="$(s) M - 1 Da"))
        return Plot(ls, Layout(; yaxis_title="abundance"))
    end
    ls = [scatter(x=ones(2) * tg.start, y=[-ε, ε], mode="lines", line_dash="dash", line_color="red", name="RT start")]
    push!(ls, scatter(x=ones(2) * tg.stop, y=[-ε, ε], mode="lines", line_dash="dash", line_color="red", name="RT stop"))
    append!(ls, [scatter(;
        x=dfs_m2[1].rt[i:i], y=[MesMS.error_rel(tg.mz, dfs_m2[1].mz[i])], name="A MS2#$(dfs_m2[1].id[i])",
        customdata=[(:A, i)], mode="markers", marker_size=4 * (1 + length(dfs_m2[1].psm[i])), marker_symbol="circle"
    ) for i in tg.m2_all_a_])
    append!(ls, [scatter(
        x=dfs_m2[2].rt[i:i], y=[MesMS.error_rel(tg.mz, dfs_m2[2].mz[i])], name="B MS2#$(dfs_m2[2].id[i])",
        customdata=[(:B, i)], mode="markers", marker_size=4 * (1 + length(dfs_m2[2].psm[i])), marker_symbol="diamond"
    ) for i in tg.m2_all_b_])
    append!(ls, [scatter(;
        x=dfs_m2[1].rt[i:i], y=[0], name="A Ext. MS2#$(dfs_m2[1].id[i])",
        customdata=[(:A, i)], mode="markers", marker_size=4 * (1 + length(dfs_m2[1].psm[i])), marker_symbol="circle"
    ) for i in tg.m2_ext_all_a_])
    append!(ls, [scatter(
        x=dfs_m2[2].rt[i:i], y=[0], name="B Ext. MS2#$(dfs_m2[2].id[i])",
        customdata=[(:B, i)], mode="markers", marker_size=4 * (1 + length(dfs_m2[2].psm[i])), marker_symbol="diamond"
    ) for i in tg.m2_ext_all_b_])
    append!(ls, [scatter(;
        x=[dfs_ft[1][i, :rtime_start], dfs_ft[1][i, :rtime_stop]], y=ε*2/3 .* ones(2), name="A FT#$(i)",
        mode="lines", line_color="green",
    ) for i in tg.ft_a_])
    append!(ls, [scatter(;
        x=[dfs_ft[2][i, :rtime_start], dfs_ft[2][i, :rtime_stop]], y=-ε*2/3 .* ones(2), name="B FT#$(i)",
        mode="lines", line_color="blue",
    ) for i in tg.ft_b_])
    p3 = Plot(ls, Layout(; xaxis_title="retention time (s)", yaxis_title="m/z error"))
    p = [ps[1]; ps[2]; p3]
    relayout!(p, Layout(Subplots(rows=3, cols=1, vertical_spacing=0.005, shared_xaxes=true), clickmode="event+select"), height=600)
    return p
end

plot_psm(iden, spec, ε, ele_plink, aa_plink, mod_plink, xl_plink) = begin
    seqs = (iden.pep_a, iden.pep_b)
    modss = (iden.mod_a, iden.mod_b)
    linker = xl_plink[Symbol(iden.linker)]
    sites = (iden.site_a, iden.site_b)
    ionss = MesMS.build_ions_crosslink(spec.peaks, seqs, modss, linker, sites, ε, ele_plink, aa_plink, mod_plink)
    p_seq = MesMS.Plotly.seq_crosslink(seqs, modss, sites, ionss)
    p_psm = MesMS.Plotly.spec(spec.peaks, filter(i -> i.peak > 0, vcat(ionss...)))
    relayout!(p_seq, Dict(:margin=>Dict(:l=>0, :r=>0, :t=>0, :b=>0, :pad=>0)))
    relayout!(p_psm, Dict(:margin=>Dict(:l=>0, :r=>0, :t=>0, :b=>0, :pad=>0)))
    return p_seq, p_psm
end

build_app(df_tg, df_xl, dfs_ft, dfs_m1, dfs_m2, dfs_psm, M2Is, ele_plink, aa_plink, mod_plink, xl_plink, ε) = begin
    df_tg_tab = DataFrames.select(df_tg, filter(n -> !endswith(n, "_"), names(df_tg)))
    app = dash(; assets_folder=DIR_DATA)
    app.layout = html_div() do
        html_h1("TargetXDualView", style=Dict("text-align"=>"center")),
        dash_datatable(
            id="tg_table",
            style_table=Dict("min-width"=>"100%", "overflow-x"=>"auto"),
            style_cell=Dict("overflow"=>"hidden", "text-overflow"=>"ellipsis", "min-width"=>"64px", "max-width"=>"256px"),
            columns=[(; name=i, id=i) for i in names(df_tg_tab)],
            data=Dict.(pairs.(eachrow(string.(df_tg_tab)))),
            filter_action="native",
            sort_action="native",
            sort_mode="multi",
            row_selectable="single",
            page_action="native",
            page_size=10,
            export_format="csv",
            export_headers="display",
        ),
        html_h4("Candidate Crosslink List of Selected Target"),
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
        dcc_graph(id="lc_graph", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"TargetXDualView_LC"))),
        html_h4("PSM List of Selected Data A MS2(s) (may differ from the selected target)"),
        dash_datatable(
            id="psm_table_a",
            style_table=Dict("min-width"=>"100%", "overflow-x"=>"auto"),
            style_cell=Dict("overflow"=>"hidden", "min-width"=>"64px",),
            columns=[(; name=i, id=i) for i in names(dfs_psm[1])],
            filter_action="native",
            sort_action="native",
            sort_mode="multi",
            row_selectable="single",
            page_action="native",
            page_size=4,
            export_format="csv",
            export_headers="display",
        ),
        html_h4("PSM List of Selected Data B MS2(s) (may differ from the selected target)"),
        dash_datatable(
            id="psm_table_b",
            style_table=Dict("min-width"=>"100%", "overflow-x"=>"auto"),
            style_cell=Dict("overflow"=>"hidden", "min-width"=>"64px",),
            columns=[(; name=i, id=i) for i in names(dfs_psm[2])],
            filter_action="native",
            sort_action="native",
            sort_mode="multi",
            row_selectable="single",
            page_action="native",
            page_size=4,
            export_format="csv",
            export_headers="display",
        ),
        html_div(className="opts") do
            html_div(className="opt") do
                html_span("frag. error (A):", className="name"),
                html_div(className="value") do
                    dcc_input(id="error_frag_a", className="input", value=20.0, type="number", placeholder="error"),
                    html_span("ppm", className="unit")
                end
            end,
            html_div(className="opt") do
                html_span("frag. error (B):", className="name"),
                html_div(className="value") do
                    dcc_input(id="error_frag_b", className="input", value=20.0, type="number", placeholder="error"),
                    html_span("ppm", className="unit")
                end
            end
        end,
        dcc_graph(id="seq_graph_a", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"TargetXDualView_SEQ_A"))),
        dcc_graph(id="seq_graph_b", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"TargetXDualView_SEQ_B"))),
        dcc_graph(id="psm_graph_a", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"TargetXDualView_PSM_A"))),
        dcc_graph(id="psm_graph_b", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"TargetXDualView_PSM_B")))
    end

    callback!(app,
        Output("xl_table", "data"),
        Output("lc_graph", "figure"),
        Input("tg_table", "derived_virtual_data"),
        Input("tg_table", "derived_virtual_selected_rows"),
    ) do v1, v2
        id = parse(Int, v1[v2[begin] + 1].id)
        tg = df_tg[id, :]
        table_data = Dict.(pairs.(eachrow(string.(df_xl[tg.xl_, :]))))
        fig = plot_lc(tg, dfs_ft, dfs_m1, dfs_m2, ε)
        return table_data, fig
    end

    callback!(app,
        Output("psm_table_a", "data"),
        Output("psm_table_b", "data"),
        Input("lc_graph", "selectedData"),
    ) do v
        ids = map(p -> p.customdata, v.points)
        ids_a = map(i -> i[2], filter(i -> i[1] == "A", ids))
        ids_b = map(i -> i[2], filter(i -> i[1] == "B", ids))
        psms_a = isempty(ids_a) ? DataFrames.DataFrame() : dfs_psm[1][vcat(dfs_m2[1][ids_a, :psm]...), :]
        psms_b = isempty(ids_b) ? DataFrames.DataFrame() : dfs_psm[2][vcat(dfs_m2[2][ids_b, :psm]...), :]
        return Dict.(pairs.(eachrow(string.(psms_a)))), Dict.(pairs.(eachrow(string.(psms_b))))
    end

    callback!(app,
        Output("seq_graph_a", "figure"),
        Output("psm_graph_a", "figure"),
        Input("psm_table_a", "derived_virtual_data"),
        Input("psm_table_a", "derived_virtual_selected_rows"),
        Input("error_frag_a", "value"),
    ) do v1, v2, ε2
        ε2 = ε2 * 1.0e-6
        id = parse(Int, v1[v2[begin] + 1].id)
        iden = dfs_psm[1][id, :]
        spec = dfs_m2[1][M2Is[1][iden.scan], :]
        p_seq, p_psm = plot_psm(iden, spec, ε2, ele_plink, aa_plink, mod_plink, xl_plink)
        relayout!(p_seq, Dict(:title=>"Data A"))
        relayout!(p_psm, Dict(:title=>"Data A"))
        return p_seq, p_psm
    end

    callback!(app,
        Output("seq_graph_b", "figure"),
        Output("psm_graph_b", "figure"),
        Input("psm_table_b", "derived_virtual_data"),
        Input("psm_table_b", "derived_virtual_selected_rows"),
        Input("error_frag_b", "value"),
    ) do v1, v2, ε2
        ε2 = ε2 * 1.0e-6
        id = parse(Int, v1[v2[begin] + 1].id)
        iden = dfs_psm[2][id, :]
        spec = dfs_m2[2][M2Is[2][iden.scan], :]
        p_seq, p_psm = plot_psm(iden, spec, ε2, ele_plink, aa_plink, mod_plink, xl_plink)
        relayout!(p_seq, Dict(:title=>"Data B"))
        relayout!(p_psm, Dict(:title=>"Data B"))
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
    fmt = args["fmt"] |> Symbol
    linker = Symbol(args["linker"])
    ε = parse(Float64, args["error"]) * 1.0e-6
    fdr = parse(Float64, args["fdr"]) / 100
    cfg = args["cfg"]
    host = parse(IPAddr, args["host"])
    port = parse(Int, args["port"])
    return (; path_ms, path_psm, out, path_xl, path_ft, fmt, linker, ε, fdr, cfg, host, port)
end

process(path; path_ms, path_psm, out, path_xl, path_ft, fmt, linker, ε, fdr, cfg, host, port) = begin
    Ms = map(MesMS.read_ms, path_ms)
    dfs_m1 = map(Ms) do M
        map(m -> (; m.id, rt=m.retention_time, m.peaks), M.MS1)
    end .|> DataFrames.DataFrame
    dfs_m2 = map(Ms) do M
        map(m -> (; m.id, mz=m.activation_center, rt=m.retention_time, m.peaks), M.MS2)
    end .|> DataFrames.DataFrame
    M2Is = map(df -> map(x -> x[2] => x[1], enumerate(df.id)), dfs_m2) .|> Dict

    dfs_psm = map(path_psm, dfs_m2, M2Is) do p, df_m2, M2I
        df_psm = pLink.read_psm_full(p).xl
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
        df_psm.id = Vector(1:size(df_psm, 1))
        DataFrames.select!(df_psm, :id, DataFrames.Not([:id]))
        df_psm.rt = [df_m2[M2I[r.scan], :rt] for r in eachrow(df_psm)]
        df_psm.mz = round.(df_psm.mz, digits=6)
        return df_psm
    end

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

    dfs_m2[1].psm = [dfs_psm[1][dfs_psm[1].scan .== i, :id] for i in dfs_m2[1].id]
    dfs_m2[2].psm = [dfs_psm[2][dfs_psm[2].scan .== i, :id] for i in dfs_m2[2].id]

    !isempty(path_xl) && @info "XL Candidtes loading from " * path_xl
    df_xl = DataFrames.DataFrame(isempty(path_xl) ? [] : CSV.File(path_xl))
    df_xl.id = Vector(1:size(df_xl, 1))
    DataFrames.select!(df_xl, :id, DataFrames.Not([:id]))

    dfs_ft = map(path_ft) do p
        !isempty(p) && @info "Feature loading from " * p
        df_ft = DataFrames.DataFrame(isempty(p) ? [] : CSV.File(p))
        df_ft.id = Vector(1:size(df_ft, 1))
        DataFrames.select!(df_ft, :id, DataFrames.Not([:id]))
        return df_ft
    end

    @info "Target loading from " * path
    df_tg = DataFrames.DataFrame(CSV.File(path))
    df_tg.id = Vector(1:size(df_tg, 1))
    parse_target_list!(df_tg, fmt)
    DataFrames.select!(df_tg, [:id, :mz, :z, :start, :stop], DataFrames.Not([:id, :mz, :z, :start, :stop]))
    df_tg.start = round.(df_tg.start; digits=2)
    df_tg.stop = round.(df_tg.stop; digits=2)

    @info "XL Candidtes mapping"
    tmp = sort!([(x.mz::Float64, x.id::Int) for x in eachrow(df_xl)])
    mzs = map(x -> x[1], tmp)
    ids = map(x -> x[2], tmp)
    df_tg.xl_ = [sort(filter(x -> df_xl[x, :z] == r.z, ids[MesMS.argquery_ε(mzs, r.mz, ε)])) for r in eachrow(df_tg)]

    @info "Feature mapping"
    df_tg.ft_a_, df_tg.ft_b_ = map(dfs_ft) do df_ft
        tmp = sort!([(x.mz::Float64, x.id::Int) for x in eachrow(df_ft)])
        mzs = map(x -> x[1], tmp)
        ids = map(x -> x[2], tmp)
        return [sort(filter(x -> df_ft[x, :z] == r.z, ids[MesMS.argquery_ε(mzs, r.mz, ε)])) for r in eachrow(df_tg)]
    end
    df_tg.n_ft_a = length.(df_tg.ft_a_)
    df_tg.n_ft_b = length.(df_tg.ft_b_)

    @info "PSM mapping"
    df_tg.psm_all_a_, df_tg.psm_all_b_ = map(dfs_psm) do df_psm
        tmp = sort!([(x.mz::Float64, x.id::Int) for x in eachrow(df_psm)])
        mzs = map(x -> x[1], tmp)
        ids = map(x -> x[2], tmp)
        return [sort(filter(x -> df_psm[x, :z] == r.z, ids[MesMS.argquery_ε(mzs, r.mz, ε)])) for r in eachrow(df_tg)]
    end
    df_tg.psm_a_ = [filter(i -> r.start ≤ dfs_psm[1].rt[i] ≤ r.stop, r.psm_all_a_) for r in eachrow(df_tg)]
    df_tg.psm_b_ = [filter(i -> r.start ≤ dfs_psm[2].rt[i] ≤ r.stop, r.psm_all_b_) for r in eachrow(df_tg)]
    df_tg.n_psm_all_a = length.(df_tg.psm_all_a_)
    df_tg.n_psm_all_b = length.(df_tg.psm_all_b_)
    df_tg.n_psm_a = length.(df_tg.psm_a_)
    df_tg.n_psm_b = length.(df_tg.psm_b_)

    @info "MS2 mapping"
    df_tg.m2_all_a_, df_tg.m2_all_b_ = map(dfs_m2, M2Is) do df_m2, M2I
        tmp = sort!([(x.mz::Float64, x.id::Int) for x in eachrow(df_m2)])
        mzs = map(x -> x[1], tmp)
        ids = map(x -> x[2], tmp)
        return [map(x -> M2I[x], sort(ids[MesMS.argquery_ε(mzs, r.mz, ε)])) for r in eachrow(df_tg)]
    end
    df_tg.m2_a_ = [filter(i -> r.start ≤ dfs_m2[1].rt[i] ≤ r.stop, r.m2_all_a_) for r in eachrow(df_tg)]
    df_tg.m2_b_ = [filter(i -> r.start ≤ dfs_m2[2].rt[i] ≤ r.stop, r.m2_all_b_) for r in eachrow(df_tg)]
    df_tg.n_m2_all_a = length.(df_tg.m2_all_a_)
    df_tg.n_m2_all_b = length.(df_tg.m2_all_b_)
    df_tg.n_m2_a = length.(df_tg.m2_a_)
    df_tg.n_m2_b = length.(df_tg.m2_b_)
    df_tg.m2_ext_all_a_ = [map(x -> M2Is[1][x], sort(unique(dfs_psm[1][r.psm_all_a_, :scan]))) for r in eachrow(df_tg)]
    df_tg.m2_ext_all_b_ = [map(x -> M2Is[2][x], sort(unique(dfs_psm[2][r.psm_all_b_, :scan]))) for r in eachrow(df_tg)]
    df_tg.m2_ext_all_a_ = setdiff.(df_tg.m2_ext_all_a_, df_tg.m2_all_a_)
    df_tg.m2_ext_all_b_ = setdiff.(df_tg.m2_ext_all_b_, df_tg.m2_all_b_)
    df_tg.m2_ext_a_ = [filter(i -> r.start ≤ dfs_m2[1].rt[i] ≤ r.stop, r.m2_ext_all_a_) for r in eachrow(df_tg)]
    df_tg.m2_ext_b_ = [filter(i -> r.start ≤ dfs_m2[2].rt[i] ≤ r.stop, r.m2_ext_all_b_) for r in eachrow(df_tg)]
    df_tg.n_m2_ext_all_a = length.(df_tg.m2_ext_all_a_)
    df_tg.n_m2_ext_all_b = length.(df_tg.m2_ext_all_b_)
    df_tg.n_m2_ext_a = length.(df_tg.m2_ext_a_)
    df_tg.n_m2_ext_b = length.(df_tg.m2_ext_b_)

    df_tg.m2_all_a_id_ = [map(x -> dfs_m2[1].id[x], r.m2_all_a_) for r in eachrow(df_tg)]
    df_tg.m2_all_b_id_ = [map(x -> dfs_m2[2].id[x], r.m2_all_b_) for r in eachrow(df_tg)]
    df_tg.m2_a_id_ = [map(x -> dfs_m2[1].id[x], r.m2_a_) for r in eachrow(df_tg)]
    df_tg.m2_b_id_ = [map(x -> dfs_m2[2].id[x], r.m2_b_) for r in eachrow(df_tg)]
    df_tg.m2_ext_all_a_id_ = [map(x -> dfs_m2[1].id[x], r.m2_ext_all_a_) for r in eachrow(df_tg)]
    df_tg.m2_ext_all_b_id_ = [map(x -> dfs_m2[2].id[x], r.m2_ext_all_b_) for r in eachrow(df_tg)]
    df_tg.m2_ext_a_id_ = [map(x -> dfs_m2[1].id[x], r.m2_ext_a_) for r in eachrow(df_tg)]
    df_tg.m2_ext_b_id_ = [map(x -> dfs_m2[2].id[x], r.m2_ext_b_) for r in eachrow(df_tg)]

    ns = filter(n -> !endswith(n, '_'), names(df_tg))
    DataFrames.select!(df_tg, ns, DataFrames.Not(ns))

    MesMS.safe_save(p -> CSV.write(p, df_tg), joinpath(out, "$(basename(splitext(path)[1])).TargetXDualView.csv"))

    @async begin
        sleep(4)
        MesMS.open_url("http://$(host):$(port)")
    end
    app = build_app(df_tg, df_xl, dfs_ft, dfs_m1, dfs_m2, dfs_psm, M2Is, ele_plink, aa_plink, mod_plink, xl_plink, ε)
    run_server(app, host, port)
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="TargetXDualView")
    ArgParse.@add_arg_table! settings begin
        "target"
            help = "target list"
            required = true
        "--ms"
            help = "two .mes or .ms1/2 files; .ms2/1 files should be in the same directory for .ms1/2"
            required = true
            nargs = 2
        "--psm"
            help = "two pLink PSM paths"
            required = true
            nargs = 2
        "--out"
            help = "output directory"
            default = "./out/"
            metavar = "./out/"
        "--xl"
            help = "corresponding candidate xl list of targets"
            default = ""
        "--ft"
            help = "two feature lists"
            default = ["", ""]
            nargs = 2
        "--fmt", "-f"
            help = "target list format: auto, TW, TmQE, TmFu"
            metavar = "auto|TW|TmQE|TmFu"
            default = "auto"
        "--linker"
            help = "default linker"
            metavar = "DSSO"
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
