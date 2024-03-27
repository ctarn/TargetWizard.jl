module TargetViewX

using Sockets

import ArgParse
import CSV
import DataFrames
import ProgressMeter: @showprogress
import RelocatableFolders: @path
import UniMS
import UniMSUtil: pFind, pLink

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
    ys = map(n -> map(p -> UniMS.max_inten_ε(p, tg.mz + n * Δ / tg.z, ε), df_m1.peaks), -1:2)
    ls = [scatter(x=x, y=ys[2], mode="lines", name="M")]
    push!(ls, scatter(x=x, y=ys[3], mode="lines", name="M + 1 Da"))
    push!(ls, scatter(x=x, y=ys[4], mode="lines", name="M + 2 Da"))
    push!(ls, scatter(x=x, y=ys[1], mode="lines", name="M - 1 Da"))
    p1 = Plot(ls, Layout(; xaxis_title="retention time (s)", yaxis_title="abundance"))
    ls = [scatter(; get_rect(df_ft[i, :], ε/2)..., mode="lines", line_dash="dash", line_color="green", name="FT#$(i)") for i in tg.ft_]
    push!(ls, scatter(x=ones(2) * tg.start, y=[-ε, ε], mode="lines", line_dash="dash", line_color="red", name="RT start"))
    push!(ls, scatter(x=ones(2) * tg.stop, y=[-ε, ε], mode="lines", line_dash="dash", line_color="red", name="RT stop"))
    append!(ls, [scatter(x=df_m2.rt[i:i], y=[UniMS.error_rel(tg.mz, df_m2.mz[i])], name="MS2#$(df_m2.id[i])", customdata=[i], mode="markers", marker_size=4 * (1 + length(df_m2.psm[i]))) for i in tg.m2_all_])
    p2 = Plot(ls, Layout(; yaxis_title="m/z error"))
    p = [p1; p2; p_hit]
    relayout!(p, Layout(Subplots(rows=3, cols=1, row_heights=[2, 1, 2], vertical_spacing=0.02, shared_xaxes=true), clickmode="event+select"), height=600)
    return p
end

build_app(df_tg, df_xl, df_ft, df_m1, df_m2, df_psm, M2I, ele_plink, aa_plink, mod_plink, xl_plink, ele_pfind, aa_pfind, mod_pfind, ε) = begin
    df_tg_tb = DataFrames.select(df_tg, filter(n -> !endswith(n, "_"), names(df_tg)))
    app = dash(; assets_folder=DIR_DATA)
    app.layout = html_div() do
        html_h1("TargetViewX", style=Dict("text-align"=>"center")),
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
        dcc_graph(id="lc_graph", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"TargetViewX_LC"))),
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
        dcc_graph(id="seq_graph", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"TargetViewX_SEQ"))),
        dcc_graph(id="psm_graph", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"TargetViewX_PSM")))
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
            ionss = UniMS.build_ions_crosslink(m2.peaks, seqs, modss, linker, sites, ε, ele_plink, aa_plink, mod_plink)
            p_seq = UniMS.Plotly.seq_crosslink(seqs, modss, sites, ionss)
            p_psm = UniMS.Plotly.spec(m2.peaks, filter(i -> i.peak > 0, vcat(ionss...)))
        elseif r.engine == :pFind
            ions = UniMS.build_ions(m2.peaks, r.pep_a, r.mod_a, ε, ele_pfind, aa_pfind, mod_pfind)
            p_seq = UniMS.Plotly.seq(r.pep_a, r.mod_a, ions)
            p_psm = UniMS.Plotly.spec(m2.peaks, filter(i -> i.peak > 0, ions))
        end
        return p_seq, p_psm
    end
    return app
end

prepare(args) = begin
    path_ms = args["ms"]
    @info "file path of selected data:"
    println("\t$(path_ms)")
    paths_ms_old = reduce(vcat, UniMS.match_path.(args["ms_old"], ".mes")) |> unique |> sort
    @info "file paths of selected original data:"
    foreach(x -> println("$(x[1]):\t$(x[2])"), enumerate(paths_ms_old))
    path_psm = args["psm"]
    out = mkpath(args["out"])
    path_xl = args["xl"]
    path_ft = args["ft"]
    path_psm_pf = args["psm_pf"]
    fmt = args["fmt"] |> Symbol
    linker = Symbol(args["linker"])
    ε = parse(Float64, args["error"]) * 1.0e-6
    fdr = parse(Float64, args["fdr"]) / 100
    decoy = args["decoy"]::Bool
    τ_ms_sim = parse(Float64, args["ms_sim_thres"])
    cfg = args["cfg"]
    cfg_pf = args["cfg_pf"]
    host = parse(IPAddr, args["host"])
    port = parse(Int, args["port"])
    return (; path_ms, paths_ms_old, path_psm, out, path_xl, path_ft, path_psm_pf, fmt, linker, ε, fdr, decoy, τ_ms_sim, cfg, cfg_pf, host, port)
end

psmstr_linear(x) = "$(x.scan)($(round(x.cov; digits=2))):$(x.pep)($(x.mod))"
psmstr_mono(x) = "$(x.scan)($(round(x.cov; digits=2))):$(x.pep)($(x.mod))@$(x.site)"
psmstr_loop(x) = "$(x.scan)($(round(x.cov; digits=2))):$(x.pep)($(x.mod))@$(x.site_a)-$(x.site_b)"
psmstr_link(x) = "$(x.scan)($(round(x.cov_a; digits=2))|$(round(x.cov_b; digits=2))):$(x.pep_a)($(x.mod_a))@$(x.site_a)-$(x.pep_b)($(x.mod_b))@$(x.site_b)"
is_same_xl(a, b) = (a.pep_a == b.pep_a) && (a.pep_b == b.pep_b) && (a.mod_a == b.mod_a) && (a.mod_b == b.mod_b) && (a.site_a == b.site_a) && (a.site_b == b.site_b)
is_same_xl_pepmod(a, b) = (a.pep_a == b.pep_a) && (a.pep_b == b.pep_b) && (a.mod_a == b.mod_a) && (a.mod_b == b.mod_b)

unify_mods_str(s) = (!ismissing(s) && startswith(s, "Any[") && endswith(s, "]")) ? s[5:end-1] : s

process(path; path_ms, paths_ms_old, path_psm, out, path_xl, path_ft, path_psm_pf, fmt, linker, ε, fdr, decoy, τ_ms_sim, cfg, cfg_pf, host, port) = begin
    M = UniMS.read_ms(path_ms)
    df_m1 = map(m -> (; m.id, rt=m.retention_time, m.peaks), M.MS1) |> DataFrames.DataFrame
    df_m2 = map(m -> (; m.id, mz=m.activation_center, rt=m.retention_time, m.peaks), M.MS2) |> DataFrames.DataFrame
    M2I = map(x -> x[2] => x[1], enumerate(df_m2.id)) |> Dict

    M_old = map(p -> splitext(basename(p))[1] => UniMS.dict_by_id(UniMS.read_ms(p).MS2), paths_ms_old) |> Dict

    if isempty(cfg)
        ele_plink = pLink.read_element() |> NamedTuple
        aa_plink = map(x -> UniMS.mass(x, ele_plink), pLink.read_amino_acid() |> NamedTuple)
        mod_plink = UniMS.mapvalue(x -> x.mass, pLink.read_modification())
        xl_plink = pLink.read_linker() |> NamedTuple
    else
        ele_plink = pLink.read_element(joinpath(cfg, "element.ini")) |> NamedTuple
        aa_plink = map(x -> UniMS.mass(x, ele_plink), pLink.read_amino_acid(joinpath(cfg, "aa.ini")) |> NamedTuple)
        mod_plink = UniMS.mapvalue(x -> x.mass, pLink.read_modification(joinpath(cfg, "modification.ini")))
        xl_plink = pLink.read_linker(joinpath(cfg, "xlink.ini")) |> NamedTuple
    end

    if isempty(cfg_pf)
        ele_pfind = pFind.read_element() |> NamedTuple
        aa_pfind = map(x -> UniMS.mass(x, ele_pfind), pFind.read_amino_acid() |> NamedTuple)
        mod_pfind = UniMS.mapvalue(x -> x.mass, pFind.read_modification())
    else
        ele_pfind = pFind.read_element(joinpath(cfg_pf, "element.ini")) |> NamedTuple
        aa_pfind = map(x -> UniMS.mass(x, ele_pfind), pFind.read_amino_acid(joinpath(cfg_pf, "aa.ini")) |> NamedTuple)
        mod_pfind = UniMS.mapvalue(x -> x.mass, pFind.read_modification(joinpath(cfg_pf, "modification.ini")))
    end

    dfs = pLink.read_psm_full(path_psm)
    df_psm = dfs.xl
    df_linear = dfs.linear
    df_mono = dfs.mono
    df_loop = dfs.loop

    for df in [df_psm, df_linear, df_mono, df_loop]
        df.engine .= :pLink
        filter!(r -> r.fdr .≤ fdr, df)
        !decoy && filter!(r -> r.td == :TT || r.td == :T, df)
        ("linker" ∉ names(df)) && (df.linker .= linker)
    end
    
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

    ion_syms = ["b", "y"]
    ion_types = map(i -> getfield(UniMS, Symbol("ion_$(i)")), ion_syms)
    M_ = [splitext(basename(path_ms))[1] => UniMS.dict_by_id(M.MS2)] |> Dict
    calc_cov_crosslink!(df_psm, M_, ε, ion_syms, ion_types, ele_plink, aa_plink, mod_plink, xl_plink)
    calc_cov_linear!(df_linear, M_, ε, ion_syms, ion_types, ele_plink, aa_plink, mod_plink)
    calc_cov_monolink!(df_mono, M_, ε, ion_syms, ion_types, ele_plink, aa_plink, mod_plink, xl_plink)
    calc_cov_looplink!(df_loop, M_, ε, ion_syms, ion_types, ele_plink, aa_plink, mod_plink, xl_plink)

    df_psm.cov_min = min.(df_psm.cov_a, df_psm.cov_b)
    df_psm.credible = map(eachrow(df_psm)) do r
        r.cov_a_ion_y ≥ 0.6 && r.cov_b_ion_y ≥ 0.6 && r.cov_a_ion_b ≥ 0.4 && r.cov_b_ion_b ≥ 0.4
    end
    df_linear.credible = map(r -> r.cov_ion_y ≥ 0.6 && r.cov_ion_b ≥ 0.4, eachrow(df_linear))
    df_mono.credible = map(r -> r.cov_ion_y ≥ 0.6 && r.cov_ion_b ≥ 0.4, eachrow(df_mono))
    df_loop.credible = map(r -> r.cov_ion_y ≥ 0.6 && r.cov_ion_b ≥ 0.4, eachrow(df_loop))

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

    for df in [df_psm, df_linear, df_mono, df_loop]
        df.id = Vector(1:size(df, 1))
        DataFrames.select!(df, :id, DataFrames.Not([:id]))
        df.rt = [df_m2[M2I[r.scan], :rt] for r in eachrow(df)]
    end

    df_m2.psm = [df_psm[df_psm.scan .== r.id, :id] for r in eachrow(df_m2)]
    df_m2.psm_linear = [df_linear[df_linear.scan .== r.id, :id] for r in eachrow(df_m2)]
    df_m2.psm_mono = [df_mono[df_mono.scan .== r.id, :id] for r in eachrow(df_m2)]
    df_m2.psm_loop = [df_loop[df_loop.scan .== r.id, :id] for r in eachrow(df_m2)]

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
    "mod_a" ∈ names(df_tg) && (df_tg.mod_a = parse.(Array{UniMS.Mod}, unify_mods_str.(df_tg.mod_a)))
    "mod_b" ∈ names(df_tg) && (df_tg.mod_b = parse.(Array{UniMS.Mod}, unify_mods_str.(df_tg.mod_b)))

    @info "XL Candidtes mapping"
    tmp = sort!([(x.mz::Float64, x.id::Int) for x in eachrow(df_xl)])
    mzs = map(x -> x[1], tmp)
    ids = map(x -> x[2], tmp)
    df_tg.xl_ = [sort(filter(x -> df_xl[x, :z] == r.z, ids[UniMS.argquery_ε(mzs, r.mz, ε)])) for r in eachrow(df_tg)]

    @info "Feature mapping"
    tmp = sort!([(x.mz::Float64, x.id::Int) for x in eachrow(df_ft)])
    mzs = map(x -> x[1], tmp)
    ids = map(x -> x[2], tmp)
    df_tg.ft_ = [sort(filter(x -> df_ft[x, :z] == r.z, ids[UniMS.argquery_ε(mzs, r.mz, ε)])) for r in eachrow(df_tg)]
    df_tg.n_ft = length.(df_tg.ft_)

    Ks = ["", "_allsim", "_all"]

    calc_sim(dda, tda) = map(p -> !isempty(UniMS.argquery_ε(tda.peaks, p.mz, ε)), M_old[dda.file][dda.scan].peaks) |> mean

    @info "MS2 mapping"
    tmp = sort!([(x.mz::Float64, x.id::Int) for x in eachrow(df_m2)])
    mzs = map(x -> x[1], tmp)
    ids = map(x -> x[2], tmp)
    df_tg.m2_all_ = [map(x -> M2I[x], sort(ids[UniMS.argquery_ε(mzs, r.mz, ε)])) for r in eachrow(df_tg)]
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
    filter_plink(x) = filter(i -> df_psm.engine[i] == :pLink, x)
    filter_pfind(x) = filter(i -> df_psm.engine[i] == :pFind, x)
    for K in Ks
        df_tg[!, "psm$(K)_"] = [vcat(df_m2[r["m2$(K)_"], :psm]...) for r in eachrow(df_tg)]
        df_tg[!, "n_psm$(K)"] = length.(df_tg[!, "psm$(K)_"])

        df_tg[!, "psm_linear$(K)_"] = [vcat(df_m2[r["m2$(K)_"], :psm_linear]...) for r in eachrow(df_tg)]
        df_tg[!, "psm_mono$(K)_"] = [vcat(df_m2[r["m2$(K)_"], :psm_mono]...) for r in eachrow(df_tg)]
        df_tg[!, "psm_loop$(K)_"] = [vcat(df_m2[r["m2$(K)_"], :psm_loop]...) for r in eachrow(df_tg)]

        df_tg[!, "n_psm_plink$(K)"] = length.(map(filter_plink, df_tg[!, "psm$(K)_"]))
        df_tg[!, "n_psm_pfind$(K)"] = length.(map(filter_pfind, df_tg[!, "psm$(K)_"]))
    end

    ns = filter(n -> !endswith(n, '_'), names(df_tg))
    DataFrames.select!(df_tg, ns, DataFrames.Not(ns))
    
    df_tg_ext = DataFrames.DataFrame(df_tg)
    for (k, s, d) in zip(["", "_linear", "_mono", "_loop"], [psmstr_link, psmstr_linear, psmstr_mono, psmstr_loop], [df_psm, df_linear, df_mono, df_loop])
        for K in Ks
            df_tg_ext[!, "iden$(k)$(K)"] = map(df_tg_ext[!, "psm$(k)$(K)_"]) do psms
                map(s, eachrow(d[psms, :])) |> xs -> join(xs, ";")
            end
            df_tg_ext[!, "have_iden$(k)$(K)"] = .!isempty.(df_tg_ext[!, "iden$(k)$(K)"])
            df_tg_ext[!, "iden$(k)$(K)_credible"] = map(df_tg_ext[!, "psm$(k)$(K)_"]) do psms
                map(s, eachrow(d[filter(i -> d.credible[i], psms), :])) |> xs -> join(xs, ";")
            end
            df_tg_ext[!, "have_iden$(k)$(K)_credible"] = .!isempty.(df_tg_ext[!, "iden$(k)$(K)_credible"])
        end
    end

    for K in Ks
        for (f, n) in zip([is_same_xl, is_same_xl_pepmod], ["iden", "pepmod"])
            df_tg_ext[!, "same_$(n)$(K)"] = map(eachrow(df_tg_ext)) do r
                filter(eachrow(df_psm[r["psm$(K)_"], :])) do s
                    (s.engine != :pLink) && return false
                    return f(r, s)
                end .|> psmstr_link |> xs -> join(xs, ";")
            end
            df_tg_ext[!, "have_same_$(n)$(K)"] = .!isempty.(df_tg_ext[!, "same_$(n)$(K)"])

            df_tg_ext[!, "same_$(n)$(K)_credible"] = map(eachrow(df_tg_ext)) do r
                filter(eachrow(df_psm[r["psm$(K)_"], :])) do s
                    (s.engine != :pLink) && return false
                    (!s.credible) && return false
                    return f(r, s)
                end .|> psmstr_link |> xs -> join(xs, ";")
            end
            df_tg_ext[!, "have_same_$(n)$(K)_credible"] = .!isempty.(df_tg_ext[!, "same_$(n)$(K)_credible"])

            vs = map(eachrow(df_tg_ext)) do r
                b = nothing
                for s in eachrow(df_psm[r["psm$(K)_"], :])
                    (s.engine != :pLink) && continue
                    isnothing(b) && (b = s)
                    (s.cov_min ≤ b.cov_min) && continue
                    b = s
                end
                return isnothing(b) ? ("", false) : (psmstr_link(b), f(r, b))
            end
            df_tg_ext[!, "best_iden$(K)"] = first.(vs)
            df_tg_ext[!, "best_iden$(K)_is_same_$(n)"] = last.(vs)
        end
    end

    UniMS.safe_save(p -> CSV.write(p, df_tg_ext), joinpath(out, "$(basename(splitext(path_ms)[1])).tg.TargetViewX.csv"))
    UniMS.safe_save(p -> CSV.write(p, df_psm), joinpath(out, "$(basename(splitext(path_ms)[1])).crosslink.TargetViewX.csv"))
    UniMS.safe_save(p -> CSV.write(p, df_linear), joinpath(out, "$(basename(splitext(path_ms)[1])).linear.TargetViewX.csv"))
    UniMS.safe_save(p -> CSV.write(p, df_mono), joinpath(out, "$(basename(splitext(path_ms)[1])).monolink.TargetViewX.csv"))
    UniMS.safe_save(p -> CSV.write(p, df_loop), joinpath(out, "$(basename(splitext(path_ms)[1])).looplink.TargetViewX.csv"))

    @async begin
        sleep(4)
        UniMS.open_url("http://$(host):$(port)")
    end
    app = build_app(df_tg, df_xl, df_ft, df_m1, df_m2, df_psm, M2I, ele_plink, aa_plink, mod_plink, xl_plink, ele_pfind, aa_pfind, mod_pfind, ε)
    run_server(app, host, port)
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="TargetViewX")
    ArgParse.@add_arg_table! settings begin
        "target"
            help = "target list"
            required = true
        "--ms"
            help = ".mes or .ms1/2 file; .ms2/1 file should be in the same directory for .ms1/2"
            required = true
        "--ms_old"
            help = "origianl .mes or .ms1/2 files; .ms2/1 files should be in the same directory for .ms1/2"
            nargs = '+'
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
        "--decoy"
            help = "preserve decoy identifications"
            action = :store_true
        "--ms_sim_thres"
            help = "threshold of MS similarity"
            default = "0.5"
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
