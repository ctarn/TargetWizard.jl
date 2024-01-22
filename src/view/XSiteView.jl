module XSiteView

using Printf
using Statistics
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

build_ions(seq, mods, tab_ele, tab_aa, tab_mod; types=[(MesMS.ion_b, 1:3), (MesMS.ion_y, 1:3)]) = begin
    ts = [(; t..., mass=MesMS.mass(t.Δ, tab_ele), zs) for (t, zs) in types]
    base = MesMS.calc_ion_base(MesMS.lsum_seq(seq, mods, tab_aa, tab_mod))
    ions = [MesMS.calc_ion(i, t.mass, z, t.type, t.sym, t.color) for t in ts for z in t.zs for i in base[t.part]]
    return ions
end

_plot_seq!(ls, x, y, seq, mods, ions, font) = begin
    colors = [:black for _ in seq]
    foreach(m -> colors[m.site] = :red, mods)
    for (i, aa) in enumerate(seq)
        push!(ls, scatter(x=[x + i - 0.5], y=[y], mode="text", name="", text=string(aa), textposition="middle center", textfont=attr(size=font*2, color=colors[i])))
    end
    n = Dict()
    for i in ions
        n[(i.part, i.loc)] = get(n, (i.part, i.loc), 0) + 1
        if i.part == :l
            push!(ls, scatter(x=[x + i.loc], y=[y - (n[(i.part, i.loc)] / 2 + 1.5)], mode="text", name="", text="┛", hovertext=i.text, customdata=[(i.seq, i.text, i.mz, i.z)], textposition="top left", textfont=attr(size=font, color=i.color)))
        elseif i.part == :r
            push!(ls, scatter(x=[x + i.loc], y=[y + (n[(i.part, i.loc)] / 2 + 2.0)], mode="text", name="", text="┏", hovertext=i.text, customdata=[(i.seq, i.text, i.mz, i.z)], textposition="bottom right", textfont=attr(size=font, color=i.color)))
        end
    end
end

seq_crosslink(seqs, modss, site_pairs, ionss; font=18) = begin
    seqs = collect.(seqs)
    lα, lβ = mean(extrema(first.(site_pairs))), mean(extrema(last.(site_pairs)))
    xs = max(lα, lβ) .- [lα, lβ] .+ [0.25, -0.25]
    ys = [4, -4]
    ls = [scatter(x=[a, b] .+ xs .- 0.5, y=copysign.(abs.(ys) .- 2, ys), mode="lines", name="") for (a, b) in site_pairs]
    for (x, y, seq, mods, ions) in zip(xs, ys, seqs, modss, ionss)
        _plot_seq!(ls, x, y, seq, mods, ions, font)
    end
    return Plot(ls, Layout(; showlegend=false, yaxis=attr(showticklabels=false, range=(-8, 8)), xaxis=attr(showticklabels=false), height=300))
end

get_inten(mz, ps, ε) = maximum(p -> p.inten, MesMS.query_ε(ps, mz, ε); init=0.0)

smooth(x, k) = begin
    n = length(k)
    x_ = vcat(zeros(eltype(x), n ÷ 2), x, zeros(eltype(x), n - n ÷ 2 - 1))
    return map(i -> sum(x_[i:i+n-1] .* k), eachindex(x)) ./ sum(k)
end

split_lc(x, τ, ext=true) = begin
    started = false
    a, b = 0, 0
    ranges = Vector{Tuple{Int, Int}}()
    x_ = vcat(x, zeros(eltype(x), 1))
    for i in eachindex(x)
        if started
            if x_[i] <= τ # stop
                b = i - 1
                started = false
                push!(ranges, (a, b))
            end
        else
            if x_[i] > τ # start
                a = i
                b = length(x)
                started = true
            end
        end
    end
    if ext
        ranges = map(r -> (max(1, r[1] - 1), min(length(x), r[2] + 1)), ranges)
    end
    return ranges
end

calc_range(df, M1, τ, ε, k) = begin
    z = df.z[end]
    mz = df.mz[end]
    ms1s = M1[df.file[end]]

    x = map(s -> s.retention_time, values(ms1s))
    ys = map(n -> smooth(map(s -> get_inten(mz + n * Δ / z, s.peaks, ε), values(ms1s)), k), 0:2)
    y = ys[1] .* (ys[1] .> τ) .* (ys[2] .> τ) .* (ys[3] .> τ)
    df_range = DataFrames.DataFrame([(; id, start=x[r[1]], stop=x[r[2]], start_i=r[1], stop_i=r[2]) for (id, r) in enumerate(split_lc(y, τ))])
    return x, ys, y, df_range
end

plot_range(x, y, gd_site, df_site, df_range) = begin
    ls = [scatter(x=x, y=y, mode="lines", name="base LC", line_dash="dash", line_color="grey")]
    for (i, g) in enumerate(gd_site)
        push!(ls, scatter(x=g.rt, y=g.abu, name=df_site.site[i], customdata=g.id, mode="markers", marker_size=8))
    end
    for r in eachrow(df_range)
        push!(ls, scatter(x=x[r.start_i:r.stop_i], y=y[r.start_i:r.stop_i], name="range#$(r.id)", mode="lines"))
    end
    return Plot(ls, Layout(; xaxis_title="retention time (s)", yaxis_title="abundance", clickmode="event+select"))
end

group_site(df) = begin
    pep_a = df.pep_a[end]
    pep_b = df.pep_b[end]

    gd_site = DataFrames.groupby(df, [:site_a, :site_b]; sort=false)
    df_site = DataFrames.unique(df[:, [:site_a, :site_b]])
    df_site.site = ["$(pep_a[r.site_a])-$(pep_b[r.site_b])" for r in eachrow(df_site)]
    df_site.is_kk = [r.site == "K-K" for r in eachrow(df_site)]
    return gd_site, df_site
end

match_range!(df, df_range) = begin
    df_range.n_psm .= 0
    df_range.same_site .= true
    for row in eachrow(df)
        for r in eachrow(df_range)
            if r.start <= row.rt <= r.stop
                row.range_id = r.id
                r.n_psm += 1
            end
        end
    end
    for g in DataFrames.groupby(df, :range_id; sort=false)
        if g.range_id[begin] < 0
            println(g)
            continue
        end
        df_range.same_site[g.range_id[begin]] = all(g.site_a .== g.site_a[begin]) && all(g.site_b .== g.site_b[begin])
    end
    return df_range
end

build_app(gd_grp, df_grp, df_psm, M1, M2D, τ, ε, smooth_k, tab_ele, tab_aa, tab_mod, tab_xl) = begin
    app = dash(; assets_folder=DIR_DATA)
    app.layout = html_div() do
        html_h1("XSiteView", style=Dict("text-align"=>"center")),
        dash_datatable(
            id="group_table",
            style_table=Dict("min-width"=>"100%", "overflow-x"=>"auto"),
            style_cell=Dict("overflow"=>"hidden", "text-overflow"=>"ellipsis", "min-width"=>"64px", "max-width"=>"256px"),
            columns=[(; name=i, id=i) for i in names(df_grp)],
            data=Dict.(pairs.(eachrow(string.(df_grp)))),
            filter_action="native",
            sort_action="native",
            sort_mode="multi",
            row_selectable="single",
            page_action="native",
            page_size=10,
            export_format="csv",
            export_headers="display",
        ),
        dcc_graph(id="group_graph", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"XSiteView_LC"))),
        dcc_tabs() do
            dcc_tab(; label="PSM List") do
                dash_datatable(
                    id="psm_table",
                    style_table=Dict("min-width"=>"100%", "overflow-x"=>"auto"),
                    style_cell=Dict("overflow"=>"hidden", "text-overflow"=>"ellipsis", "min-width"=>"64px", "max-width"=>"256px"),
                    columns=[(; name="psm", id="psm"),],
                    filter_action="native",
                    sort_action="native",
                    sort_mode="multi",
                    row_selectable="single",
                    page_action="native",
                    page_size=10,
                    export_format="csv",
                    export_headers="display",
                ),
                dcc_graph(id="seq_graph", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"XSiteView_SEQ"))),
                dcc_graph(id="psm_graph", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"XSiteView_PSM")))
            end,
            dcc_tab(; label="Fragment Ion") do
                dcc_graph(id="ion_graph", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"XSiteView_ION"))),
                dcc_graph(id="ion_lc_graph", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"XSiteView_ION_LC")))
            end
        end
    end

    callback!(app,
        Output("group_graph", "figure"),
        Input("group_table", "derived_virtual_data"),
        Input("group_table", "derived_virtual_selected_rows"),
    ) do v1, v2
        id = parse(Int, v1[v2[begin] + 1].id)
        df = gd_grp[id]
        x, ys, y, df_range = calc_range(df, M1, τ, ε, smooth_k)
        gd_site, df_site = group_site(df)
        df_range = match_range!(df, df_range)
        return plot_range(x, y, gd_site, df_site, df_range)
    end

    callback!(app,
        Output("psm_table", "columns"),
        Output("psm_table", "data"),
        Output("ion_graph", "figure"),
        Input("group_table", "derived_virtual_data"),
        Input("group_table", "derived_virtual_selected_rows"),
        Input("group_graph", "selectedData"),
    ) do v1, v2, v3
        id = parse(Int, v1[v2[begin] + 1].id)
        df = gd_grp[id]
        if !isnothing(v3)
            selected = map(p -> p.customdata, v3.points)
            s = [id ∈ selected for id in df.id]
            df = df[s, :]
        end
        r = df[1, :]
        site_pairs = [(r.site_a, r.site_b) for r in eachrow(df)]
        ions_a = [(; i..., seq="α") for i in build_ions(r.pep_a, r.mod_a, tab_ele, tab_aa, tab_mod)]
        ions_b = [(; i..., seq="β") for i in build_ions(r.pep_b, r.mod_b, tab_ele, tab_aa, tab_mod)]
        p_ion = seq_crosslink((r.pep_a, r.pep_b), (r.mod_a, r.mod_b), site_pairs, [ions_a, ions_b])
        return [(; name=i, id=i) for i in names(df)], Dict.(pairs.(eachrow(string.(df)))), p_ion
    end

    callback!(app,
        Output("seq_graph", "figure"),
        Output("psm_graph", "figure"),
        Input("psm_table", "derived_virtual_data"),
        Input("psm_table", "derived_virtual_selected_rows"),
    ) do v1, v2
        id = parse(Int, v1[v2[begin] + 1].id)
        r = df_psm[id, :]
        m2 = M2D[r.file][r.scan]
        seqs = (r.pep_a, r.pep_b)
        modss = (r.mod_a, r.mod_b)
        linker = getproperty(tab_xl, Symbol(r.linker))
        sites = (r.site_a, r.site_b)
        ionss = MesMS.build_ions_crosslink(m2.peaks, seqs, modss, linker, sites, ε, tab_ele, tab_aa, tab_mod)
        p_seq = MesMS.Plotly.seq_crosslink(seqs, modss, sites, ionss)
        p_psm = MesMS.Plotly.spec(m2.peaks, filter(i -> i.peak > 0, vcat(ionss...)))
        return p_seq, p_psm
    end

    callback!(app,
        Output("ion_lc_graph", "figure"),
        Input("group_table", "derived_virtual_data"),
        Input("group_table", "derived_virtual_selected_rows"),
        Input("ion_graph", "selectedData"),
    ) do v1, v2, v3
        id = parse(Int, v1[v2[begin] + 1].id)
        df = gd_grp[id]
        ions = isnothing(v3) ? [] : filter(p -> hasproperty(p, :customdata), v3.points)
        ions = map(i -> i.customdata, ions)

        ms2s = [M2D[df.file[end]][r.scan] for r in eachrow(df)]
        ms2s = sort(ms2s, by=s -> s.retention_time)
        linker = getproperty(tab_xl, Symbol(df[1, :linker]))
        xs = map(s -> s.retention_time, ms2s)
        ls = AbstractTrace[]
        for ion in ions
            s, name, mz, z = ion
            for (sym, δ, sty) in zip(["", "'", "''"], [0, linker.masses...], ["solid", "dot", "dash"])
                mz_ = mz + δ / z
                ys = map(s -> get_inten(mz_, s.peaks, ε), values(ms2s))
                push!(ls, scatter(x=xs, y=ys, mode="lines+markers", line_dash=sty, name=@sprintf("%s (%.4f Th)", s * sym * name, mz_)))
            end
        end
        p_ion = Plot(ls, Layout(; yaxis_title="abundance"))

        mz = df.mz[end]
        ms1s = M1[df.file[end]]
        xs = map(s -> s.retention_time, values(ms1s))
        ys = map(s -> get_inten(mz, s.peaks, ε), values(ms1s))
        ls = [scatter(x=xs, y=ys, mode="lines", name=@sprintf("%s (%.4f Th)", "M", mz), line_color="grey")]
        p_m = Plot(ls, Layout(; xaxis_title="retention time (s)", yaxis_title="abundance"))
        p = [p_ion; p_m]
        relayout!(p, Layout(Subplots(rows=2, cols=1, row_heights=[1, 3], shared_xaxes=true)))
        return p
    end
    return app
end

prepare(args) = begin
    path_psm = args["psm"]
    path_ms = args["ms"]
    out = mkpath(args["out"])
    linker = Symbol(args["linker"])
    ε = parse(Float64, args["error"]) * 1.0e-6
    τ = parse(Float64, args["inten"])
    smooth_size = parse(Int, args["smooth"])
    cfg = args["cfg"]
    host = parse(IPAddr, args["host"])
    port = parse(Int, args["port"])
    return (; path_psm, path_ms, out, linker, ε, τ, smooth_size, cfg, host, port)
end

process(; path_psm, path_ms, out, linker,  ε, τ, smooth_size, cfg, host, port) = begin
    smooth_k = vcat(Vector(1:smooth_size), Vector(smooth_size-1:-1:1))
    M = MesMS.read_all(MesMS.read_ms, path_ms)
    M1 = MesMS.mapvalue(m -> m.MS1, M)
    M2 = MesMS.mapvalue(m -> m.MS2, M)
    M1I = map(((k, v),) -> k => (map(m -> m[2].id => m[1], enumerate(v)) |> Dict), collect(M1)) |> Dict
    M2I = map(((k, v),) -> k => (map(m -> m[2].id => m[1], enumerate(v)) |> Dict), collect(M2)) |> Dict
    M1D = map(((k, v),) -> k => (map(m -> m.id => m, v) |> Dict), collect(M1)) |> Dict
    M2D = map(((k, v),) -> k => (map(m -> m.id => m, v) |> Dict), collect(M2)) |> Dict
    df_psm = pLink.read_psm(path_psm)

    df_psm.linker .= linker
    df_psm.group_id .= -1
    df_psm.range_id .= -1
    df_psm.rt = [M2D[r.file][r.scan].retention_time for r in eachrow(df_psm)]
    df_psm.abu .= 0.0

    keys = [:file, :pep_a, :pep_b, :mod_a, :mod_b, :z]
    gd_grp = DataFrames.groupby(df_psm, keys; sort=false)
    df_grp = DataFrames.unique(df_psm[:, keys])

    df_grp.n_psm .= -1
    df_grp.n_site .= -1
    df_grp.n_site_not_kk .= -1
    df_grp.n_range .= -1
    df_grp.n_range_iden .= -1
    df_grp.n_range_diff .= -1
    df_grp.id = Vector(1:size(df_grp, 1))
    DataFrames.select!(df_grp, :id, DataFrames.Not([:id]))

    @showprogress for (id, df) in enumerate(gd_grp)
        df.group_id .= id
        x, ys, y, df_range = calc_range(df, M1, τ, ε, smooth_k)
        df.abu = [y[M1I[r.file][M2D[r.file][r.scan].pre]] for r in eachrow(df)]
        gd_site, df_site = group_site(df)
        df_range = match_range!(df, df_range)
        df_grp.n_psm[id] = size(df, 1)
        df_grp.n_site[id] = size(df_site, 1)
        df_grp.n_site_not_kk[id] = sum(.!df_site.is_kk)
        df_grp.n_range[id] = size(df_range, 1)
        df_grp.n_range_iden[id] = sum(df_range.n_psm .> 0)
        df_grp.n_range_diff[id] = sum(.!df_range.same_site)
    end

    MesMS.safe_save(p -> CSV.write(p, df_grp), joinpath(out, splitext(basename(path_psm))[1] * ".grp.csv"))
    MesMS.safe_save(p -> CSV.write(p, df_psm), joinpath(out, splitext(basename(path_psm))[1] * ".psm.csv"))

    if length(cfg) == 0
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

    @async begin
        sleep(4)
        MesMS.open_url("http://$(host):$(port)")
    end
    app = build_app(gd_grp, df_grp, df_psm, M1, M2D, τ, ε, smooth_k, tab_ele, tab_aa, tab_mod, tab_xl)
    run_server(app, host, port)
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="XSiteView")
    ArgParse.@add_arg_table! settings begin
        "--psm"
            help = "PSM path"
            required = true
        "--ms"
            help = ".mes or .ms1/2 file; .ms2/1 file should be in the same directory for .ms1/2"
            required = true
        "--out"
            help = "output directory"
            default = "./out/"
        "--linker"
            help = "default linker"
            default = "DSSO"
        "--error"
            help = "m/z error"
            metavar = "ppm"
            default = "20.0"
        "--inten"
            help = "intensity threshold"
            metavar = "intensity"
            default = "0.0"
        "--smooth"
            help = "smoothing kernel size"
            metavar = "smooth"
            default = "16"
        "--cfg"
            help = "config directory"
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
    process(; prepare(args)...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

julia_main()::Cint = begin
    main()
    return 0
end

end
