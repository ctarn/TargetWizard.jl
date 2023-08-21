module XSiteView

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
        dcc_graph(id="group_graph"),
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
        dcc_graph(id="seq_graph"),
        dcc_graph(id="psm_graph")
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
        return [(; name=i, id=i) for i in names(df)], Dict.(pairs.(eachrow(string.(df))))
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
        ionss = MesMS.Plot.build_ions_xl(m2.peaks, seqs, modss, linker, sites, ε, tab_ele, tab_aa, tab_mod)
        p_seq = MesMS.Plotly.seq_xl(seqs, modss, sites, ionss)
        p_psm = MesMS.Plotly.spec(m2.peaks, filter(i -> i.peak > 0, vcat(ionss...)))
        return p_seq, p_psm
    end
    return app
end

prepare(args) = begin
    path_psm = args["psm"]
    path_ms = args["ms"]
    out = mkpath(args["out"])
    ε = parse(Float64, args["error"]) * 1.0e-6
    τ = parse(Float64, args["inten"])
    smooth_size = parse(Int, args["smooth"])
    cfg = args["cfg"]
    host = parse(IPAddr, args["host"])
    port = parse(Int, args["port"])
    return (; path_psm, path_ms, out, ε, τ, smooth_size, cfg, host, port)
end

process(; path_psm, path_ms, out, ε, τ, smooth_size, cfg, host, port) = begin
    smooth_k = vcat(Vector(1:smooth_size), Vector(smooth_size-1:-1:1))
    M = MesMS.read_all(MesMS.read_ms, path_ms)
    M1 = MesMS.mapvalue(m -> m.MS1, M)
    M2 = MesMS.mapvalue(m -> m.MS2, M)
    M1I = map(((k, v),) -> k => (map(m -> m[2].id => m[1], enumerate(v)) |> Dict), collect(M1)) |> Dict
    M2I = map(((k, v),) -> k => (map(m -> m[2].id => m[1], enumerate(v)) |> Dict), collect(M2)) |> Dict
    M1D = map(((k, v),) -> k => (map(m -> m.id => m, v) |> Dict), collect(M1)) |> Dict
    M2D = map(((k, v),) -> k => (map(m -> m.id => m, v) |> Dict), collect(M2)) |> Dict
    df_psm = pLink.read_psm(path_psm)

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
        tab_mod = MesMS.mapvalue(x -> x.mass, pLink.read_mod())
        tab_xl = pLink.read_linker() |> NamedTuple
    else
        tab_ele = pLink.read_element(joinpath(cfg, "element.ini")) |> NamedTuple
        tab_aa = map(x -> MesMS.mass(x, tab_ele), pLink.read_amino_acid(joinpath(cfg, "aa.ini")) |> NamedTuple)
        tab_mod = MesMS.mapvalue(x -> x.mass, pLink.read_mod(joinpath(cfg, "modification.ini")))
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
