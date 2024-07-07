"""
Exhaustive Search View for Crosslink
"""
module ExhaustiveSearchXLView

"""
- candidate crosslink pair list
- (multiple) mass spectrometry data, e.g., `TMS.raw`.
  - The raw data should be converted into an open-source format such as MS1/MS2. [ThermoRawRead](http://thermorawread.ctarn.io) is recommended.
"""
require = true

"""
Once finished, TargetWizard will start a web server for user to review each target.
([example](../assets/manual/ExhaustiveSearchXLView.pdf))
"""
output = true

"""
# Max. MS1 Mass Error
mass error used to match targets, PSMs, and MS scans.

![Exhaustive Search View for Crosslink](../assets/manual/ExhaustiveSearchXLView.png)
"""
usage = true

"""
# Select Crosslink Pair
![Crosslink Pair List](../assets/manual/ExhaustiveSearchXLView_target.png)

# Select MS File
![MS File List](../assets/manual/ExhaustiveSearchXLView_file.png)

# View MS Acquisition
![MS Acquisition](../assets/manual/ExhaustiveSearchXLView_ms.png)

# View Candidate MS List
![MS List](../assets/manual/ExhaustiveSearchXLView_ms_list.png)

# View Candidate PSM
![Candidate PSM](../assets/manual/ExhaustiveSearchXLView_psm.png)
"""
example = true

using Sockets

import ArgParse
import CSV
import DataFrames
import ProgressMeter: @showprogress
import RelocatableFolders: @path
import UniMZ
import UniMZUtil: pFind, pLink

using Dash
using PlotlyBase

include("../util.jl")

const DIR_DATA = @path joinpath(@__DIR__, "../../data/dash")

Δ = 1.00335

plot_lc(tg, name, dfs_m1, dfs_m2, ε) = begin
    df_m1, df_m2 = dfs_m1[name], dfs_m2[name]
    x = df_m1.rt
    ys = map(n -> map(p -> UniMZ.max_inten_ε(p, tg.mz + n * Δ / tg.z, ε), df_m1.peaks), -1:2)
    ls = [scatter(x=x, y=ys[2], mode="lines", name="M")]
    push!(ls, scatter(x=x, y=ys[3], mode="lines", name="M + 1 Da"))
    push!(ls, scatter(x=x, y=ys[4], mode="lines", name="M + 2 Da"))
    push!(ls, scatter(x=x, y=ys[1], mode="lines", name="M - 1 Da"))
    p1 = Plot(ls, Layout(; xaxis_title="retention time (s)", yaxis_title="abundance"))
    ls = [
        scatter(x=ones(2) * tg.start * 60, y=[-ε, ε], mode="lines", line_dash="dash", line_color="red", name="RT start"),
        scatter(x=ones(2) * tg.stop * 60, y=[-ε, ε], mode="lines", line_dash="dash", line_color="red", name="RT stop"),
    ]
    append!(ls, [scatter(x=df_m2.rt[i:i], y=[UniMZ.error_rel(tg.mz, df_m2.mz[i])], name="MS2#$(df_m2.id[i])", customdata=[i], mode="markers") for i in tg.m2_[name]])
    p2 = Plot(ls, Layout(; yaxis_title="m/z error"))
    p = [p1; p2]
    relayout!(p, Layout(Subplots(rows=2, cols=1, row_heights=[1, 2], vertical_spacing=0.02, shared_xaxes=true), clickmode="event+select"), height=500)
    return p
end

build_app(df, dfs_m1, dfs_m2, tab_ele, tab_aa, tab_mod, tab_xl, ε) = begin
    tb = DataFrames.select(df, filter(n -> !endswith(n, "_"), names(df)))
    tbs_m2 = UniMZ.mapvalue(df -> DataFrames.select(df, DataFrames.Not(:peaks)), dfs_m2)
    app = dash(; assets_folder=DIR_DATA)
    app.layout = html_div() do
        html_h1("Cross-link Exhaustive Search View", style=Dict("text-align"=>"center")),
        dash_datatable(
            id="tg_table",
            style_table=Dict("min-width"=>"100%", "overflow-x"=>"auto"),
            style_cell=Dict("overflow"=>"hidden", "text-overflow"=>"ellipsis", "min-width"=>"64px", "max-width"=>"256px"),
            columns=[(; name=i, id=i) for i in names(tb)],
            data=Dict.(pairs.(eachrow(string.(tb)))),
            filter_action="native",
            sort_action="native",
            sort_mode="multi",
            row_selectable="single",
            page_action="native",
            page_size=10,
            export_format="csv",
            export_headers="display",
        ),
        dcc_dropdown(
            id="run",
            options=[(; label=string(k), value=k) for k in keys(dfs_m1)],
            value=first(keys(dfs_m1)),
        ),
        dcc_graph(id="lc_graph", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"XSeekView_LC"))),
        dash_datatable(
            id="ms_table",
            style_table=Dict("min-width"=>"100%", "overflow-x"=>"auto"),
            style_cell=Dict("overflow"=>"hidden", "min-width"=>"64px",),
            columns=[(; name=i, id=i) for i in names(first(values(tbs_m2)))],
            filter_action="native",
            sort_action="native",
            sort_mode="multi",
            row_selectable="single",
            page_action="native",
            page_size=10,
            export_format="csv",
            export_headers="display",
        ),
        dcc_graph(id="seq_graph", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"XSeekView_SEQ"))),
        dcc_graph(id="psm_graph", config=PlotConfig(toImageButtonOptions=Dict(:format=>"svg", :filename=>"XSeekView_PSM")))
    end

    callback!(app,
        Output("lc_graph", "figure"),
        Input("tg_table", "derived_virtual_data"),
        Input("tg_table", "derived_virtual_selected_rows"),
        Input("run", "value"),
    ) do v1, v2, run
        run = Symbol(run)
        id = parse(Int, v1[v2[begin] + 1].id)
        tg = df[id, :]
        fig = plot_lc(tg, run, dfs_m1, dfs_m2, ε)
        return fig
    end

    callback!(app,
        Output("ms_table", "data"),
        Input("lc_graph", "selectedData"),
        Input("run", "value"),
    ) do v, run
        run = Symbol(run)
        ids = map(p -> p.customdata, v.points)
        return Dict.(pairs.(eachrow(string.(tbs_m2[run][ids, :]))))
    end

    callback!(app,
        Output("seq_graph", "figure"),
        Output("psm_graph", "figure"),
        Input("tg_table", "derived_virtual_data"),
        Input("tg_table", "derived_virtual_selected_rows"),
        Input("ms_table", "derived_virtual_data"),
        Input("ms_table", "derived_virtual_selected_rows"),
        Input("run", "value"),
    ) do v1, v2, v3, v4, run
        run = Symbol(run)
        id = parse(Int, v1[v2[begin] + 1].id)
        r = df[id, :]
        idx = parse(Int, v3[v4[begin] + 1].idx)
        m2 = dfs_m2[run][idx, :]
        seqs = (r.pep_a, r.pep_b)
        modss = (r.mod_a, r.mod_b)
        linker = tab_xl[Symbol(r.linker)]
        sites = (r.site_a, r.site_b)
        ionss = UniMZ.build_ions_crosslink(m2.peaks, seqs, modss, linker, sites, ε, tab_ele, tab_aa, tab_mod)
        p_seq = UniMZ.Plotly.seq_crosslink(seqs, modss, sites, ionss)
        p_psm = UniMZ.Plotly.spec(m2.peaks, filter(i -> i.peak > 0, vcat(ionss...)))
        return p_seq, p_psm
    end
    return app
end

pepmass(pep, mod, tab_ele, tab_aa, tab_mod) = begin
    H2O = UniMZ.mass(UniMZ.Formula(H=2, O=1), tab_ele)
    return sum(x -> tab_aa[Symbol(x)], pep) + sum(x -> tab_mod[x.id], mod; init=0.0) + H2O
end

prepare(args) = begin
    out = mkpath(args["out"])
    linker = Symbol(args["linker"])
    ε = parse(Float64, args["error"]) * 1.0e-6
    tab_ele, tab_aa, tab_mod, tab_xl = pLink.read_mass_table(args["cfg"])
    host = parse(IPAddr, args["host"])
    port = parse(Int, args["port"])
    return (; out, linker, ε, tab_ele, tab_aa, tab_mod, tab_xl, host, port)
end

process(path, path_ms; out, linker, ε, tab_ele, tab_aa, tab_mod, tab_xl, host, port) = begin
    Ms = map(UniMZ.read_ms, path_ms)
    runs = [splitext(basename(p))[1] for p in path_ms] .|> Symbol
    dfs_m1 = map(Ms, runs) do M, p
        p => [(; idx, m.id, rt=m.retention_time, m.peaks) for (idx, m) in enumerate(M.MS1)] |> DataFrames.DataFrame
    end |> Dict
    dfs_m2 = map(Ms, runs) do M, p
        p => [(; idx, m.id, mz=m.activation_center, rt=m.retention_time, m.peaks) for (idx, m) in enumerate(M.MS2)] |> DataFrames.DataFrame
    end |> Dict

    df = DataFrames.DataFrame(CSV.File(path)) |> unique
    df.id = Vector(1:size(df, 1))
    src = [:peptide, :modification]
    dst = [:pep_a, :mod_a, :site_a, :prot_a, :pep_b, :mod_b, :site_b, :prot_b]
    DataFrames.transform!(df, src => DataFrames.ByRow(pLink.parse_pair) => dst)

    df.linker .= linker
    DataFrames.rename!(df, :charge => :z)
    df.m = map(eachrow(df)) do r
        pepmass(r.pep_a, r.mod_a, tab_ele, tab_aa, tab_mod) + pepmass(r.pep_b, r.mod_b, tab_ele, tab_aa, tab_mod) + tab_xl[r.linker].mass
    end
    df.mz = UniMZ.m_to_mz.(df.m, df.z)

    df.m = round.(df.m; digits=4)
    df.mz = round.(df.mz; digits=4)

    DataFrames.select!(df, :id, :m, :mz, :z, dst..., :linker, DataFrames.Not(src))
    DataFrames.select!(df, DataFrames.Not([:prot_a, :prot_b]))
    ("start" ∉ names(df)) && (df.start .= 0.0)
    ("stop" ∉ names(df)) && (df.stop .= 0.0)

    @info "MS2 mapping"
    # TODO: detect precursor ions
    df.m2_ = [Dict{eltype(runs), Vector{Int}}() for _ in eachrow(df)]
    for k in runs
        df_m2 = dfs_m2[k]
        tmp = sort!([(x.mz::Float64, x.idx::Int) for x in eachrow(df_m2)])
        mzs = map(x -> x[1], tmp)
        ids = map(x -> x[2], tmp)
        for r in eachrow(df)
            r.m2_[k] = sort(ids[UniMZ.argquery_ε(mzs, r.mz, ε)])
        end
    end

    @async begin
        sleep(4)
        UniMZ.open_url("http://$(host):$(port)")
    end
    app = build_app(df, dfs_m1, dfs_m2, tab_ele, tab_aa, tab_mod, tab_xl, ε)
    run_server(app, host, port)
end

main() = begin
    settings = ArgParse.ArgParseSettings(prog="XSeekView")
    ArgParse.@add_arg_table! settings begin
        "list"
            help = "target list"
            required = true
        "--ms"
            help = ".umz or .ms1/2 file; .ms2/1 file should be in the same directory for .ms1/2"
            required = true
            nargs = '+'
        "--out"
            help = "output directory"
            default = "./out/"
            metavar = "./out/"
        "--linker"
            help = "default linker"
            default = "DSSO"
        "--error"
            help = "m/z error"
            metavar = "ppm"
            default = "20.0"
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
    paths = reduce(vcat, UniMZ.match_path.(args["ms"], ".umz")) |> unique |> sort
    @info "file paths of selected MS data:"
    foreach(x -> println("$(x[1]):\t$(x[2])"), enumerate(paths))
    process(args["list"], paths; prepare(args)...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

julia_main()::Cint = begin
    main()
    return 0
end

end
