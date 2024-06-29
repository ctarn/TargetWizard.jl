import UniMZ

_nearbymax(M, i, mz, ε, τ, δ) = begin
    skip = 0
    v = -Inf
    i_max = i
    while 1 ≤ (i + δ) ≤ length(M) && skip ≤ τ
        v_ = UniMZ.max_inten_ε(M[i].peaks, mz, ε)
        if v_ > v
            skip = 0
            v = v_
            i_max = i
        else
            skip += 1
        end
        i += δ
    end
    return i_max, v
end

nearbymax(M, i, mz, ε, τ=2) = begin
    li, lv = _nearbymax(M, i, mz, ε, τ, -1)
    ri, rv = _nearbymax(M, i, mz, ε, τ, 1)
    return lv ≥ rv ? (li, lv) : (ri, rv)
end

unify_mods_str(s) = (!ismissing(s) && startswith(s, "Any[") && endswith(s, "]")) ? s[5:end-1] : s
psmstr(x) = "$(x.scan)($(round(x.cov; digits=2))):$(x.pep)($(x.mod))"
psmstr_mono(x) = "$(x.scan)($(round(x.cov; digits=2))):$(x.pep)($(x.mod))@$(x.site)"
psmstr_loop(x) = "$(x.scan)($(round(x.cov; digits=2))):$(x.pep)($(x.mod))@$(x.site_a)-$(x.site_b)"
psmstr_link(x) = "$(x.scan)($(round(x.cov_a; digits=2))|$(round(x.cov_b; digits=2))):$(x.pep_a)($(x.mod_a))@$(x.site_a)-$(x.pep_b)($(x.mod_b))@$(x.site_b)"
is_same_pepmod(a, b) = (a.pep == b.pep) && (a.mod == b.mod)
is_same_xl(a, b) = (a.pep_a == b.pep_a) && (a.pep_b == b.pep_b) && (a.mod_a == b.mod_a) && (a.mod_b == b.mod_b) && (a.site_a == b.site_a) && (a.site_b == b.site_b)
is_same_xl_pepmod(a, b) = (a.pep_a == b.pep_a) && (a.pep_b == b.pep_b) && (a.mod_a == b.mod_a) && (a.mod_b == b.mod_b)

build_target(df, Ms, paths, out, name, ε, batch_size, rt, lc, fmt) = begin
    M1 = map((p, M) -> splitext(basename(p))[1] => UniMZ.dict_by_id(M.MS1), paths, Ms) |> Dict
    M2 = map((p, M) -> splitext(basename(p))[1] => UniMZ.dict_by_id(M.MS2), paths, Ms) |> Dict
    M1V = map((p, M) -> splitext(basename(p))[1] => M.MS1, paths, Ms) |> Dict
    M1I = map((p, M) -> splitext(basename(p))[1] => [m.id => i for (i, m) in enumerate(M.MS1)] |> Dict, paths, Ms) |> Dict

    df.rt = [M2[r.file][r.scan].retention_time for r in eachrow(df)]

    df.inten = map(eachrow(df)) do r
        m2 = M2[r.file][r.scan]
        m1 = M1[r.file][m2.pre]
        return UniMZ.max_inten_ε(m1.peaks, r.mz, ε)
    end
    vs = map(eachrow(df)) do r
        m2 = M2[r.file][r.scan]
        i = M1I[r.file][m2.pre]
        i, v = nearbymax(M1V[r.file], i, r.mz, ε, 2)
        return v, (M1V[r.file][i].retention_time - m2.retention_time)
    end
    df.inten_max = first.(vs)
    df.inten_max_delta_rt = last.(vs)

    df.start = min.(lc * 60, max.(0, df.rt .- (rt / 2)))
    df.stop = min.(lc * 60, max.(0, df.rt .+ (rt / 2)))

    n_batch = isinf(batch_size) ? 1 : ceil(Int, size(df, 1) / batch_size)
    df = sort(df, :rt)
    df.id = 1:size(df, 1)
    df.batch = (df.id .- 1) .% n_batch .+ 1

    @info "$(size(df, 1)) features splitting into $(n_batch) batches"

    tw = :TW ∈ fmt
    tmqe = :TmQE ∈ fmt
    tmfu = :TmFu ∈ fmt
    p = joinpath(out, name)
    tw && UniMZ.safe_save(p -> CSV.write(p, df), "$(p).all.TW.target.csv", "list")

    for i in 1:n_batch
        df_ = df[df.batch .== i, :]
        @info "batch $(i): $(size(df_, 1))"
        tw && UniMZ.safe_save(p -> CSV.write(p, df_), "$(p).batch$(i).TW.target.csv", "list")
        tmqe && UniMZ.safe_save(p -> CSV.write(p, TMS.build_target_TmQE(df_)), "$(p).batch$(i).TmQE.target.csv", "list (Thermo Q Exactive)")
        tmfu && UniMZ.safe_save(p -> CSV.write(p, TMS.build_target_TmFu(df_)), "$(p).batch$(i).TmFu.target.csv", "list (Thermo Fusion)")
    end
end
