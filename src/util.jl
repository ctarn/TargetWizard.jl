using Statistics

import DataFrames
import ProgressMeter: @showprogress
import UniMS

parse_target_list!(df, fmt) = begin
    if fmt == :auto
        cols_TW = ["mz", "z", "start", "stop"]
        cols_TmQE = ["Mass [m/z]", "CS [z]", "Start [min]", "End [min]"]
        cols_TmFu = ["m/z", "z", "t start (min)", "t stop (min)"]
        if all(n -> n ∈ names(df), cols_TW)
            @info "list treated as `TargetWizard` format based on following columns: $(cols_TW)"
            fmt = :TW
        elseif all(n -> n ∈ names(df), cols_TmQE)
            @info "list treated as `Thermo Q Exactive` format based on following columns: $(cols_TmQE)"
            fmt = :TmQE
        elseif all(n -> n ∈ names(df), cols_TmFu)
            @info "list treated as `Thermo Fusion` format based on following columns: $(cols_TmFu)"
            fmt = :TmFu
        else
            error("failed to detect list format")
            fmt = :unknown
        end
    end
    if fmt == :TmQE
        DataFrames.rename!(df, "Mass [m/z]" => "mz", "CS [z]" => "z", "Start [min]" => "start", "End [min]" => "stop")
        df.start = df.start .* 60
        df.stop = df.stop .* 60
    elseif fmt == :TmFu
        DataFrames.rename!(df, "m/z" => "mz", "t start (min)" => "start", "t stop (min)" => "stop")
        df.start = df.start .* 60
        df.stop = df.stop .* 60
    end
    df.m = UniMS.mz_to_m.(df.mz, df.z)
    df.rt = (df.start .+ df.stop) ./ 2
    return df
end

calc_cov_crosslink!(df, M, ε, ion_syms, ion_types, tab_ele, tab_aa, tab_mod, tab_xl) = begin
    @info "cross-linked peptide processing"
    @info "fragment ion calculating"
    df.ion_ = @showprogress map(eachrow(df)) do r
        peaks = M[r.file][r.scan].peaks
        seqs = (r.pep_a, r.pep_b)
        modss = (r.mod_a, r.mod_b)
        sites = (r.site_a, r.site_b)
        types = [(i, 1:(r.z-1)) for i in ion_types]
        ionss = UniMS.build_ions_crosslink(peaks, seqs, modss, tab_xl[r.linker], sites, ε, tab_ele, tab_aa, tab_mod; types)
        return filter.(i -> i.peak > 0 && i.loc > 0, ionss)
    end

    @info "coverage calculating"
    df.ion_a_ = first.(df.ion_)
    df.ion_b_ = last.(df.ion_)
    foreach(r -> filter!(i -> i.loc < length(r.pep_a), r.ion_a_), eachrow(df))
    foreach(r -> filter!(i -> i.loc < length(r.pep_b), r.ion_b_), eachrow(df))

    match_a = [falses(length(r.pep_a)-1) for r in eachrow(df)]
    match_b = [falses(length(r.pep_b)-1) for r in eachrow(df)]
    for idx in 1:size(df, 1)
        foreach(i -> match_a[idx][i.loc] = true, df.ion_a_[idx])
        foreach(i -> match_b[idx][i.loc] = true, df.ion_b_[idx])
    end

    df.cov = round.(mean.(vcat.(match_a, match_b)); digits=4)
    df.cov_a = round.(mean.(match_a); digits=4)
    df.cov_b = round.(mean.(match_b); digits=4)

    @info "coverage of each type of fragment ion calculating"
    @showprogress for (sym, ion_type) in zip(ion_syms, ion_types)
        ion_a = map(r -> filter(i -> i.type == ion_type.type, r.ion_a_), eachrow(df))
        ion_b = map(r -> filter(i -> i.type == ion_type.type, r.ion_b_), eachrow(df))
        match_a = [falses(length(r.pep_a)-1) for r in eachrow(df)]
        match_b = [falses(length(r.pep_b)-1) for r in eachrow(df)]
        for idx in 1:size(df, 1)
            foreach(i -> match_a[idx][i.loc] = true, ion_a[idx])
            foreach(i -> match_b[idx][i.loc] = true, ion_b[idx])
        end
        df[!, "cov_ion_$(sym)"] = round.(mean.(vcat.(match_a, match_b)); digits=4)
        df[!, "cov_a_ion_$(sym)"] = round.(mean.(match_a); digits=4)
        df[!, "cov_b_ion_$(sym)"] = round.(mean.(match_b); digits=4)
    end

    df.ion_a = map(ions -> join(getfield.(ions, :text_abbr), ','), df.ion_a_)
    df.ion_b = map(ions -> join(getfield.(ions, :text_abbr), ','), df.ion_b_)
    return df
end

calc_cov_linear!(df, M, ε, ion_syms, ion_types, tab_ele, tab_aa, tab_mod) = begin
    @info "linear peptide processing"
    @info "fragment ion calculating"
    df.ion = @showprogress map(eachrow(df)) do r
        peaks = M[r.file][r.scan].peaks
        types = [(i, 1:(r.z-1)) for i in ion_types]
        ions = UniMS.build_ions(peaks, r.pep, r.mod, ε, tab_ele, tab_aa, tab_mod; types)
        return filter(i -> i.peak > 0 && 0 < i.loc < length(r.pep), ions)
    end

    @info "coverage calculating"
    match = [falses(length(r.pep)-1) for r in eachrow(df)]
    for idx in 1:size(df, 1)
        foreach(i -> match[idx][i.loc] = true, df.ion[idx])
    end

    df.cov = round.(mean.(match); digits=4)

    @info "coverage of each type of fragment ion calculating"
    @showprogress for (sym, ion_type) in zip(ion_syms, ion_types)
        ion = map(r -> filter(i -> i.type == ion_type.type, r.ion), eachrow(df))
        match = [falses(length(r.pep)-1) for r in eachrow(df)]
        for idx in 1:size(df, 1)
            foreach(i -> match[idx][i.loc] = true, ion[idx])
        end
        df[!, "cov_ion_$(sym)"] = round.(mean.(match); digits=4)
    end

    DataFrames.select!(df, DataFrames.Not([:ion]), :ion)
    df.ion = map(ions -> join(getfield.(ions, :text_abbr), ','), df.ion)

    return df
end

calc_cov_monolink!(df, M, ε, ion_syms, ion_types, tab_ele, tab_aa, tab_mod, tab_xl) = begin
    @info "mono-linked peptide processing"
    @info "fragment ion calculating"
    df.ion = @showprogress map(eachrow(df)) do r
        peaks = M[r.file][r.scan].peaks
        types = [(i, 1:(r.z-1)) for i in ion_types]
        ions = UniMS.build_ions_monolink(peaks, r.pep, r.mod, tab_xl[r.linker], r.site, ε, tab_ele, tab_aa, tab_mod; types)
        return filter(i -> i.peak > 0 && 0 < i.loc < length(r.pep), ions)
    end

    @info "coverage calculating"
    match = [falses(length(r.pep)-1) for r in eachrow(df)]
    for idx in 1:size(df, 1)
        foreach(i -> match[idx][i.loc] = true, df.ion[idx])
    end

    df.cov = round.(mean.(match); digits=4)

    @info "coverage of each type of fragment ion calculating"
    @showprogress for (sym, ion_type) in zip(ion_syms, ion_types)
        ion = map(r -> filter(i -> i.type == ion_type.type, r.ion), eachrow(df))
        match = [falses(length(r.pep)-1) for r in eachrow(df)]
        for idx in 1:size(df, 1)
            foreach(i -> match[idx][i.loc] = true, ion[idx])
        end
        df[!, "cov_ion_$(sym)"] = round.(mean.(match); digits=4)
    end

    DataFrames.select!(df, DataFrames.Not([:ion]), :ion)
    df.ion = map(ions -> join(getfield.(ions, :text_abbr), ','), df.ion)

    return df
end

calc_cov_looplink!(df, M, ε, ion_syms, ion_types, tab_ele, tab_aa, tab_mod, tab_xl) = begin
    @info "loop-linked peptide processing"
    @info "fragment ion calculating"
    df.ion = @showprogress map(eachrow(df)) do r
        peaks = M[r.file][r.scan].peaks
        types = [(i, 1:(r.z-1)) for i in ion_types]
        ions = UniMS.build_ions_looplink(peaks, r.pep, r.mod, tab_xl[r.linker], r.site_a, r.site_b, ε, tab_ele, tab_aa, tab_mod; types)
        return filter(i -> i.peak > 0 && 0 < i.loc < length(r.pep), ions)
    end

    @info "coverage calculating"
    match = [falses(length(r.pep)-1) for r in eachrow(df)]
    for idx in 1:size(df, 1)
        foreach(i -> match[idx][i.loc] = true, df.ion[idx])
    end

    df.cov = round.(mean.(match); digits=4)

    @info "coverage of each type of fragment ion calculating"
    @showprogress for (sym, ion_type) in zip(ion_syms, ion_types)
        ion = map(r -> filter(i -> i.type == ion_type.type, r.ion), eachrow(df))
        match = [falses(length(r.pep)-1) for r in eachrow(df)]
        for idx in 1:size(df, 1)
            foreach(i -> match[idx][i.loc] = true, ion[idx])
        end
        df[!, "cov_ion_$(sym)"] = round.(mean.(match); digits=4)
    end

    DataFrames.select!(df, DataFrames.Not([:ion]), :ion)
    df.ion = map(ions -> join(getfield.(ions, :text_abbr), ','), df.ion)

    return df
end

_nearbymax(M, i, mz, ε, τ, δ) = begin
    skip = 0
    v = -Inf
    i_max = i
    while 1 ≤ (i + δ) ≤ length(M) && skip ≤ τ
        v_ = UniMS.max_inten_ε(M[i].peaks, mz, ε)
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
psmstr_linear(x) = "$(x.scan)($(round(x.cov; digits=2))):$(x.pep)($(x.mod))"
psmstr_mono(x) = "$(x.scan)($(round(x.cov; digits=2))):$(x.pep)($(x.mod))@$(x.site)"
psmstr_loop(x) = "$(x.scan)($(round(x.cov; digits=2))):$(x.pep)($(x.mod))@$(x.site_a)-$(x.site_b)"
psmstr_link(x) = "$(x.scan)($(round(x.cov_a; digits=2))|$(round(x.cov_b; digits=2))):$(x.pep_a)($(x.mod_a))@$(x.site_a)-$(x.pep_b)($(x.mod_b))@$(x.site_b)"
is_same_xl(a, b) = (a.pep_a == b.pep_a) && (a.pep_b == b.pep_b) && (a.mod_a == b.mod_a) && (a.mod_b == b.mod_b) && (a.site_a == b.site_a) && (a.site_b == b.site_b)
is_same_xl_pepmod(a, b) = (a.pep_a == b.pep_a) && (a.pep_b == b.pep_b) && (a.mod_a == b.mod_a) && (a.mod_b == b.mod_b)
