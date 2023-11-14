using Statistics

import DataFrames
import MesMS
import ProgressMeter: @showprogress

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
    df.m = MesMS.mz_to_m.(df.mz, df.z)
    df.rt = (df.start .+ df.stop) ./ 2
    return df
end

calc_cov_crosslink!(df, M, ε, ion_syms, ion_types, tab_ele, tab_aa, tab_mod, tab_xl) = begin
    @info "cross-linked peptide processing"
    @info "fragment ion calculating"
    df.ion = @showprogress map(eachrow(df)) do r
        peaks = M[r.file][r.scan].peaks
        seqs = (r.pep_a, r.pep_b)
        modss = (r.mod_a, r.mod_b)
        sites = (r.site_a, r.site_b)
        types = [(i, 1:(r.z-1)) for i in ion_types]
        ionss = MesMS.build_ions_crosslink(peaks, seqs, modss, tab_xl[r.linker], sites, ε, tab_ele, tab_aa, tab_mod; types)
        return filter.(i -> i.peak > 0 && i.loc > 0, ionss)
    end

    @info "coverage calculating"
    df.ion_a = first.(df.ion)
    df.ion_b = last.(df.ion)
    foreach(r -> filter!(i -> i.loc < length(r.pep_a), r.ion_a), eachrow(df))
    foreach(r -> filter!(i -> i.loc < length(r.pep_b), r.ion_b), eachrow(df))

    match_a = [falses(length(r.pep_a)-1) for r in eachrow(df)]
    match_b = [falses(length(r.pep_b)-1) for r in eachrow(df)]
    for idx in 1:size(df, 1)
        foreach(i -> match_a[idx][i.loc] = true, df.ion_a[idx])
        foreach(i -> match_b[idx][i.loc] = true, df.ion_b[idx])
    end

    df.cov = round.(mean.(vcat.(match_a, match_b)); digits=4)
    df.cov_a = round.(mean.(match_a); digits=4)
    df.cov_b = round.(mean.(match_b); digits=4)

    @info "coverage of each type of fragment ion calculating"
    @showprogress for (sym, ion_type) in zip(ion_syms, ion_types)
        ion_a = map(r -> filter(i -> i.type == ion_type.type, r.ion_a), eachrow(df))
        ion_b = map(r -> filter(i -> i.type == ion_type.type, r.ion_b), eachrow(df))
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

    df.prot_a = pFind.protstr.(df.prot_a)
    df.prot_b = pFind.protstr.(df.prot_b)

    DataFrames.select!(df, DataFrames.Not([:ion, :ion_a, :ion_b]), :ion_a, :ion_b)
    df.ion_a = map(ions -> join(getfield.(ions, :text_abbr), ','), df.ion_a)
    df.ion_b = map(ions -> join(getfield.(ions, :text_abbr), ','), df.ion_b)

    return df
end

calc_cov_linear!(df, M, ε, ion_syms, ion_types, tab_ele, tab_aa, tab_mod) = begin
    @info "linear peptide processing"
    @info "fragment ion calculating"
    df.ion = @showprogress map(eachrow(df)) do r
        peaks = M[r.file][r.scan].peaks
        types = [(i, 1:(r.z-1)) for i in ion_types]
        ions = MesMS.build_ions(peaks, r.pep, r.mod, ε, tab_ele, tab_aa, tab_mod; types)
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
        ions = MesMS.build_ions_monolink(peaks, r.pep, r.mod, tab_xl[r.linker], r.site, ε, tab_ele, tab_aa, tab_mod; types)
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
        ions = MesMS.build_ions_looplink(peaks, r.pep, r.mod, tab_xl[r.linker], r.site_a, r.site_b, ε, tab_ele, tab_aa, tab_mod; types)
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
