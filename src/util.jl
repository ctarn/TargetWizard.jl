import DataFrames
import MesMS

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
            @error "failed to detect list format"
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
