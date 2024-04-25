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
