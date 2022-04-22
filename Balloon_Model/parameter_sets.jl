@with_kw struct ballonModelParameters{R}
    E_0::R = 0.34
    κ::R = 0.65
    γ::R = 0.41
    τ::R = 0.098
    α::R = 0.32
    V_0::R = 0.02
    TE::R = 0.03
    v_0::R = 84.795
    ϵ::R = 0.47
    r_0::R = 110.0
    k_1::R = 4.3 *v_0 * E_0 * TE
    k_2::R = ϵ * r_0 * E_0 * TE
    k_3::R = 1.0 - ϵ
end

@with_kw struct ballonModelParameters2{R}
    E_0::R = 0.34
    κ::R = 0.65
    γ::R = 0.41
    τ::R = 2.0
    α::R = 0.32
    V_0::R = 0.02
    TE::R = 0.03
    v_0::R = 84.795
    ϵ::R = 0.47
    r_0::R = 110.0
    k_1::R = 4.3 *v_0 * E_0 * TE
    k_2::R = ϵ * r_0 * E_0 * TE
    k_3::R = 1.0 - ϵ
end