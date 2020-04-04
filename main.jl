# This code generates figures 1, 2, 3, 4 and D.1 from 
# Aguiar, Mark and Manuel Amador, "Self-fulfilling Debt Dilution: Maturity and 
# Multiplicity in Sovereign Debt Models". 
#
# The code uses the analytical closed-form solutions in the paper 
# to make the figures. 
#
# Julia 1.4 


using Parameters
using PGFPlotsX
using Roots

# basic parameter values
Model = @with_kw (
    r = 1.0,
    ρ = 2.0, 
    λ = 2.0, 
    y = 1.0, 
    v_low = 0.8 / 2.0, # low outside option
    v_high = 0.95 / 2.0,  # high outside option
    c_low = 0.0, # minimum consumption 
    c_high = 1.2, # maximum consumption 
    δ = 10.0 
)


function get_pstar_s(m)
    @unpack r, λ, δ, y, ρ, v_high, v_low, c_high, c_low = m
    p_star_safe = function (v)
        (y - c_high + (c_high - ρ * v_high) ^ ((ρ - r) / ρ) * 
                (c_high - ρ * v)^(r / ρ)
        ) / r
    end
    p_vbar = (y - ρ * v_high) / r 
    p_hat = function (v)
        (
            y - c_low + (c_low - y + (r + λ) * p_vbar) * 
            (
                (c_low + λ * v_high - (ρ + λ) * v) /
                (c_low + λ * v_high - (ρ + λ) * v_high) 
            ) ^ ((r + λ) / (ρ + λ))
        ) / (r + λ)
    end 
    p_b = function (v)
        y - c_high + 
        (c_high + λ * v_high - (ρ + λ) * v_low)^((ρ - r) / (ρ + λ)) * 
        (c_high + λ * v_high - (ρ + λ) * v) ^ ((r + λ) / (ρ + λ))
    end
    if p_hat(v_low) < p_b(v_low) && p_hat(v_high) > p_b(v_high)
        v_i = find_zero(v -> p_b(v) - p_hat(v), (v_low, v_high))
    else 
        v_i = v_low
    end 
    p_s = function (v)
        if v ≥ v_high
            p_star_safe(v)
        elseif v_i ≤ v ≤ v_high
            p_hat(v)
        else 
            p_b(v)
        end 
    end
    return (p_s = p_s)
end 


function get_pstar_b(m)
    @unpack r, λ, δ, y, ρ, v_high, v_low, c_high, c_low = m
    p_crisis = function (v) 
        (1 / (r + λ)) * (
            y - c_high + 
            (c_high + λ * v_high - (ρ + λ) * v_low) ^ ((ρ - r) / (ρ + λ)) * 
            (c_high + λ * v_high - (ρ + λ) * v) ^ ( (r + λ) / (ρ + λ))
        )
    end
    p_vbar = p_crisis(v_high)
    p_star_safe = function (v)
        (1/r) * (
            y - c_high + 
            (c_high - y + r * p_vbar) * 
            (
                (c_high - ρ * v) / 
                (c_high - ρ * v_high)
            ) ^ (r / ρ)
        )
    end 
    p_b = function (v) 
        if v ≥ v_high 
            p_star_safe(v)
        else 
            p_crisis(v)
        end 
    end 
    return (p_b = p_b)
end 


function get_savings_eqm(m, q_low, b_max_b, b_b, b_s)
    
    @unpack r, λ, δ, y, ρ, v_high, v_low, 
    c_high, c_low = m

    v_safe = function (b)
        (1/ρ) * (
            c_high - 
            (c_high + r * b - y) ^ (ρ / r) * 
            (c_high - ρ * v_high) ^ ( (r - ρ) / r)
        )
    end 
    q_hat = function (b) 
        (
            r + δ + 
            (λ + ρ - r) * (b_s / b) ^ ( (ρ + δ + λ) / δ)
        ) / 
        (ρ + δ + λ)
    end
    c_hat = function (b)
        y - (r + δ * (1 - q_hat(b))) * b - 
        q_hat(b) * δ * b * (
            (q_hat(b) - q_low) / 
            (q_hat(b) - q_low + (ρ - r) * q_hat(b) / (r + δ + λ))
        )
    end
    v_save = function (b) 
        (y - (r + δ  *(1 - q_hat(b))) * b + λ * v_high) / (ρ + λ)
    end
    v_borrow = function (b)
        (1 / (λ + ρ)) * (
            c_high + λ * v_high - 
            (c_high + λ * v_high - (λ + ρ) * v_low) ^ ( (r - ρ) / (r + λ)) *
            (c_high - y + (r + λ) * q_low * b) ^ ( (ρ + λ) / (r + λ))
        )
    end 
    if b_b > b_s 
        @warn "No good equilibrium exists"
    end 
    b_max_s = find_zero(b -> v_save(b) - v_low, b_s)

    if b_max_s < b_max_b && b_s ≥ b_b 
        b_i = find_zero(b -> v_borrow(b) - v_save(b), (b_s, b_max_s))
        b_max_s = b_max_b
    else 
        b_i = b_max_s
    end 

    v_s = function (b) 
        if b ≤ b_s 
            v_safe(b)
        elseif b_i ≥ b ≥ b_s 
            v_save(b)
        else 
            v_borrow(b)
        end 
    end 
    q_s = function (b) 
        if b ≤ b_s 
            1.0
        elseif b_i ≥ b ≥ b_s 
            q_hat(b)
        else 
            q_low
        end 
    end
    c_s = function (b) 
        if b < b_s 
            c_high
        elseif b == b_s 
            y - r * b_s 
        elseif b_i ≥ b > b_s 
            c_hat(b)
        elseif b_max_s > b > b_i 
            c_high
        else
            y - (r + δ * (1 - q_low)) * b_max_s
        end 
    end

    return (
        c_s = c_s, q_s = q_s, v_s = v_s, b_i = b_i,
        b_max_s = b_max_s
    )
end 


function get_borrowing_eqm(m, q_low, b_max_b, b_b, b_s)
    @unpack r, λ, δ, y, ρ, v_high, v_low, c_high, 
            c_low = m

    v_crisis = function (b)
        (1 / (ρ + λ))  * (
            c_high + λ * v_high - 
            (c_high + λ * v_high - (λ + ρ) * v_low) ^ ((r - ρ) / (r + λ)) * 
            (c_high - y + (r + λ) * q_low * b) ^ ( (ρ + λ) / (r + λ))
        )
    end 
    q_safe(b) = find_zero( (q_low, 1.0) ) do (q)
        c_high + r * q * b - y - (c_high + r * q * b_b - y) * 
        (
            (1 - q) / (λ / (r + δ + λ))
        ) ^ (r / (r + δ))
    end
    v_safe = function (b) 
        1/ρ * (
            c_high - (c_high - ρ * v_high) * ( 
                (c_high + r * q_safe(b) * b - y) / 
                (c_high + r * q_low * b_b - y)
            ) ^ (ρ / r)
        )
    end 
    v_b = function (b) 
        if b ≤ b_b 
            v_safe(b)
        else
            v_crisis(b)
        end 
    end 
    q_b = function (b) 
        if b < 0.0 
            1.0
        elseif 0 ≤ b ≤ b_b 
            q_safe(b)
        else 
            q_low 
        end
    end
    return (q_b = q_b, v_b = v_b)
end


function get_hybrid_eqm(m, q_low, b_max_b, b_max_s, b_b, b_s, v_b, q_b)
    @unpack r, λ, δ, y, ρ, v_high, v_low, 
            c_high, c_low = m

    b_h = find_zero(b -> v_b(b) - (y - r * b) / ρ, (0.0, max(b_max_b, b_max_s)))
    v_h = function (b) 
        if b < b_h 
            (1/ρ) * (
                c_high - (c_high + r * b_h - y) * 
                (
                    (c_high + r * b - y) / 
                    (c_high + r * b_h - y)
                ) ^ (ρ / r)
            )
        else 
            v_b(b)
        end 
    end 
    q_h = function (b) 
        if b < b_h 
            1.0 
        else
            q_b(b)
        end 
    end 
    c_h = function (b) 
        if b < b_h 
            c_high 
        elseif b == b_h 
            y - r * b_h 
        elseif b < b_max_b 
            c_high 
        else 
            y - (r + δ * (1 - q_b(b_max_b))) * b_max_b
        end
    end 
    return (v_h = v_h, q_h = q_h, c_h = c_h, b_h = b_h)
end


function get_all(m)
    @unpack r, λ, δ, y, ρ, v_high, v_low, 
            c_high, c_low = m
    
    q_low = (r + δ) / (r + δ + λ)
    b_max_b = (
        (y + λ * v_high - (ρ + λ) * v_low) / 
        ((r + λ) * q_low)
    )
    b_b = (1 / ((r + λ) * q_low)) * (
        (c_high - ρ * v_high) ^ ((r + λ) / (ρ + λ)) *
        (c_high + λ * (v_high - v_low) - ρ * v_low) ^ ( (ρ - r) / (ρ + λ)) - 
        (c_high - y)
    )
    b_s = (y - ρ * v_high) / r

    savings = get_savings_eqm(m, q_low, b_max_b, b_b, b_s)
    @unpack c_s, q_s, v_s, b_i, b_max_s = savings 

    borrowing = get_borrowing_eqm(m, q_low, b_max_b, b_b, b_s)
    @unpack q_b, v_b = borrowing 

    hybrid = get_hybrid_eqm(m, q_low, b_max_b, b_max_s, b_b, b_s, v_b, q_b)
    @unpack v_h, q_h, c_h, b_h = hybrid

    c_ss = (b, q) -> y - (r + δ * (1 - q)) * b
    v_ss = function (b, q, v)
        if v ≥ v_high 
            return c_ss(b, q) / ρ
        else 
            (c_ss(b, q) + λ * v_high) / (ρ + λ)
        end
    end

    return ( 
        m = m, 
        q_low = q_low, 
        b_max_b = b_max_b,
        b_b = b_b,
        b_s = b_s, 
        c_ss = c_ss,
        v_ss = v_ss,
        c_ss_s = b -> c_ss(b, q_s(b)),
        c_ss_b = b -> c_ss(b, q_b(b)),
        c_ss_h = b -> c_ss(b, q_h(b)),
        v_ss_s = b -> v_ss(b, q_s(b), v_s(b)),
        v_ss_b = b -> v_ss(b, q_b(b), v_b(b)),
        v_ss_h = b -> v_ss(b, q_h(b), v_h(b)),
        c_b = b -> (b < b_max_b) ? c_high : c_ss(b_max_b, q_b(b_max_b)), 
        mkt_v_b = b -> b * q_b(b), 
        mkt_v_h = b -> b * q_s(b),
        p_star_s = get_pstar_s(m),
        p_star_b = get_pstar_b(m),
        savings..., borrowing..., hybrid...
    )
end 


function make_fig_borrowing_v(s)
    f = @pgf Axis(
        {
            width="3in",
            height="2.2in",
            xtick = [0.0, s.b_b, s.b_max_b] ,
            xticklabels=["0", "\$\\underline{b}_B\$", "\$\\overline{b}_B\$"],
            ytick = [s.m.v_low, s.m.v_high],
            yticklabels=["\$\\underline{V}\$", "\$\\overline{V}\$"]
        },
        VLine({style="lightgray, thick"}, 0.0),
        VLine({style="lightgray, thick"}, s.b_b),
        VLine({style="lightgray, thick"}, s.b_max_b),
        HLine({style="lightgray, thick"}, s.m.v_low),
        HLine({style="lightgray, thick"}, s.m.v_high),
        PlotInc(
            {
                style="blue, line width=3pt", mark="none"
            },
            Coordinates(
                [(b, s.v_b(b)) for b in range(0.0, s.b_max_b, length=20)]
            )
        ),
        PlotInc(
            {
                style="gray, ultra thick, dashed, dash pattern=on 5pt off 3pt", 
                mark="none"
            },
            Coordinates(
                [(b, s.v_ss_b(b)) for b in range(0.0, s.b_max_b, length=200)]
            )
        ),
        [
            "\\node [anchor=south west]",
            " at ",
            Coordinate(s.b_b, s.v_b(s.b_b)),
            "{\\textcolor{blue}{\$V_B(b)\$}};"
        ]
    )
end 


function make_fig_savings_v(s)
    f = @pgf Axis(
        {
            width="3in",
            height="2.2in",
            xtick = [0.0, s.b_s, s.b_i, s.b_max_s] ,
            xticklabels=[
                    "0", "\$\\underline{b}_S\$", "\$b^I\$", 
                    "\$\\overline{b}_S\$"],
            ytick = [s.m.v_low, s.m.v_high],
            yticklabels=["\$\\underline{V}\$", "\$\\overline{V}\$"]
        },
        VLine({style="lightgray, thick"}, 0.0),
        VLine({style="lightgray, thick"}, s.b_s),
        VLine({style="lightgray, thick"}, s.b_max_s),
        VLine({style="lightgray, thick"}, s.b_i),
        HLine({style="lightgray, thick"}, s.m.v_low),
        HLine({style="lightgray, thick"}, s.m.v_high),
        PlotInc(
            {
                style="blue, line width=3pt", mark="none"
            },
            Coordinates(
                [(b, s.v_s(b)) for b in range(0.0, s.b_max_s, length=20)]
            )
        ),
        PlotInc(
            {
                style="gray, ultra thick, dashed, dash pattern=on 5pt off 3pt", 
                mark="none"
            },
            Coordinates(
                [(b, s.v_ss_s(b)) for b in range(0.0, s.b_max_s, length=200)]
            )
        ),
        [
            "\\node [anchor=south west]",
            " at ",
            Coordinate(s.b_s / 2, s.v_s(s.b_s / 2)),
            "{\\textcolor{blue}{\$V_S(b)\$}};"
        ]
    )
end 


function make_fig_hybrid_v(s)
    f = @pgf Axis(
        {
            width="3in",
            height="2.2in",
            xtick = [0.0, s.b_h, s.b_b, s.b_max_b] ,
            xticklabels=[
                    "0", "\$b_H\$", "\$\\underline{b}_B\$", 
                    "\$\\overline{b}_B\$"],
            ytick = [s.m.v_low, s.m.v_high],
            yticklabels=["\$\\underline{V}\$", "\$\\overline{V}\$"]
        },
        VLine({style="lightgray, thick"}, 0.0),
        VLine({style="lightgray, thick"}, s.b_h),
        VLine({style="lightgray, thick"}, s.b_max_b),
        VLine({style="lightgray, thick"}, s.b_b),
        HLine({style="lightgray, thick"}, s.m.v_low),
        HLine({style="lightgray, thick"}, s.m.v_high),
        PlotInc(
            {
                style="blue, line width=3pt", mark="none"
            },
            Coordinates(
                [(b, s.v_h(b)) for b in range(0.0, s.b_max_b, length=20)]
            )
        ),
        PlotInc(
            {
                style="gray, ultra thick, dashed, dash pattern=on 5pt off 3pt", 
                mark="none"
            },
            Coordinates(
                [(b, s.v_ss_h(b)) for b in range(0.0, s.b_max_b, length=200)]
            )
        ),
        [
            "\\node [anchor=south west]",
            " at ",
            Coordinate(
                (s.b_b + s.b_max_b)/ 2, 
                s.v_h((s.b_b + s.b_max_b)/ 2)
            ),
            "{\\textcolor{blue}{\$V_H(b)\$}};"
        ]
    )
end 


function make_fig_borrowing_q(s)
    f = @pgf Axis(
        {
            width="3in",
            height="2.2in",
            xtick = [0.0, s.b_b, s.b_max_b] ,
            xticklabels=["0", "\$\\underline{b}_B\$", "\$\\overline{b}_B\$"],
            ytick = [1.0, s.q_low],
            yticklabels=["1", "\$\\underline{q}\$"]
        },
        VLine({style="lightgray, thick"}, 0.0),
        VLine({style="lightgray, thick"}, s.b_b),
        VLine({style="lightgray, thick"}, s.b_max_b),
        HLine({style="lightgray, thick"}, 1.0),
        HLine({style="lightgray, thick"}, s.q_low),
        PlotInc(
            {
                style="blue, line width=3pt", mark="none"
            },
            Coordinates(
                [(b, s.q_b(b)) for b in range(0.0, s.b_max_b, length=300)]
            )
        ),
        PlotInc(
            {
                style="blue, line width=3pt", mark="none"
            },
            Coordinates(
                [(b, s.q_b(b)) for b in range(-0.05, -0.000001, length=100)]
            )
        ),
        # white marker
        PlotInc(
            {
                style="blue", mark="*", 
                mark_options= "{line width=1pt, fill=white, scale=1.2}"
            },
            Coordinates(
                [(0.0, 1.0)]
            )
        ),
        # solid marker
        PlotInc(
            {
                style="blue", mark="*", 
                mark_options= "{line width=1pt, solid, scale=1.2}"
            },
            Coordinates(
                [(0.0, s.q_b(0.0))]
            )
        ),
        [
            "\\node [anchor=south]",
            " at ",
            Coordinate((s.b_b + s.b_max_b)/2, s.q_b(s.b_b)),
            "{\\textcolor{blue}{\$q_B(b)\$}};"
        ]
    )
end 


function make_fig_savings_q(s)
    f = @pgf Axis(
        {
            width="3in",
            height="2.2in",
            xtick = [0.0, s.b_s, s.b_i, s.b_max_s] ,
            xticklabels=["0", "\$\\underline{b}_S\$", 
                "\$b^I\$", "\$\\overline{b}_S\$"],
            ytick = [1.0, s.q_low],
            yticklabels=["1", "\$\\underline{q}\$"]
        },
        VLine({style="lightgray, thick"}, 0.0),
        VLine({style="lightgray, thick"}, s.b_s),
        VLine({style="lightgray, thick"}, s.b_i),
        VLine({style="lightgray, thick"}, s.b_max_s),
        HLine({style="lightgray, thick"}, 1.0),
        HLine({style="lightgray, thick"}, s.q_low),
        PlotInc(
            {
                style="blue, line width=3pt", mark="none"
            },
            Coordinates(
                [(b, s.q_s(b)) for b in range(-0.02, s.b_i-0.001, length=300)]
            )
        ),
        PlotInc(
            {
                style="blue, line width=3pt", mark="none"
            },
            Coordinates(
                [(b, s.q_s(b)) for b in range(s.b_i+0.001, s.b_max_s, length=10)]
            )
        ),
        PlotInc(
            {
                style="blue", mark="*", 
                mark_options= "{line width=1pt, fill=white, scale=1.2}"
            },
            Coordinates(
                [(s.b_i, s.q_low)]
            )
        ),
        PlotInc(
            {
                style="blue", mark="*", 
                mark_options= "{line width=1pt, solid, scale=1.2}"
            },
            Coordinates(
                [(s.b_i, s.q_s(s.b_i - 0.001))]
            )
        ),
        [
            "\\node [anchor=south]",
            " at ",
            Coordinate((s.b_i + s.b_max_s)/2, s.q_low),
            "{\\textcolor{blue}{\$q_S(b)\$}};"
        ]
    )
end 


function make_fig_hybrid_q(s)
    f = @pgf Axis(
        {
            width="3in",
            height="2.2in",
            xtick = [0.0, s.b_h, s.b_b, s.b_max_b] ,
            xticklabels=["0", "\$b_H\$", "\$\\underline{b}_B\$", 
                "\$\\overline{b}_B\$"],
            ytick = [1.0, s.q_low],
            yticklabels=["1", "\$\\underline{q}\$"]
        },
        VLine({style="lightgray, thick"}, 0.0),
        VLine({style="lightgray, thick"}, s.b_b),
        VLine({style="lightgray, thick"}, s.b_h),
        VLine({style="lightgray, thick"}, s.b_max_b),
        HLine({style="lightgray, thick"}, 1.0),
        HLine({style="lightgray, thick"}, s.q_low),
        PlotInc(
            {
                style="blue, line width=3pt", mark="none"
            },
            Coordinates(
                [(b, s.q_h(b)) for b in range(-0.02, s.b_h-0.001, length=300)]
            )
        ),
        PlotInc(
            {
                style="blue, line width=3pt", mark="none"
            },
            Coordinates(
                [(b, s.q_h(b)) for b in range(s.b_h+0.001, s.b_max_b, length=300)]
            )
        ),
        PlotInc(
            {
                style="blue", mark="*", 
                mark_options= "{line width=1pt, solid, scale=1.2}"
            },
            Coordinates(
                [(s.b_h, 1.0)]
            )
        ),
        PlotInc(
            {
                style="blue", mark="*", 
                mark_options= "{line width=1pt, fill=white, scale=1.2}"
            },
            Coordinates(
                [(s.b_h, s.q_h(s.b_h))]
            )
        ),
        [
            "\\node [anchor=south]",
            " at ",
            Coordinate((s.b_b + s.b_max_b)/2, s.q_low),
            "{\\textcolor{blue}{\$q_H(b)\$}};"
        ]
    )
end 


function make_fig_borrowing_c(s)
    f = @pgf Axis(
        {
            width="3in",
            height="2.2in",
            xtick = [0.0, s.b_b, s.b_max_b] ,
            xticklabels=["0", "\$\\underline{b}_B\$", "\$\\overline{b}_B\$"],
            ytick = [1],
            yticklabels=["\$y\$"]
        },
        VLine({style="lightgray, thick"}, 0.0),
        VLine({style="lightgray, thick"}, s.b_b),
        VLine({style="lightgray, thick"}, s.b_max_b),
        PlotInc(
            {
                style="blue, line width=3pt", mark="none"
            },
            Coordinates(
                [(b, s.c_b(b)) for b in range(0.0, s.b_max_b -0.0001, length=200)]
            )
        ),
        PlotInc(
            {
                style="gray, ultra thick, dashed, dash pattern=on 5pt off 3pt", 
                mark="none"
            },
            Coordinates(
                [(b, s.c_ss_b(b)) for b in range(0.0, s.b_max_b, length=200)]
            )
        ),
        # solid marker
        PlotInc(
            {
                style="blue", mark="*", 
                mark_options= "{line width=1pt, solid, scale=1.2}"
            },
            Coordinates(
                [(s.b_max_b, s.c_b(s.b_max_b))]
            )
        ),
        # white marker
        PlotInc(
            {
                style="blue", mark="*", 
                mark_options= "{line width=1pt, fill=white, scale=1.2}"
            },
            Coordinates(
                [(s.b_max_b, s.c_b(s.b_max_b - 0.0001))]
            )
        ),        
        [
            "\\node [anchor=north west]",
            " at ",
            Coordinate(s.b_b, s.c_b(s.b_b)),
            "{\\textcolor{blue}{\$\\mathcal{C}_B(b)\$}};"
        ]
    )
end 


function make_fig_savings_c(s)
    f = @pgf Axis(
        {
            width="3in",
            height="2.2in",
            xtick = [0.0, s.b_s, s.b_i, s.b_max_s],
            xticklabels=[
                    "0", "\$\\underline{b}_S\$", "\$b^I\$", 
                    "\$\\overline{b}_S\$"],
            ytick = [1],
            yticklabels=["\$y\$"]
        },
        VLine({style="lightgray, thick"}, 0.0),
        VLine({style="lightgray, thick"}, s.b_s),
        VLine({style="lightgray, thick"}, s.b_max_s),
        VLine({style="lightgray, thick"}, s.b_i),
        # consumption
        PlotInc(
            {
                style="blue, line width=3pt", mark="none"
            },
            Coordinates(
                [(b, s.c_s(b)) for b in range(0.0, s.b_s -0.0001, length=200)]
            )
        ),
        PlotInc(
            {
                style="blue, line width=3pt", mark="none"
            },
            Coordinates(
                [(b, s.c_s(b)) for b in range(s.b_s + 0.0001, s.b_i -0.0001, length=200)]
            )
        ),
        PlotInc(
            {
                style="blue, line width=3pt", mark="none"
            },
            Coordinates(
                [(b, s.c_s(b)) for b in range(s.b_i + 0.00001, s.b_max_s -0.0001, length=200)]
            )
        ),
        # stationary consumption
        PlotInc(
            {
                style="gray, ultra thick, dashed, dash pattern=on 5pt off 3pt", 
                mark="none"
            },
            Coordinates(
                [(b, s.c_ss_s(b)) for b in range(0.0, s.b_max_s, length=500)]
            )
        ),
        # solid markers
        PlotInc(
            {
                style="blue", mark="*", only_marks,
                mark_options= "{line width=1pt, solid, scale=1.2}"
            },
            Coordinates(
                [
                 (s.b_max_s, s.c_s(s.b_max_s)),
                 (s.b_s, s.c_s(s.b_s)),
                 (s.b_i, s.c_s(s.b_i))
                ]
            )
        ),
        # white markers
        PlotInc(
            {
                style="blue", mark="*", only_marks,
                mark_options= "{line width=1pt, solid, fill=white, scale=1.2}"
            },
            Coordinates(
                [
                    (s.b_max_s, s.c_s(s.b_max_s - 0.0001)),
                    (s.b_s, s.c_s(s.b_s - 0.0001)), 
                    (s.b_s, s.c_s(s.b_s + 0.00001)),
                    (s.b_i, s.c_s(s.b_i + 0.0001))
                ]
            )
        ),
        [
            "\\node [anchor=north west]",
            " at ",
            Coordinate(s.b_s/3, s.c_s(s.b_s /3)),
            "{\\textcolor{blue}{\$\\mathcal{C}_S(b)\$}};"
        ]
    )
end 


function make_fig_hybrid_c(s)
    f = @pgf Axis(
        {
            width="3in",
            height="2.2in",
            xtick = [0.0, s.b_h, s.b_b, s.b_max_b] ,
            xticklabels=["0", "\$b_H\$", "\$\\underline{b}_B\$", 
                "\$\\overline{b}_B\$"],
            ytick = [1.0],
            yticklabels=["\$y\$"]
        },
        VLine({style="lightgray, thick"}, 0.0),
        VLine({style="lightgray, thick"}, s.b_b),
        VLine({style="lightgray, thick"}, s.b_h),
        VLine({style="lightgray, thick"}, s.b_max_b),
        # consumption
        PlotInc(
            {
                style="blue, line width=3pt", mark="none"
            },
            Coordinates(
                [(b, s.c_h(b)) for b in range(0.0, s.b_h -0.0001, length=200)]
            )
        ),
        PlotInc(
            {
                style="blue, line width=3pt", mark="none"
            },
            Coordinates(
                [(b, s.c_h(b)) for b in range(s.b_h + 0.0001, s.b_max_b - 0.0001, length=200)]
            )
        ),
        # stationary consumption
        PlotInc(
            {
                style="gray, ultra thick, dashed, dash pattern=on 5pt off 3pt", 
                mark="none"
            },
            Coordinates(
                [(b, s.c_ss_h(b)) for b in range(0.0, s.b_max_b, length=500)]
            )
        ),
        # solid markers
        PlotInc(
            {
                style="blue", mark="*", only_marks,
                mark_options= "{line width=1pt, solid, scale=1.2}"
            },
            Coordinates(
                [
                 (s.b_max_b, s.c_h(s.b_max_b)),
                 (s.b_h, s.c_h(s.b_h))
                ]
            )
        ),
        # white markers
        PlotInc(
            {
                style="blue", mark="*", only_marks,
                mark_options= "{line width=1pt, solid, fill=white, scale=1.2}"
            },
            Coordinates(
                [
                    (s.b_max_b, s.c_h(s.b_max_b - 0.0001)),
                    (s.b_h, s.c_h(s.b_h - 0.0001))
                ]
            )
        ),
        [
            "\\node [anchor=north west]",
            " at ",
            Coordinate(s.b_max_b/2, s.c_h(s.b_max_b /2)),
            "{\\textcolor{blue}{\$\\mathcal{C}_H(b)\$}};"
        ]
    )
end 


function make_joint_surplus_savings(s)

    f = @pgf Axis(
        {
            width="3in",
            height="2.2in",
            xtick = [s.m.v_low, s.m.v_high, s.v_s(s.b_i)],
            xlabel = "Government Value, \$ v\$",
            ylabel = "Lender Value, \$ P\$",
            xticklabels=["\$\\underline{V}\$", 
                "\$\\overline{V}\$", 
                "\$V_S(b^I)\$"],
            ytick = [],
            yticklabels=[]
        },
        VLine({style="lightgray, thick"}, s.m.v_low),
        VLine({style="lightgray, thick"}, s.m.v_high),
        VLine({style="lightgray, thick"}, s.v_s(s.b_i)),
        # equilibriun
        PlotInc(
            {
                style="blue, line width=3pt", mark="none"
            },
            Coordinates(
                [(s.v_s(b), s.q_s(b) * b) for b in range(0.0, s.b_max_s, length=300)]
            )
        ),
        # PS
        PlotInc(
            {
                style="gray, ultra thick, dashed, dash pattern=on 5pt off 3pt", 
                mark="none"
            },
            Coordinates(
                [(s.v_s(b), s.p_star_s(s.v_s(b))) for b in range(0.0, s.b_max_s, length=200)]
            )
        ),
        # PB
        PlotInc(
            {
                style="gray, ultra thick, dashed, dash pattern=on 5pt off 3pt", 
                mark="none"
            },
            Coordinates(
                [(s.v_s(b), s.p_star_b(s.v_s(b))) for b in range(0.0, s.b_max_s, length=200)]
            )
        ),
        [
            "\\node [anchor=south west]",
            " at ",
            Coordinate(
                s.v_s(s.b_max_s)  + 0.015, 
                s.p_star_s(s.v_s(s.b_max_s)) - 0.02
            ),
            "{\\textcolor{gray}{\$P_S(v)\$}};"
        ],
        [
            "\\node [anchor=north east]",
            " at ",
            Coordinate(
                s.v_b(0.0), 
                s.p_star_b(s.v_b(0.0)) + .005
            ),
            "{\\textcolor{gray}{\$P_B(v)\\, \$}};"
        ],
        [
            "\\node [anchor=south]",
            " at ",
            Coordinate(
                s.v_b(0.0), 
                s.p_star_s(s.v_b(0.0)) + .012
            ),
            "{\\textcolor{blue}{\$q_S(b) b\\, \$}};"
        ]
    )
end


function make_joint_surplus_borrowing(s)

    f = @pgf Axis(
        {
            width="3in",
            height="2.2in",
            xtick = [s.m.v_low, s.m.v_high],
            xlabel = "Government Value, \$ v\$",
            ylabel = "Lender Value, \$ P\$",
            xticklabels=["\$\\underline{V}\$", 
                "\$\\overline{V}\$"],
            ytick = [],
            yticklabels=[]
        },
        VLine({style="lightgray, thick"}, s.m.v_low),
        VLine({style="lightgray, thick"}, s.m.v_high),
        # equilibriun
        PlotInc(
            {
                style="blue, line width=3pt", mark="none"
            },
            Coordinates(
                [(s.v_b(b), s.q_b(b) * b) for b in range(0.0, s.b_max_b, length=300)]
            )
        ),
        # PS
        PlotInc(
            {
                style="gray, ultra thick, dashed, dash pattern=on 5pt off 3pt", 
                mark="none"
            },
            Coordinates(
                [(s.v_b(b), s.p_star_s(s.v_b(b))) for b in range(0.0, s.b_max_b, length=200)]
            )
        ),
        # PB
        PlotInc(
            {
                style="gray, ultra thick, dashed, dash pattern=on 5pt off 3pt", 
                mark="none"
            },
            Coordinates(
                [(s.v_b(b), s.p_star_b(s.v_b(b))) for b in range(0.0, s.b_max_b, length=200)]
            )
        ),
        [
            "\\node [anchor=south west]",
            " at ",
            Coordinate(
                s.v_s(s.b_max_s/3), 
                s.p_star_s(s.v_s(s.b_max_s/3))
            ),
            "{\\textcolor{gray}{\$P_S(v)\$}};"
        ],
        [
            "\\node [anchor=north east]",
            " at ",
            Coordinate(
                s.v_b(s.b_max_b /2), 
                s.p_star_b(s.v_b(s.b_max_b /2))
            ),
            "{\\textcolor{blue}{\$P_B(v)\\, \$}};"
        ]
    )
end


function save_figures(; m=Model(), dir="Figures")

    f(x) = joinpath("Figures", x)

    sol = get_all(m)
    pgfsave(f("plotVB.pdf"), make_fig_borrowing_v(sol))
    pgfsave(f("plotqB.pdf"), make_fig_borrowing_q(sol))
    pgfsave(f("plotCB.pdf"), make_fig_borrowing_c(sol))
 
    pgfsave(f("plotVS.pdf"), make_fig_savings_v(sol))
    pgfsave(f("plotqS.pdf"), make_fig_savings_q(sol))
    pgfsave(f("plotCS.pdf"), make_fig_savings_c(sol))

    pgfsave(f("plotVH.pdf"), make_fig_hybrid_v(sol))
    pgfsave(f("plotqH.pdf"), make_fig_hybrid_q(sol))
    pgfsave(f("plotCH.pdf"), make_fig_hybrid_c(sol))

    pgfsave(f("plotPB.pdf"), make_joint_surplus_borrowing(sol))
    pgfsave(f("plotPS.pdf"), make_joint_surplus_savings(sol))
end 


push!(
    PGFPlotsX.CUSTOM_PREAMBLE, 
    "\\usepackage{libertine}\\usepackage[libertine]{newtxmath}"
)
save_figures()