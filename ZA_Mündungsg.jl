using  Parameters, Plots, LinearAlgebra, GeometryBasics

# ------Code zum berechnen der Mündungsgeschwindigkeit des Projektils
# Also snippets genommen aus FieldSim
# ─────────────────────────── physical constants & geometry ───────────────────────────
μ0       = 4π * 1e-7        # vacuum permeability (N A⁻²)
I_total  = 400.0            # current per rail (A)

l_prj = 0.01                # width of projectile, aka just the gap lol
m_prj = 0.010                # mass of prj 10 grams
# Rail cross‑section (looking into the barrel)
rail_h   = 0.025             # height  (y‑direction)  = 2.5 cm
rail_len = 0.04             # length  (x‑direction)  = 4 cm
gap      = 0.01             # air gap between inner faces = 1 cm
rail_depth = 0.8             # length of rails = 80cm
rail_x1  = 0.0                               # inner face of left rail
rail_x2  = rail_x1 + rail_len + gap          # inner face of right rail
rail_y0, rail_y1 = 0.0, rail_h               # bottom / top of rails

# ───────────────────────── discretise each rail into Ny×Nx filaments / Biot Savart Prep ────────────────
Nx, Ny = 30, 20                               # along x   ×   along y (resolution)
xs1 = range(rail_x1 + rail_len/(2Nx), rail_x1 + rail_len - rail_len/(2Nx); length = Nx)
xs2 = xs1 .+ (rail_x2 - rail_x1)              # shift for right rail
ys  = range(rail_y0 + rail_h /(2Ny), rail_y1 - rail_h /(2Ny); length = Ny)
dI  = I_total / (Nx*Ny)                       # current per filament (A)

filaments = vcat(
    [(x, y, +1.0) for x in xs1, y in ys],    # left rail  (+z)
    [(x, y, -1.0) for x in xs2, y in ys],    # right rail (−z)
)
p_gap = Point2(rail_x1 + rail_len + gap/2,   # x-Koordinate
               rail_h/2)                     # y-Koordinate
# ─────────────────────────────── B‑field function ───────────────────────────────────
function B(p::Point2)
    bx = 0.0;  by = 0.0
    for (x0, y0, s) in filaments
        dx = p[1] - x0
        dy = p[2] - y0
        r2 = dx^2 + dy^2 + 1e-20
        # ----- lineares (sehr langes) Leitersegment -----
        fac = s * μ0 * dI / (2π * r2)   # 1/r
        bx -= fac * dy                  #   φ-Komponente
        by += fac * dx
    end
    return Point2f(bx, by)
end
# ─────────────────────────────── Lorentz Force function ───────────────────────────────────
function LorentzF(rail_depth)
    By = 2* B(p_gap)[2]  
    return I_total*l_prj*By
end
FL = LorentzF(rail_depth)   # this is the part where i make a typical assumption of no friction because i engineered the rails and projectile perfect
# ─────────────────────────────── Basic Mechanics ───────────────────────────────────

function acc(LorentzF, m_prj)
    return LorentzF / m_prj
end
a = acc(FL, m_prj)


# ─────────────────────────────── OUTPUT ───────────────────────────────────
println("Magnetfeld Stärke B=", 2* B(p_gap)[2])
println("Lorentz Kraft F_L=", LorentzF(rail_depth))
println("Beschleunigung a=", acc(FL, m_prj))


