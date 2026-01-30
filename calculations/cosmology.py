#!python3

import math
import argparse
from dataclasses import dataclass

# ==========================================
# 0. EXTERNAL REFERENCE VALUES (Source of Truth)
# ==========================================
@dataclass
class MeasuredVal:
    value: float
    uncertainty: float
    units: str
    citation: str

# Citations:
# [1] Planck Collaboration, "Planck 2018 results. VI. Cosmological parameters," Astron. Astrophys. 641, A6 (2020).
# [2] A. G. Riess et al. (SH0ES), "A Comprehensive Measurement of the Local Value of the Hubble Constant...", Astrophys. J. Lett. 934, L7 (2022).
# [3] T. M. C. Abbott et al. (DES), "Dark Energy Survey Year 3 Results...", Phys. Rev. D 105, 023520 (2022).
# [4] PDG 2024, R. L. Workman et al., "Review of Particle Physics," Prog. Theor. Exp. Phys. 2022, 083C01 (2022) and 2023 update.

REFS = {
    # Primordial Parameters (Planck 2018)
    "ns": MeasuredVal(0.9649, 0.0042, "", "Planck2018"),
    "eta": MeasuredVal(6.12e-10, 0.04e-10, "", "Planck2018"),
    "zeq": MeasuredVal(3387, 21, "", "Planck2018"),
    
    # Tension Parameters
    "H0_early": MeasuredVal(67.4, 0.5, "km/s/Mpc", "Planck2018"),
    "H0_late": MeasuredVal(73.0, 1.0, "km/s/Mpc", "Riess2022"),
    "S8_planck": MeasuredVal(0.832, 0.013, "", "Planck2018"),
    "S8_lensing": MeasuredVal(0.776, 0.017, "", "DES-Y3"), # DES Y3 value
    
    # Dark Sector Parameters (Planck 2018)
    "omega_ratio": MeasuredVal(5.36, 0.05, "", "Planck2018"), # Omega_c / Omega_b = 0.120 / 0.0224
    "neutrino_sum_cosmo": MeasuredVal(0.12, 0.0, "eV", "Planck2018"), # This is an upper limit, handled specially
    "neutrino_sum_osc": MeasuredVal(0.06, 0.0, "eV", "PDG2024"), # This is a lower limit, handled specially
}


# ==========================================
# 1. HELPER FUNCTIONS
# ==========================================

def format_float_latex(num, precision=9):
    return f"{num:.{precision}f}".rstrip('0').rstrip('.')

def to_latex_sci(num, precision=4):
    if num == 0: return "0"
    exponent = int(math.floor(math.log10(abs(num))))
    mantissa = num / (10**exponent)
    if -3 <= exponent < 6:
        if abs(num - round(num)) < 1e-9:
            return f"{int(num)}"
        return f"{num:.{precision}f}".rstrip('0').rstrip('.')
    return f"{mantissa:.{precision}f} \\times 10^{{{exponent}}}"

def to_latex_sci_with_err(val, err, precision=4):
    if val == 0: return f"0 \\pm {err}"
    exponent = int(math.floor(math.log10(abs(val))))
    mantissa_val = val / (10**exponent)
    mantissa_err = err / (10**exponent)
    return f"({mantissa_val:.{precision}f} \\pm {mantissa_err:.{precision}f}) \\times 10^{{{exponent}}}"

def print_latex_tag(tag, content):
    print(f"%<*{tag}>{content}%</{tag}>")

def mag(z):
    """Safe magnitude for complex or float"""
    return abs(z)

def print_derivation(name, tag, latex_eq, calc_val, ref_key, unit="", latex_mode=False):
    # If unit is empty but ref has one, use ref's unit
    if ref_key and ref_key in REFS and unit == "":
        unit = REFS[ref_key].units

    if latex_mode:
        # Value
        val_str = format_float_latex(calc_val) if 0.001 < abs(calc_val) < 1000 else to_latex_sci(calc_val, 5)
        if unit: val_str += f" \\unit{{{unit}}}"
        print_latex_tag(tag + "Val", val_str)

        # Equation
        if latex_eq:
            print_latex_tag(tag + "Eq", latex_eq)

        # Experimental comparison
        if ref_key and ref_key in REFS:
            ref = REFS[ref_key]
            cite_str = f"~\\cite{{{ref.citation}}}" if ref.citation else ""
            if 0.001 < abs(ref.value) < 1000:
                exp_str = f"\\qty{{{format_float_latex(ref.value)} \\pm {format_float_latex(ref.uncertainty)}}}{{{unit}}}{cite_str}"
            else:
                exp_str = f"${to_latex_sci_with_err(ref.value, ref.uncertainty, 5)}$ {cite_str}"
            
            print_latex_tag(tag + "ExperimentalValue", exp_str)

            diff = calc_val - ref.value
            sigma = abs(diff / ref.uncertainty) if ref.uncertainty > 0 else 0.0
            
            print_latex_tag(tag + "Diff", f"{diff:.2e}")
            print_latex_tag(tag + "Sigma", f"{sigma:.2f}")

            if sigma < 1.0:
                acc_text = f"matches the experimental consensus to within ${sigma:.2f}\\sigma$."
            elif sigma < 3.0:
                acc_text = f"lies within ${sigma:.2f}\\sigma$ of the observed value."
            else:
                acc_text = f"deviates by ${sigma:.2f}\\sigma$ from the observed value."
            print_latex_tag(tag + "AccText", acc_text)
        return

    # Console Mode
    print(f"--- {name} ---")
    if latex_eq: print(f"Formula:    {latex_eq}")
    print(f"Calculated: {calc_val:.12f} {unit}")
    if ref_key and ref_key in REFS:
        ref = REFS[ref_key]
        sigma = (calc_val - ref.value) / ref.uncertainty if ref.uncertainty > 0 else 0
        print(f"Target:     {ref.value:.6f} +/- {ref.uncertainty:.6f}")
        print(f"Deviation:  {sigma:+.2f}σ")
    print("")

# ==========================================
# 2. MAIN EXECUTION
# ==========================================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--latex', action='store_true')
    args = parser.parse_args()

    # --- INVARIANTS (from Paper I) ---
    D, DELTA, SIGMA, NU, CHI = 4, 43, 5, 16, 2
    ALPHA_INV = 137.035999212
    ALPHA = 1.0 / ALPHA_INV
    JARLSKOG = 3.08490e-5 # Derived in Paper I/III
    H_SYS = NU + SIGMA + CHI # Systemic Channel = 23
    ME_EV = 510998.95 # Electron mass in eV
    PI = math.pi

    if not args.latex:
        print("\n" + "="*60)
        print("  PAPER IV: COSMOLOGY CALCULATIONS")
        print("="*60 + "\n")

    # --- 1. PRIMORDIAL UNIVERSE ---
    if not args.latex: print("--- SECTION II: PRIMORDIAL UNIVERSE ---")
    
    # Spectral Index (ns)
    ns_geo = 1.0 - (1.0 / (2 * NU))
    print_derivation("Spectral Index (ns)", "SpectralIndex", r"1 - \frac{1}{2\nu}", ns_geo, "ns", "", args.latex)

    # Baryon-to-Photon Ratio (eta)
    epsilon = JARLSKOG / NU
    zeta = 1.0 / ((PI * DELTA) * (H_SYS + CHI / D))
    eta_geo = epsilon * zeta
    print_derivation("Baryon-to-Photon Ratio (eta)", "EtaBaryon", r"\frac{J}{\nu} \cdot \frac{1}{(\pi\Delta)(H_{sys} + \chi/D)}", eta_geo, "eta", "", args.latex)

    # Matter-Radiation Equality (zeq)
    zeq_geo = (PI * DELTA) * (H_SYS + CHI)
    print_derivation("Matter-Radiation Equality (zeq)", "zeq", r"(\pi\Delta)(H_{sys} + \chi)", zeq_geo, "zeq", "", args.latex)

    # --- 2. LATE UNIVERSE TENSIONS ---
    if not args.latex: print("\n--- SECTIONS III-V: LATE UNIVERSE TENSIONS ---")
    
    # Lattice Load Factor (L)
    sq = 1.0 / (2.0 * NU)
    st = CHI / DELTA
    load_factor = 1.0 + sq + st
    print_derivation("Lattice Load Factor (L)", "LoadFactor", r"1 + \frac{1}{2\nu} + \frac{\chi}{\Delta}", load_factor, None, "", args.latex)

    # Hubble Tension
    h0_early_ref = REFS["H0_early"].value
    h0_late_geo = h0_early_ref * load_factor
    print_derivation("Hubble Constant (Late)", "HubbleLate", r"H_0^{early} \times \mathcal{L}", h0_late_geo, "H0_late", "km/s/Mpc", args.latex)

    # S8 Tension
    s8_planck_ref = REFS["S8_planck"].value
    s8_geo = s8_planck_ref / load_factor
    print_derivation("S8 Parameter (Late)", "S8Geo", r"S_8^{Planck} / \mathcal{L}", s8_geo, "S8_lensing", "", args.latex)

    # --- 3. THE DARK SECTOR ---
    if not args.latex: print("\n--- SECTION VII: THE DARK SECTOR ---")
    
    # Dark Matter Abundance Ratio
    dim_e8 = 248.0
    active_channels = 3.0 * NU
    binding_corr = CHI / DELTA
    charge_corr = ALPHA
    denom_corr = active_channels * (1.0 - binding_corr + charge_corr)
    omega_ratio_geo = dim_e8 / denom_corr
    print_derivation("DM/Baryon Ratio", "OmegaRatio", r"\frac{\text{dim}(E_8)}{3\nu(1 - \chi/\Delta + \alpha)}", omega_ratio_geo, "omega_ratio", "", args.latex)

    # Neutrino Mass Paradox
    m_min_ev = (ME_EV * PI * ALPHA) / (NU * (DELTA**3))
    sum_mnu_rest = m_min_ev * (1.0 + (1.0 / DELTA) + (1.0 / (DELTA**2)))
    eta_w = SIGMA + 1.0
    sum_mnu_weak = sum_mnu_rest * eta_w

    if args.latex:
        print_latex_tag("MnuRestVal", f"{sum_mnu_rest * 1000:.1f}") # in meV
        print_latex_tag("MnuWeakVal", f"{sum_mnu_weak * 1000:.0f}") # in meV
    else:
        print("--- Neutrino Mass Paradox ---")
        print(f"Gravitational Rest Mass (Cosmology): {sum_mnu_rest:.6f} eV ({sum_mnu_rest*1000:.1f} meV)")
        print(f"Target (Cosmology): < {REFS['neutrino_sum_cosmo'].value:.2f} eV (Planck Limit)")
        print(f"Interaction Mass (Oscillations):   {sum_mnu_weak:.6f} eV ({sum_mnu_weak*1000:.0f} meV)")
        print(f"Target (Oscillations): > {REFS['neutrino_sum_osc'].value:.2f} eV (from splittings)")
        print("")


if __name__ == "__main__":
    main()
