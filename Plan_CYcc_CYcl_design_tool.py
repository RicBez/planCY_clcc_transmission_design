#!/usr/bin/env python3
"""
Combined Cycloidal–Planetary Reducer Interactive Design Tool
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


# ======================================
# ---- BASIC GEOMETRY UTILITIES ----
# ======================================

def make_circle(cx, cy, radius, npts=200):
    """Return circle coordinates centered at (cx, cy) with radius."""
    ang = np.linspace(0, 2*np.pi, npts)
    return cx + radius*np.cos(ang), cy + radius*np.sin(ang)

def make_cycloidal_profile(eps, roller_r, nrollers, pitch_R, resolution=2000):
    """Generate a cycloidal cam disk profile."""
    pts = []
    theta = np.linspace(0, 2*np.pi, resolution)
    for t in theta:
        delta = np.arctan(np.sin((1 - nrollers)*t) /
                          ((pitch_R / (eps*nrollers)) - np.cos((1 - nrollers)*t)))
        x = eps + pitch_R*np.cos(t) - roller_r*np.cos(t + delta) - eps*np.cos(nrollers*t)
        y = -pitch_R*np.sin(t) + roller_r*np.sin(t + delta) + eps*np.sin(nrollers*t)
        pts.append((x, y))
    pts.append(pts[0])
    return np.array(pts)

def compute_Dmin(N, Dr):
    return N*Dr/math.pi

# =========================
# ---- PARAM SUMMARY -----
# =========================

def print_summary(params, stage=None):
    """
    Print a compact summary of parameters chosen so far.
    If stage provided, print only params relevant to that stage; otherwise print all.
    """
    print("\n--- RIEPILOGO PARAMETRI ---")
    if not params:
        print("Nessun parametro impostato ancora.")
        print("----------------------------\n")
        return
    if stage is None or stage == "planetary":
        p = params.get("planetary", {})
        if p:
            print("Planetario:")
            for k, v in p.items():
                print(f"  {k}: {v}")
    if stage is None or stage.startswith("cycloidal"):
        c = params.get("cycloidal", {})
        if c:
            print("Cycloidale (opzione: {}):".format(c.get("variant", "compact-cam")))
            for k, v in c.items():
                if k == "variant":
                    continue
                print(f"  {k}: {v}")
    print("----------------------------\n")


# =========================
# ---- PLANETARY STAGE ----
# =========================
def design_planetary_stage(existing_params=None):
    """
    Interactive design of planetary stage.
    existing_params: dict to prefill or modify
    Returns a dict with planetary parameters.
    """
    import math
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    params = existing_params.copy() if existing_params else {}
    params_planet = params.get("planetary", {})

    steps = [
        "R_available",
        "ratio_req",
        "mod",
        "ring_teeth_max_and_choice",
        "select_solution",
        "select_N_planets"
    ]

    solutions = []
    step_index = 0
    while step_index < len(steps):
        step = steps[step_index]

        if step == "R_available":
            R_available = float(input("Enter the maximum allowed encumbrance radius for the planetary stage [mm]: "))
            params_planet["R_available_mm"] = R_available
            step_index += 1

        elif step == "ratio_req":
            ratio_req = float(input("Enter desired reduction ratio (sun→carrier, ring fixed) [3:9]: "))
            if not (3 <= ratio_req <= 9):
                print("ERROR: Ratio must be between 3 and 9. Please re-enter.")
                continue
            params_planet["ratio_req"] = ratio_req
            step_index += 1

        elif step == "mod":
            mod = float(input("Enter gear module m [mm]: "))
            params_planet["m_mm"] = mod
            margin = 2.25 * mod
            params_planet["margin_mm"] = margin
            ring_teeth_max = math.floor(2 * (params_planet["R_available_mm"] - margin) / mod)
            params_planet["ring_teeth_max"] = ring_teeth_max
            step_index += 1

        elif step == "ring_teeth_max_and_choice":
            print(f"\n→ Maximum feasible ring teeth (internal gear): {params_planet['ring_teeth_max']}")
            use_auto = input("Use this value? (y/n): ").strip().lower()
            if use_auto != "y":
                ring_teeth = int(input("Enter custom ring teeth count: "))
            else:
                ring_teeth = params_planet["ring_teeth_max"]
            params_planet["ring_teeth"] = ring_teeth
            step_index += 1

        elif step == "select_solution":
            solutions = []
            z_ring = params_planet["ring_teeth"]
            mod = params_planet["m_mm"]
            ratio_req = params_planet["ratio_req"]

            for z_sun in range(8, z_ring - 8):
                z_planet = (z_ring - z_sun) / 2
                if z_planet.is_integer():
                    ratio = z_ring / z_sun
                    if abs(ratio - ratio_req) < 0.2:
                        valid_N = []
                        for npl in range(3, 8):
                            if (z_ring - z_sun) % npl == 0:
                                valid_N.append(npl)
                        if valid_N:
                            r_sun = z_sun * mod / 2
                            r_ring = z_ring * mod / 2
                            r_carrier = (r_sun + r_ring) / 2
                            solutions.append({
                                "Zs": z_sun,
                                "Zr": z_ring,
                                "Zp": int(z_planet),
                                "ratio": ratio,
                                "r_sun_mm": r_sun,
                                "r_carrier_mm": r_carrier,
                                "valid_planet_counts": valid_N
                            })

            if not solutions:
                print("\n No valid combinations found with the current parameters.")
                choice = input("Do you want to modify parameters? (y/n): ").strip().lower()
                if choice == "y":
                    print("\nWhat would you like to modify?")
                    print(" 1) Module (m)")
                    print(" 2) Ring teeth number (Zr)")
                    print(" 3) Both")
                    sub_choice = input("Select option (1/2/3): ").strip()

                    if sub_choice == "1":
                        # Go back to step 'mod'
                        step_index = steps.index("mod")
                    elif sub_choice == "2":
                        # Go back to step 'ring_teeth_max_and_choice'
                        step_index = steps.index("ring_teeth_max_and_choice")
                    elif sub_choice == "3":
                        # Go back to module step (it will recompute everything)
                        step_index = steps.index("mod")
                    else:
                        print("Invalid choice. Returning to module input.")
                        step_index = steps.index("mod")
                    continue  # restart from selected step
                else:
                    print("Exiting planetary design — no valid configuration found.")
                    return None

            # Display possible solutions
            print("\nFeasible configurations (index: Zs, Zr, Zp, ratio, r_sun, r_carrier, possible planet counts):")
            for i, s in enumerate(solutions, 1):
                print(f"{i}: {s}")
            sel = int(input("\nSelect desired configuration index: ")) - 1
            chosen = solutions[sel]
            params_planet.update(chosen)
            step_index += 1


        elif step == "select_N_planets":
            N_list = params_planet["valid_planet_counts"]
            print(f"Possible numbers of planets: {N_list}")
            idxN = int(input(f"Select number of planets index from 1 to {len(N_list)}): ")) - 1
            n_planets = N_list[idxN]
            params_planet["N_planets"] = n_planets

            # compute diameters (foot, pitch, head radii)
            def diameters(z, mod_local):
                d_pitch = z * mod_local
                d_add = d_pitch + 2 * mod_local
                d_ded = d_pitch - 2.5 * mod_local
                return d_ded / 2, d_pitch / 2, d_add / 2

            r_sun_ded, r_sun_pitch, r_sun_add = diameters(params_planet["Zs"], params_planet["m_mm"])
            r_ring_ded, r_ring_pitch, r_ring_add = diameters(params_planet["Zr"], params_planet["m_mm"])
            r_planet_ded, r_planet_pitch, r_planet_add = diameters(params_planet["Zp"], params_planet["m_mm"])

            params_planet.update({
                "r_sun_ded_mm": r_sun_ded,
                "r_sun_pitch_mm": r_sun_pitch,
                "r_sun_add_mm": r_sun_add,
                "r_ring_ded_mm": r_ring_ded,
                "r_ring_pitch_mm": r_ring_pitch,
                "r_ring_add_mm": r_ring_add,
                "r_planet_ded_mm": r_planet_ded,
                "r_planet_pitch_mm": r_planet_pitch,
                "r_planet_add_mm": r_planet_add
            })

            # Plotting
            while True:
                fig, ax = plt.subplots(figsize=(8, 8))
                ax.set_aspect("equal")
                ax.set_title("Planetary Gear Stage – Geometric Layout", fontsize=14)

                # Draw ring
                ring_circles = [
                    (params_planet["r_ring_ded_mm"], 'black', ':'),
                    (params_planet["r_ring_pitch_mm"], 'grey', '-'),
                    (params_planet["r_ring_add_mm"], 'black', '--')
                ]
                for radius, color, style in ring_circles:
                    ax.add_artist(plt.Circle((0, 0), radius, fill=False, color=color, linestyle=style, linewidth=1.5))

                # Draw sun
                sun_circles = [
                    (params_planet["r_sun_ded_mm"], 'red', '--'),
                    (params_planet["r_sun_pitch_mm"], 'orange', '-'),
                    (params_planet["r_sun_add_mm"], 'red', ':')
                ]
                for radius, color, style in sun_circles:
                    ax.add_artist(plt.Circle((0, 0), radius, fill=False, color=color, linestyle=style, linewidth=1.5))

                # Draw planets
                r_carrier = params_planet["r_carrier_mm"]
                for i in range(params_planet["N_planets"]):
                    angle = 2 * math.pi * i / params_planet["N_planets"]
                    x_c = r_carrier * math.cos(angle)
                    y_c = r_carrier * math.sin(angle)
                    planet_circles = [
                        (params_planet["r_planet_ded_mm"], 'green', '--'),
                        (params_planet["r_planet_pitch_mm"], 'lime', '-'),
                        (params_planet["r_planet_add_mm"], 'green', ':')
                    ]
                    for radius, color, style in planet_circles:
                        ax.add_artist(plt.Circle((x_c, y_c), radius, fill=False, color=color, linestyle=style, linewidth=1.2))

                # Legend
                legend_elements = [
                    Line2D([0], [0], color='black', linestyle='--', lw=1.5, label='Ring – Dedendum circle (tooth root)'),
                    Line2D([0], [0], color='grey', linestyle='-', lw=1.5, label='Ring – Pitch circle (primitive)'),
                    Line2D([0], [0], color='black', linestyle=':', lw=1.5, label='Ring – Addendum circle (tooth tip)'),
                    Line2D([0], [0], color='red', linestyle='--', lw=1.5, label='Sun – Dedendum circle'),
                    Line2D([0], [0], color='orange', linestyle='-', lw=1.5, label='Sun – Pitch circle'),
                    Line2D([0], [0], color='red', linestyle=':', lw=1.5, label='Sun – Addendum circle'),
                    Line2D([0], [0], color='green', linestyle='--', lw=1.2, label='Planet – Dedendum circle'),
                    Line2D([0], [0], color='lime', linestyle='-', lw=1.2, label='Planet – Pitch circle'),
                    Line2D([0], [0], color='green', linestyle=':', lw=1.2, label='Planet – Addendum circle')
                ]

                ax.legend(handles=legend_elements, loc='upper right', fontsize=8, frameon=True)
                ax.set_xlabel("mm")
                ax.set_ylabel("mm")
                ax.grid(True, linestyle='--', alpha=0.3)

                limit = params_planet["r_ring_add_mm"] + 10 * params_planet["m_mm"]
                ax.set_xlim(-limit, limit)
                ax.set_ylim(-limit, limit)

                plt.show()

                print("\nYou have closed the planetary stage plot window.")
                change = input("Do you want to modify any planetary parameter? (y/n): ").strip().lower()
                if change != "y":
                    # Final summary
                    print("\n=== FINAL PLANETARY STAGE DESIGN ===")
                    print(f"Sun teeth (Zs): {params_planet['Zs']}")
                    print(f"Ring teeth (Zr): {params_planet['Zr']}")
                    print(f"Planet teeth (Zp): {params_planet['Zp']}")
                    print(f"Gear module: {params_planet['m_mm']} mm")
                    print(f"Reduction ratio (sun→carrier, ring fixed): {-params_planet['ratio']:.2f}")
                    print(f"Number of planets: {params_planet['N_planets']}")
                    print(f"Sun radius: {params_planet['r_sun_mm']:.2f} mm")
                    print(f"Carrier radius: {params_planet['r_carrier_mm']:.2f} mm")
                    print(f"Ring pitch radius: {params_planet['r_ring_pitch_mm']:.2f} mm")
                    print("\nPlanetary design completed successfully.\n")

                    params["planetary"] = params_planet
                    return params
                else:
                    print("From which step do you want to resume?")
                    for idx, sname in enumerate(steps, 1):
                        print(f" {idx}) {sname}")
                    pick = int(input(f"Select index (1..{len(steps)}): "))
                    step_index = pick - 1
                    params["planetary"] = params_planet
                    break  # exit plotting loop to continue outer while

    params["planetary"] = params_planet
    return params




# ======================================
# ---- OUTER CYCLOIDAL COMPACT-CAM STAGE DESIGN ----
# ======================================

def design_cycloidal_cc_stage(inner_d_ring):
    print("\n=== OUTER COMPACT-CAM CYCLOIDAL STAGE DESIGN ===\n")

    Dr = float(input("Enter roller diameter Dr [mm]: "))
    u_target = float(input("Target transmission ratio u: "))
    u_tol = float(input("Tolerance ± on ratio: "))
    Encmax = float(input("Maximum allowed radial encumbrance for the reducer [diameter, mm]: "))
    margin = float(input("Margin determining the Planetary Sun thickness [mm]: "))

    # Ensure cycloidal encloses planetary
    D = inner_d_ring + 2*margin#+ Dr

    print(f"Given the chosen parameters, the minimum internal diameter for the bearing to assemble between the sun-shaft and the cycloidal disk is: {D}")
    Dm = float(input("Provide the expected outer diameter of the mentioend bearing [mm]: "))
    Dmin = Dm + Dr + margin
    Dmax = Encmax - Dr            

    def nearest_int(x): return math.floor(x) if x - int(x) < 0.5 else math.ceil(x)
    def estimate_maxN(D, Dr): return nearest_int((D - Dr)*math.pi / Dr)

    N1max = estimate_maxN(Dmax, Dr)
    N2max = estimate_maxN(Dmax, Dr)
    #print(f"Maximum roller count allowed stage 2: {N2max}")

    Nmin = 3
    possible = []
    for N1 in range(Nmin, N1max+1):
        for N2 in range(Nmin, N2max+1):
            n1, n2 = N1-1, N2-1
            if (n1*N2 - n2*N1) == 0:
                continue
            ucalc = (n1*N2)/(n1*N2 - n2*N1)
            if (u_target - u_tol) <= abs(ucalc) <= (u_target + u_tol):
                possible.append([ucalc, N1, N2])

    if not possible:
        print("No valid roller counts found.")
        return None

    for i, sol in enumerate(possible, 1):
        print(f"{i}: u={sol[0]:.4f}, N1={sol[1]}, N2={sol[2]}")

    print(f"Maximum rollers' number allowed: {N1max}")

    choice = int(input("Select desired configuration index: ")) - 1
    ur, N1, N2 = possible[choice]

    Dmin1 = compute_Dmin(N1, Dr)
    Dmin2 = compute_Dmin(N2, Dr)
    Dmin1 = max(Dmin1, Dmin)
    Dmin2 = max(Dmin2, Dmin)
    print(f"\nAllowed pitch diameter range for stage 1: {(Dmin1):.1f} mm → {Dmax:.1f} mm")
    print(f"\nAllowed pitch diameter range for stage 2: {(Dmin2):.1f} mm → {Dmax:.1f} mm")

    D1_pitch = float(input("Enter pitch diameter Dp1 within that range [mm]: "))
    if not (Dmin1 <= D1_pitch <= Dmax):
        print("Invalid Dp1, out of range.")
        return None
    
    D2_pitch = float(input("Enter pitch diameter Dp2 within that range [mm]: "))
    if not (Dmin2 <= D2_pitch <= Dmax):
        print("Invalid Dp2, out of range.")
        return None
    
    # Dmin = inner_d_ring + 2*margin + Dr
    #dmin = inner_d_ring + 2*margin
    while True:
        ecc = float(input("Enter eccentricity [mm]: "))
        d = min(D1_pitch, D2_pitch)
        if (d - 2*ecc) > Dmin:
            break  # valid eccentricity, exit loop
        else:
            th = (d - Dmin)/2
            print(f"Invalid eccentricity. The value must must be < {th:.2f} in order to let the sun-shaft have the previously defined margin thickness.")
            print("Please enter new values for the pitch diameters and eccentricity.")
            print(f"\nAllowed pitch diameter range for stage 1: {(Dmin1):.1f} mm → {Dmax:.1f} mm")
            print(f"\nAllowed pitch diameter range for stage 2: {(Dmin2):.1f} mm → {Dmax:.1f} mm")

            D1_pitch = float(input("Enter pitch diameter Dp1 within that range [mm]: "))
            if not (Dmin1 <= D1_pitch <= Dmax):
                print("Invalid Dp1, out of range.")
                return None
            
            D2_pitch = float(input("Enter pitch diameter Dp2 within that range [mm]: "))
            if not (Dmin2 <= D2_pitch <= Dmax):
                print("Invalid Dp2, out of range.")
                return None

    R1_pitch = D1_pitch / 2
    R2_pitch = D2_pitch / 2
    roll_r = Dr / 2

    cycloid1 = make_cycloidal_profile(ecc, roll_r, N1, R1_pitch)
    cycloid2 = make_cycloidal_profile(ecc, roll_r, N2, R2_pitch)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))
    for ax, (cyc, N, R) in zip([ax1, ax2], [(cycloid1, N1, R1_pitch), (cycloid2, N2, R2_pitch)]):
        ax.set_aspect("equal")
        ax.plot(cyc[:, 0], cyc[:, 1], lw=2)
        xb, yb = make_circle(ecc, 0, Dm/2)
        ax.fill(xb, yb, "grey")
        for i in range(N):
            ang = 2*np.pi*i/N
            x, y = R*np.cos(ang), R*np.sin(ang)
            xr, yr = make_circle(x, y, roll_r)
            ax.fill(xr, yr, "yellow")
        ax.set_title(f"Stage with N={N}")

    legend = [
        Patch(facecolor='yellow', edgecolor='yellow', linestyle='-', label='Rollers'),
        Patch(facecolor='none', edgecolor='blue', linestyle='-', label='cycloidal Disk'),
        Patch(facecolor='grey', edgecolor='grey', label='Bearing')
    ]
    ax.legend(handles=legend, loc='upper right', fontsize=8)
    plt.show()

    return {"Ratio": ur, "N1": N1, "Dp1": D1_pitch, "N2": N2, "Dp2": D2_pitch, "Dr": Dr, "ecc": ecc, "marg": margin, "bearing": Dm}


# ======================================
# ---- OUTER CYCLOIDAL CLASSIC STAGE DESIGN ----
# ======================================

def design_cycloidal_cl_stage(inner_d_ring):
    print("\n=== OUTER CLASSIC CYCLOIDAL STAGE DESIGN ===\n")

    # loop principale per parametrizzazione
    while True:
        Dr = float(input("Enter roller diameter Dr [mm]: "))
        Encmax = float(input("Maximum allowed radial encumbrance for the reducer [diameter, mm]: "))
        margin = float(input("Margin determining the Planetary Sun thickness [mm]: "))

        # Ensure cycloidal encloses planetary
        D = inner_d_ring + 2*margin
        print(f"Minimum internal diameter for the bearing to assemble between the sun-shaft and the cycloidal disk: {D:.2f} mm")

        Dm = float(input("Provide the expected outer diameter of the mentioned bearing [mm]: "))
        Dmin = Dm + Dr + margin
        Dmax = Encmax - Dr
        print(f"Calculated Dmin = {Dmin:.2f} mm, Dmax = {Dmax:.2f} mm")

        # funzione per arrotondare al più vicino intero
        def custom_round(x):
            decimal_part = x - int(x)
            if decimal_part < 0.5:
                return math.floor(x)
            else:
                return math.ceil(x)

        # funzione per stimare massimo numero di rulli
        def compute_value(D, Dr):
            raw = (D / 2) * math.pi / (Dr / 2)
            return custom_round(raw)

        rounded_N = compute_value(Dmax, Dr)
        rounded_Nmin = compute_value(Dmin, Dr)
        print(f"\nApprox. maximum number of rollers: {rounded_N}")
        print(f"Thus, the maximum allowed transmission ratio is 1:{rounded_N - 1}")
        print(f"\nApprox. minimum number of rollers: {rounded_Nmin}")
        print(f"Thus, the minimum allowed transmission ratio is 1:{rounded_Nmin - 1}")

        # loop per scegliere il rapporto desiderato
        while True:
            u_target = float(input("Target transmission ratio u: [ex: for 1:30 write 30]: "))
            if (u_target <= (rounded_N - 1)) & (u_target >= (rounded_Nmin - 1)):
                N = u_target + 1
                print(f"\nTarget ratio accepted. Number of rollers N = {N}")
                break
            else:
                print(f"Target ratio too high! Maximum allowed is 1:{rounded_N - 1}")
                print(f"The minimum allowed transmission ratio is 1:{rounded_Nmin - 1}")
                modify = input("Do you want to modify Dr, Encmax, or margin? (y/n): ").strip().lower()
                if modify == 'y':
                    print("Let's re-enter the parameters.\n")
                    break  # esce dal loop del target e ritorna al loop principale per ricalcolare Dmax
                else:
                    print("Please enter a valid target ratio.")
        
        # se il target è approvato e non si vuole modificare, esci dal loop principale
        if u_target <= (rounded_N - 1):
            break

    # stampa finale dei parametri scelti
    print("\n=== FINAL PARAMETERS FOR THE CYCLIDAL STAGE ===")
    print(f"Roller diameter Dr: {Dr} mm")
    print(f"Maximum allowed encumbrance Encmax: {Encmax} mm")
    print(f"Margin: {margin} mm")
    print(f"Internal diameter D: {D:.2f} mm")
    print(f"Bearing outer diameter Dm: {Dm:.2f} mm")
    print(f"Dmin = {Dmin:.2f} mm, Dmax = {Dmax:.2f} mm")
    print(f"Approx. max number of rollers: {rounded_N}")
    print(f"Approx. min number of rollers: {rounded_Nmin}")
    print(f"Chosen target transmission ratio: 1:{u_target}")
    print(f"Number of rollers used: N = {N}")
        
    Dmin = compute_Dmin(N, Dr)
    print(f"\nAllowed pitch diameter range: {(Dmin):.1f} mm → {Dmax:.1f} mm")
    
    D_pitch = float(input("Enter pitch diameter Dp within that range [mm]: "))
    if not (Dmin <= D_pitch <= Dmax):
        print("Invalid Dp, out of range.")
        return None
    
    # Dmin = inner_d_ring + 2*margin + Dr
    #dmin = inner_d_ring + 2*margin
    finish = False
    while finish == False:
        while True:
            ecc = float(input("Enter eccentricity [mm]: "))
            d = D_pitch
            if (d - 2*ecc) > Dmin:
                break  # valid eccentricity, exit loop
            else:
                th = (d - Dmin)/2
                print(f"Invalid eccentricity. The value must must be < {th:.2f} in order to let the sun-shaft have the previously defined margin thickness.")
                print("Please enter new values for the pitch diameters and eccentricity.")
                print(f"\nAllowed pitch diameter range: {(Dmin):.1f} mm → {Dmax:.1f} mm")
                
                D_pitch = float(input("Enter pitch diameter Dp within that range [mm]: "))
                if not (Dmin <= D_pitch <= Dmax):
                    print("Invalid Dp, out of range.")
                    return None
        max_pin_D = (((D_pitch - Dr - 2*ecc) - Dm)/2) - 2*ecc
        print(f"Given the chosen parameters, the maximum allowed diameter for the cycloidal carrier pins is {max_pin_D:.2f} mm.")
        pin = float(input("Enter the employed pins diameter for the cycloidal carrier [mm]: "))
        if not (pin <= max_pin_D):
            print("Invalid Dp, out of range.")
            return None
        RH = pin/2 + ecc
        print(f"The resulting diameter for the holes in the cycloidal disk is {RH*2:.2f} mm.")
        Dc = (((D_pitch - Dr - 2*ecc) + Dm)/2)
        print(f"The diameter for the carrier pins is {Dc:.2f} mm.")
        Rc = Dc/2
        ###
        R_pitch = D_pitch / 2
        roll_r = Dr / 2

        cycloid = make_cycloidal_profile(ecc, roll_r, N, R_pitch)

        fig, ax = plt.subplots(figsize=(8, 8))  # un solo subplot
        cyc = cycloid
        R = R_pitch

        # Se cyc è Nx2 array di punti (x, y)
        ax.plot(cyc[:, 0], cyc[:, 1], lw = 2)

        ax.set_title("Cycloidal profile")
        ax.set_aspect("equal")
        ax.grid(True)
        ax.set_xlabel("X [mm]")
        ax.set_ylabel("Y [mm]")

        xb, yb = make_circle(ecc, 0, Dm/2)
        ax.fill(xb, yb, "grey")
        for i in range(int(N)):
            ang = 2*np.pi*i/N
            x, y = R*np.cos(ang), R*np.sin(ang)
            xr, yr = make_circle(x, y, roll_r)
            ax.fill(xr, yr, "yellow")
        for i in range(int(4)):
            ang = 2*np.pi*i/4
            x, y = ecc + Rc*np.cos(ang), Rc*np.sin(ang)
            xr, yr = make_circle(x, y, RH)
            ax.fill(xr, yr, "orange")
        ax.set_title(f"Stage with N={N}")

        legend = [
            Patch(facecolor='yellow', edgecolor='yellow', linestyle='-', label='Rollers'),
            Patch(facecolor='none', edgecolor='blue', linestyle='-', label='cycloidal Disk'),
            Patch(facecolor='grey', edgecolor='grey', label='Bearing'),
            Patch(facecolor='orange', edgecolor='orange', label='Carrier pin Holes')]
        ax.legend(handles=legend, loc='upper right', fontsize=8)
        plt.show()
        Dp = D_pitch
        ###
        ok = input("Are these parameters acceptable? (y/n): ").lower()
        if ok == "y":

            finish = True
        elif ok == "n":
            continue

    return {"Ratio": u_target, "N": N, "Dp": D_pitch, "N2": N, "Dp2": D_pitch, "Dr": Dr, "ecc": ecc, "marg": margin, "bearing": Dm, "Rc": Rc, "RH": RH}

# ======================================
# ---- INTEGRATION AND VISUALIZATION ----
# ======================================

# Final summary
#                    print("\n=== FINAL PLANETARY STAGE DESIGN ===")
#                    print(f"Sun teeth (Zs): {params_planet['Zs']}")
#                    print(f"Ring teeth (Zr): {params_planet['Zr']}")
#                    print(f"Planet teeth (Zp): {params_planet['Zp']}")
#                    print(f"Gear module: {params_planet['m_mm']} mm")
#                    print(f"Reduction ratio (sun→carrier, ring fixed): {-params_planet['ratio']:.2f}")
#                    print(f"Number of planets: {params_planet['N_planets']}")
#                    print(f"Sun radius: {params_planet['r_sun_mm']:.2f} mm")
#                    print(f"Carrier radius: {params_planet['r_carrier_mm']:.2f} mm")
#                    print(f"Ring pitch radius: {params_planet['r_ring_pitch_mm']:.2f} mm")
#                    print("\nPlanetary design completed successfully.\n")

def combine_stages(inner_planet, outer_cycloid, variant):
    print("\n=== INTEGRATED DESIGN VISUALIZATION ===\n")

    Zs, Zr, m = inner_planet["planetary"]["Zs"], inner_planet["planetary"]["Zr"], inner_planet["planetary"]["m_mm"]
    Np, rc, up = inner_planet["planetary"]["N_planets"], inner_planet["planetary"]["r_carrier_mm"], inner_planet["planetary"]["ratio"]
    margin = outer_cycloid["marg"]
    uc = outer_cycloid["Ratio"]
    Db = outer_cycloid["bearing"]
    if variant == 1:
        Rc = outer_cycloid["Rc"]
        RH = outer_cycloid["RH"]
    total_ratio = -up * uc
    CV = variant

    print(f"→ Planetary ratio = {-up:.2f}")
    print(f"→ Cycloidal ratio = {uc:.2f}")
    print(f"→ TOTAL REDUCTION = {total_ratio:.2f}\n")

    def plot_planet(ax):
        r_sun = Zs*m/2
        r_ring = Zr*m/2
        r_planet = ((Zr - Zs)/2)*m/2

        ax.set_aspect("equal")
        xs, ys = make_circle(0, 0, r_sun)
        ax.plot(xs, ys, "red", linestyle='--')
        xr, yr = make_circle(0, 0, r_ring)
        xrs, yrs = make_circle(ecc, 0, r_ring + margin)
        ax.plot(xr, yr, "k", linestyle='--')
        ax.plot(xrs, yrs, "k")
        for i in range(Np):
            ang = 2*np.pi*i/Np
            xc, yc = rc*np.cos(ang), rc*np.sin(ang)
            xp, yp = make_circle(xc, yc, r_planet)
            ax.plot(xp, yp, "g", linestyle='--')

    if CV == 1:
        R_pitch = outer_cycloid["Dp"]/2
        roll_r = outer_cycloid["Dr"]/2
        ecc = outer_cycloid["ecc"]
        N = outer_cycloid["N"]
        cyc = make_cycloidal_profile(ecc, roll_r, N, R_pitch)
        R = R_pitch
        fig, ax = plt.subplots(figsize=(14, 7))        
        plot_planet(ax)
        ax.plot(cyc[:, 0], cyc[:, 1], "blue", lw=2)
        circle = plt.Circle((ecc, 0), Db/2, fill=False, color='grey', linestyle='--', linewidth=1.5)
        ax.add_artist(circle)
        for i in range(int(N)):
            ang = 2*np.pi*i/N
            x, y = R*np.cos(ang), R*np.sin(ang)
            xr, yr = make_circle(x, y, roll_r)
            ax.fill(xr, yr, "yellow", alpha=0.6)
        ax.set_title("Planetary inside Cycloidal")
        
        # Nc, Rc, RH
        for i in range(int(4)):
            ang = 2*np.pi*i/4
            xc, yc = ecc + Rc*np.cos(ang), Rc*np.sin(ang)
            circle = plt.Circle((xc, yc), RH, fill=False, color='orange', linestyle='--', linewidth=1.5)
            ax.add_artist(circle)

        legend = [
            Patch(facecolor='none', edgecolor='red', linestyle='--', label='Sun gear (primitive circle)'),
            Patch(facecolor='none', edgecolor='black', linestyle='--', label='Ring gear (primitive circle)'),
            Patch(facecolor='none', edgecolor='black', label='Eccentric Shaft'),
            Patch(facecolor='none', edgecolor='grey', linestyle='--', label='Bearing outer diameter'),
            Patch(facecolor='none', edgecolor='green', linestyle='--', label='Planet gears (primitive circles)'),
            Line2D([], [], color='blue', lw=2, label='Cycloidal profile'),
            Patch(facecolor='none', edgecolor='orange', linestyle='--', label='Disk Holes for Carrier pins'),
            Patch(facecolor='yellow', label='Rollers')]
        ax.legend(handles=legend, loc='upper right', fontsize=8)

    elif CV == 2:
            R1_pitch = outer_cycloid["Dp1"]/2
            R2_pitch = outer_cycloid["Dp2"]/2
            roll_r = outer_cycloid["Dr"]/2
            ecc = outer_cycloid["ecc"]
            N1, N2 = outer_cycloid["N1"], outer_cycloid["N2"]
            cyc1 = make_cycloidal_profile(ecc, roll_r, N1, R1_pitch)
            cyc2 = make_cycloidal_profile(ecc, roll_r, N2, R2_pitch)

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
            for ax, cyc, N, stg, R in zip([ax1, ax2], [cyc1, cyc2], [N1, N2], ["Stage 1", "Stage 2"], [R1_pitch, R2_pitch]):
                plot_planet(ax)
                ax.plot(cyc[:, 0], cyc[:, 1], "blue", lw=2)
                circle = plt.Circle((ecc, 0), Db/2, fill=False, color='grey', linestyle='--', linewidth=1.5)
                ax.add_artist(circle)
                for i in range(int(N)):
                    ang = 2*np.pi*i/N
                    x, y = R*np.cos(ang), R*np.sin(ang)
                    xr, yr = make_circle(x, y, roll_r)
                    ax.fill(xr, yr, "yellow", alpha=0.6)
                ax.set_title(f"{stg}: Planetary inside Cycloidal")

            legend = [
                Patch(facecolor='none', edgecolor='red', linestyle='--', label='Sun gear (primitive circle)'),
                Patch(facecolor='none', edgecolor='black', linestyle='--', label='Ring gear (primitive circle)'),
                Patch(facecolor='none', edgecolor='black', label='Eccentric Shaft'),
                Patch(facecolor='none', edgecolor='grey', linestyle='--', label='Bearing outer diameter'),
                Patch(facecolor='none', edgecolor='green', linestyle='--', label='Planet gears (primitive circles)'),
                Line2D([], [], color='blue', lw=2, label='Cycloidal profile'),
                Patch(facecolor='yellow', label='Rollers')
            ]
            ax1.legend(handles=legend, loc='upper right', fontsize=8)
            ax2.legend(handles=legend, loc='upper right', fontsize=8)

    plt.show()

# ======================================
# ---- MAIN ----
# ======================================

if __name__ == "__main__":
    planetary = design_planetary_stage()
    if planetary:
        proceed = input("\nProceed to outer cycloidal stage design? (y/n): ").lower()
        if proceed == "y":
            print("\nWhich cycloidal architecture do you want to use for the external stage?: ")
            CYC_v = float(input("[1 - Classic cycloidal] or [2 - Compact-Cam]: "))
            while True:
                if CYC_v == 1 or CYC_v == 2:
                    if CYC_v == 1:
                        cycloidal = design_cycloidal_cl_stage(inner_d_ring=planetary["planetary"]["Zr"] * planetary["planetary"]["m_mm"] + 2 * planetary["planetary"]["m_mm"])
                        break
                    elif CYC_v == 2:
                        cycloidal = design_cycloidal_cc_stage(inner_d_ring=planetary["planetary"]["Zr"] * planetary["planetary"]["m_mm"] + 2 * planetary["planetary"]["m_mm"])
                        break
                else:
                    CYC_v = float(input("[1 - Classic cycloidal] or [2 - Compact-Cam]: "))


        # === FINAL REPORT ===
        print("\n" + "="*50)
        print("        FINAL DESIGN PARAMETERS SUMMARY")
        print("="*50 + "\n")

        if "planetary" in planetary:
            print(" PLANETARY STAGE")
            p = planetary["planetary"]
            print(f"  Sun teeth (Zs): {p['Zs']}")
            print(f"  Ring teeth (Zr): {p['Zr']}")
            print(f"  Planet teeth (Zp): {p['Zp']}")
            print(f"  Module (m): {p['m_mm']} mm")
            print(f"  Reduction ratio (sun→carrier, ring fixed): {-p['ratio']:.2f}")
            print(f"  Number of planets: {p['N_planets']}")
            print(f"  Carrier radius: {p['r_carrier_mm']:.2f} mm")
            print("-"*50)

        if cycloidal:
            if CYC_v == 1:
                print("\n CLASSIC CYCLOIDAL STAGE")
                print(f"  Roller diameter Dr: {cycloidal['Dr']} mm")
                print(f"  Bearing outer diameter Dm: {cycloidal['bearing']} mm")
                print(f"  Eccentricity: {cycloidal['ecc']} mm")
                print(f"  Pitch diameter: Dp={cycloidal['Dp']} mm")
                print(f"  Number of rollers: N={cycloidal['N']}")
                print(f"  Transmission ratio:{cycloidal['Ratio']}")
            elif CYC_v == 2:
                print("\n COMPACT-CAM CYCLOIDAL STAGE")
                print(f"  Roller diameter Dr: {cycloidal['Dr']} mm")
                print(f"  Bearing outer diameter Dm: {cycloidal['bearing']} mm")
                print(f"  Eccentricity: {cycloidal['ecc']} mm")
                print(f"  Pitch diameter(s): Dp1={cycloidal['Dp1']} mm, Dp2={cycloidal['Dp2']} mm")
                print(f"  Number of rollers: N1={cycloidal['N1']}, N2={cycloidal['N2']}")
                print(f"  Transmission ratio:{cycloidal['Ratio']}")
            print("-"*50)

            combine_stages(planetary, cycloidal, CYC_v)

        print("\nDESIGN COMPLETE. All data summarized above.\n")