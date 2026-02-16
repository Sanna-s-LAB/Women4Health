#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import sys
from pathlib import Path
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

# ======================== CLI ========================
if len(sys.argv) != 4:
    print("Usage: python3 script.py <input_csv> <output_csv> <log_txt>")
    sys.exit(1)

input_path = sys.argv[1]
output_path = sys.argv[2]
log_path = sys.argv[3]

log = []
def logh(msg): log.append(msg)

# ======================== Load & basic filtering ========================
df = pd.read_csv(input_path)

# Hard-removal of specific codes (as in your original logic)
df = df[~df['Code'].isin(['S005_4','X012_4','X017_4','X075_4'])]

req_cols = ['Code','ID','Visit_number','cycle_index','BES','PROG','FSH','LH']
miss = [c for c in req_cols if c not in df.columns]
if miss:
    raise ValueError(f"Missing columns: {', '.join(miss)}")

df = df[req_cols].copy()

# Coercions to numeric / integer
for col in ['BES','PROG','FSH','LH']:
    df[col] = pd.to_numeric(df[col], errors='coerce')
df['Visit_number'] = pd.to_numeric(df['Visit_number'], errors='coerce').astype('Int64')
df['cycle_index']  = pd.to_numeric(df['cycle_index'],  errors='coerce').astype('Int64')

before = len(df)
df = df.dropna(subset=req_cols).reset_index(drop=True)
after = len(df)
logh(f"📥 Loaded rows: {before} → kept {after} after NA/drop filters.")

hormone_cols = ['FSH','LH','PROG','BES']

# ======================== Helper functions ========================
def get_phase_mapping(lh_peak_visit):
    """
    Mapping Visit_number -> Phase for 'regular-like' patterns based on LH peak position.

    The idea:
    - If LH peaks at visit 1:
        We interpret that the whole cycle is shifted "early".
        Mapping (Visit → Phase) is:
            1 → 2
            2 → 3
            3 → 4
            4 → NA (to be resolved with a trajectory-based rule).
    - If LH peaks at visit 3:
        We interpret that the whole cycle is shifted "late".
        Mapping:
            1 → NA (to be resolved)
            2 → 1
            3 → 2
            4 → 3
    - If LH peaks at 4:
        We treat the subject as irregular (no mapping returned).
    """
    if lh_peak_visit == 1:
        return {1:2, 2:3, 3:4, 4:pd.NA}
    elif lh_peak_visit == 3:
        return {1:pd.NA, 2:1, 3:2, 4:3}
    return None  # peak in 4 (or other unexpected cases) → no "regular" mapping

def mean_euclidean(X, Y):
    """
    Mean Euclidean distance between corresponding rows of X and Y.

    If X and Y are arrays of shape (n, d), this returns:
        (1/n) * sum_i ||X[i] - Y[i]||_2
    """
    X = np.asarray(X); Y = np.asarray(Y)
    return float(np.mean(np.linalg.norm(X - Y, axis=1)))

def save_confusion_matrix_png(df_res, out_csv_path):
    """
    Build and save the confusion matrix between Phase and Visit_number.

    Convention:
    - rows    = Phase (1..4)
    - columns = Visit_number (1..4)

    This allows us to visually see how visits are redistributed into phases.
    """
    cm = pd.crosstab(
        df_res['Phase'],
        df_res['Visit_number'],
        dropna=False
    ).reindex(index=[1,2,3,4], columns=[1,2,3,4], fill_value=0)

    png_path = Path(out_csv_path).with_suffix('.png')

    fig, ax = plt.subplots(figsize=(6,5), dpi=150)
    im = ax.imshow(cm.values, cmap='viridis', aspect='auto')

    # Cell counts
    for i in range(4):
        for j in range(4):
            ax.text(j, i, str(cm.values[i, j]), ha='center', va='center')

    ax.set_xticks(range(4)); ax.set_xticklabels([1,2,3,4])
    ax.set_yticks(range(4)); ax.set_yticklabels([1,2,3,4])
    ax.set_xlabel("Visit_number")
    ax.set_ylabel("Phase")
    ax.set_title("Confusion Matrix: Phase × Visit_number")
    plt.colorbar(im, ax=ax, label="Count")
    fig.tight_layout()

    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, bbox_inches='tight')
    plt.close(fig)
    logh(f"🧮 Confusion matrix saved: {png_path}")
    return cm

# ======================== Split regular / irregular ========================
# "Regular" women: only 1 cycle_index → assumed one full 4-visit cycle.
# "Irregular" women: more than 1 cycle_index (multiple cycles / misaligned cycles).
n_cycles = df.groupby('ID')['cycle_index'].nunique()
regular_ids   = set(n_cycles[n_cycles == 1].index)
irregular_ids = set(n_cycles[n_cycles > 1].index)

df_reg = df[df['ID'].isin(regular_ids)].copy()
df_irr = df[df['ID'].isin(irregular_ids)].copy()

logh(f"👥 Regular women (1 cycle): {len(regular_ids)}")
logh(f"👥 Irregular women (>1 cycle): {len(irregular_ids)}")

# ======================== StandardScaler (fit on regulars) ========================
# We standardize hormones using only regular women, so that:
# - all "template" trajectories live in this standardized space,
# - individual subjects (regular and irregular) can be compared consistently to them.
scaler = StandardScaler().fit(df_reg[hormone_cols])
means = dict(zip(hormone_cols, scaler.mean_))
stds  = dict(zip(hormone_cols, scaler.scale_))
logh("🧮 StandardScaler (regulars): " + ", ".join(
    [f"{k}: mean={means[k]:.3f}, sd={stds[k]:.3f}" for k in hormone_cols]
))

# ======================== TM23 / TM34 (standardized) for rule 1B ========================
# Here we build *visit-based* trajectories for regular women in standardized space:
#   TM23[h] = (mean_h_visit2_std, mean_h_visit3_std)
#   TM34[h] = (mean_h_visit3_std, mean_h_visit4_std)
#
# These represent average hormone trajectories between visits 2→3 and 3→4.
# They are used in the "local" rule for regular women with LH peak at visit 1
# to decide whether visit 4 should be mapped to phase 3 or 4, by comparing
# the subject's T34 to TM23 and TM34.
df_reg_std_TM = df_reg.copy()
df_reg_std_TM[hormone_cols] = scaler.transform(df_reg_std_TM[hormone_cols])

TM23 = {}
TM34 = {}
for h in hormone_cols:
    v2 = df_reg_std_TM.loc[df_reg_std_TM['Visit_number'] == 2, h]
    v3 = df_reg_std_TM.loc[df_reg_std_TM['Visit_number'] == 3, h]
    v4 = df_reg_std_TM.loc[df_reg_std_TM['Visit_number'] == 4, h]

    TM23[h] = None
    TM34[h] = None

    if len(v2) > 0 and len(v3) > 0:
        TM23[h] = np.array([v2.mean(), v3.mean()], dtype=float)
    if len(v3) > 0 and len(v4) > 0:
        TM34[h] = np.array([v3.mean(), v4.mean()], dtype=float)

used_23 = [h for h in hormone_cols if TM23[h] is not None]
used_34 = [h for h in hormone_cols if TM34[h] is not None]
logh(f"🧭 TM23 (standardized) available for: {', '.join(used_23) if used_23 else 'none'}")
logh(f"🧭 TM34 (standardized) available for: {', '.join(used_34) if used_34 else 'none'}")

def traj_distance_T34_vs_TM(row3_std, row4_std, TM):
    """
    Distance between the subject's T34 and a template TMxx in standardized space.

    For each hormone h:
        - Subject's trajectory T34_h = (h3_std, h4_std).
        - Template trajectory TM[h] = (mean_h_visit_a_std, mean_h_visit_b_std)
          where (a,b) is (2,3) for TM23 or (3,4) for TM34.

    We compute:
        total_distance = sum_h || T34_h - TM[h] ||_2

    and also count how many hormones had a non-None template (used).

    This gives a *trajectory-based* distance between:
        - subject's visit3→visit4 profile, and
        - average visit2→3 or 3→4 profile across regular women.
    """
    total = 0.0
    used = 0
    for h in hormone_cols:
        tm_vec = TM.get(h)
        if tm_vec is None:
            continue
        t34_vec = np.array([row3_std[h], row4_std[h]], dtype=float)
        total += np.linalg.norm(t34_vec - tm_vec)
        used += 1
    return (total, used)

# ======================== Mapping for regulars ========================
mapped = []
moved_to_irregular_ids = set()

for wid, g in df_reg.groupby('ID'):
    g = g.sort_values('Visit_number').copy()

    # LH peak: find the visit with maximum LH value
    peak_idx = g['LH'].idxmax()
    lh_peak_visit = int(g.loc[peak_idx, 'Visit_number'])
    lh_peak_val   = float(g.loc[peak_idx, 'LH'])

    # If LH peaks at visit 4 → treat this subject as irregular.
    if lh_peak_visit == 4:
        moved_to_irregular_ids.add(wid)
        df_irr = pd.concat([df_irr, g], ignore_index=True)
        logh(f"➡️  REG {wid}: LHmax={lh_peak_val:.4f} at visit 4 → moved to IRREGULAR (peak in 4).")
        continue

    mapping = get_phase_mapping(lh_peak_visit)

    if mapping is None:
        # No "regular" mapping: keep Phase = Visit_number
        g['Phase'] = g['Visit_number']
        logh(
            f"ℹ️ REG {wid}: LHmax={lh_peak_val:.4f} at visit {lh_peak_visit}; "
            f"no mapping available → Phase=Visit."
        )
    else:
        # Apply the Visit→Phase mapping based on LH peak (1-step shift model)
        g['Phase'] = g['Visit_number'].map(mapping)
        logh(
            f"🌀 REG {wid}: LHmax={lh_peak_val:.4f} at visit {lh_peak_visit} "
            f"→ applied phase mapping."
        )

    g['Phase'] = g['Phase'].astype('Int64')

    # --- Handle NAs introduced by mapping (only occur for peak 1 or 3) ---
    if g['Phase'].isna().any():
        na_visits = g.loc[g['Phase'].isna(), 'Visit_number'].astype(int).tolist()

        if lh_peak_visit == 3:
            # Rule 1A: peak at visit 3
            # Typically visit 1 is unmapped; we assign Phase=1 to the NA.
            g.loc[g['Phase'].isna(), 'Phase'] = 1
            logh(f"🔧 REG {wid}: NA at visits {na_visits} (peak=3) → NA→Phase=1 (rule 1A).")

        elif lh_peak_visit == 1:
            # Rule 1B (more complex, trajectory-based):
            #
            # We must decide how to treat visit 4 (mapped as NA by {1:2,2:3,3:4,4:NA}).
            # Intuition:
            #   - If visit 3 is missing but visit 4 exists, we keep visit 4 as Phase=4.
            #   - If both visit 3 and 4 exist, we compare the subject's T34
            #     to TM23 and TM34 in standardized space:
            #         if T34 is closer to TM23 → Phase(visit4) = 3
            #         else                     → Phase(visit4) = 4
            visits_num = g['Visit_number'].astype(int).tolist()
            has_v3 = 3 in visits_num
            has_v4 = 4 in visits_num

            if not has_v4:
                # Theoretically unexpected: NA but no visit 4 row.
                # We keep a defensive fallback Phase=4 for all NA.
                g.loc[g['Phase'].isna(), 'Phase'] = 4
                logh(
                    f"🔧 REG {wid}: NA at visits {na_visits} (peak=1) but visit 4 missing "
                    f"(unexpected) → NA→Phase=4 fallback."
                )

            elif not has_v3:
                # If visit 3 is missing but visit 4 exists:
                #   we have no T34 trajectory; we keep visit 4 as Phase=4.
                g.loc[g['Phase'].isna(), 'Phase'] = 4
                logh(
                    f"🔧 REG {wid}: NA at visits {na_visits} (peak=1) and visit 3 missing "
                    f"→ keep visit 4 as Phase=4."
                )

            else:
                # Normal case: both visit 3 and 4 exist.
                # We build the subject's T34 in standardized space and
                # compare it to TM23 and TM34.
                v3 = g[g['Visit_number'] == 3][hormone_cols]
                v4 = g[g['Visit_number'] == 4][hormone_cols]

                if len(v3) == 1 and len(v4) == 1:
                    # Standardize visit 3 and 4 values with the same scaler
                    v3_std = scaler.transform(v3[hormone_cols])
                    v4_std = scaler.transform(v4[hormone_cols])
                    row3_std = pd.Series(v3_std[0], index=hormone_cols)
                    row4_std = pd.Series(v4_std[0], index=hormone_cols)

                    d23, used23 = traj_distance_T34_vs_TM(row3_std, row4_std, TM23)
                    d34, used34 = traj_distance_T34_vs_TM(row3_std, row4_std, TM34)

                    if used23 == 0 or used34 == 0:
                        # If no template is available → fallback Phase=4
                        chosen = 4
                        g.loc[g['Phase'].isna(), 'Phase'] = chosen
                        logh(
                            f"🔧 REG {wid}: NA at visits {na_visits} (peak=1) but no usable "
                            f"TM23/TM34 → NA→Phase={chosen} fallback."
                        )
                    else:
                        # Trajectory-based decision:
                        # Compare subject T34 to TM23 vs TM34
                        #   if d(T34, TM23) < d(T34, TM34) → Phase(4)=3
                        #   else                          → Phase(4)=4
                        chosen = 3 if d23 < d34 else 4
                        g.loc[g['Phase'].isna(), 'Phase'] = chosen
                        logh(
                            f"🔧 REG {wid}: NA at visits {na_visits} (peak=1) → "
                            f"d(T34,TM23)={d23:.3f}, d(T34,TM34)={d34:.3f} "
                            f"→ NA→Phase={chosen} (rule 1B, standardized trajectories)."
                        )
                else:
                    # Unexpected situation (multiple rows per visit, etc.) → fallback Phase=4
                    g.loc[g['Phase'].isna(), 'Phase'] = 4
                    logh(
                        f"🔧 REG {wid}: NA at visits {na_visits} (peak=1) but unexpected "
                        f"rows for visit 3/4 → NA→Phase=4 fallback."
                    )

    mapped.append(g)

df_reg = pd.concat(mapped, ignore_index=True) if mapped else df_reg.iloc[0:0, :]
df_reg['Phase'] = df_reg['Phase'].astype('Int64')
if moved_to_irregular_ids:
    logh(
        f"↪️ Subjects moved to IRREGULAR (peak in 4): "
        f"{len(moved_to_irregular_ids)} ({', '.join(sorted(map(str, moved_to_irregular_ids)))})"
    )

# ======================== Phase-based templates on regulars ========================
# Now that each regular subject has a Phase label, we:
# 1) Standardize them with the same scaler
# 2) Build:
#    - pair_templates[(1,2)], (2,3), (3,4):
#        average segments (Phase p → Phase p+1)
#    - trajectory_template:
#        average point per Phase (Phase 1..4)
reg_std = df_reg.copy()
reg_std[hormone_cols] = scaler.transform(reg_std[hormone_cols])

# --- Pair templates: segments (Phase p → Phase p+1) ---
pair_bags = {(1,2):[], (2,3):[], (3,4):[]}
for wid, g in reg_std.groupby('ID'):
    for p in [1,2,3]:
        a = g[g['Phase'] == p]
        b = g[g['Phase'] == p+1]
        if len(a) == 1 and len(b) == 1:
            # For each subject, create a 2×4 matrix:
            #   [ hormones at Phase p
            #     hormones at Phase p+1 ]
            pair_bags[(p,p+1)].append(
                np.vstack([a[hormone_cols].to_numpy(), b[hormone_cols].to_numpy()])
            )

pair_templates = {}
for key, bag in pair_bags.items():
    if bag:
        # bag is a list of matrices of shape (2,4).
        # Stacking gives (N,2,4), we average over N → (2,4) template.
        M = np.stack(bag)     # (N, 2, 4)
        tmpl = M.mean(axis=0) # (2, 4)
        pair_templates[key] = tmpl
        logh(f"📚 Visit-pair template {key} built with N={len(bag)} pairs.")
    else:
        logh(f"⚠️ No data to build visit-pair template {key}.")

# --- Phase trajectory template: average point per Phase ---
# This is a 4×4 matrix where each row is the mean standardized hormone vector
# for Phase 1,2,3,4. It represents the "typical" phase trajectory of regular women.
trajectory_template = reg_std.groupby('Phase')[hormone_cols].mean().sort_index()

# ======================== Irregulars: shift-based alignment & boundary fixes ========================
# For irregular women, we:
# 1) Standardize their hormone values with the same scaler.
# 2) For each subject, consider shifts s ∈ {-1,0,1} applied to the sequence of visits.
# 3) For each shift, compute the average distance between subject's points and
#    the phase trajectory template.
# 4) Choose the shift that minimizes this distance.
# 5) Assign Phase = Visit_number + best_shift.
# 6) Apply boundary fixes when we get Phase=0 or 5, using pair_templates and
#    the same segment-based logic used for regulars with NA.
df_irr_std = df_irr.copy()
if not df_irr_std.empty:
    df_irr_std[hormone_cols] = scaler.transform(df_irr_std[hormone_cols])

    for wid, g in df_irr_std.groupby('ID'):
        g = g.sort_values('Visit_number').copy()
        X = g[hormone_cols].to_numpy()                # (n_visits, 4), standardized
        visits = g['Visit_number'].astype(int).tolist()

        # --- Step 2–4: choose the best shift s ∈ {-1,0,1} ---
        # We define for each shift s a set of valid indices I_s such that
        #   0 <= i+s < 4  (since trajectory_template has 4 phases)
        # and we compute:
        #   D(s) = mean_i∈I_s || X[i] - trajectory_template[phase_index=i+s] ||
        best_shift, best_dist = None, np.inf
        for shift in (-1, 0, 1):
            idxs = np.arange(len(X))
            valid = (idxs + shift >= 0) & (idxs + shift < len(trajectory_template))
            if not np.any(valid):
                continue

            ref = trajectory_template.iloc[idxs[valid] + shift, :].to_numpy()
            d = mean_euclidean(X[valid], ref)
            if d < best_dist:
                best_dist, best_shift = d, shift

        if best_shift is None:
            # No valid shift (very unlikely with 1–4 visits) → fallback Phase=Visit
            phases = g['Visit_number'].astype(int).to_numpy()
            logh(f"🔁 IRR {wid}: no valid shift found → Phase=Visit.")
        else:
            # Step 5: initial Phase assignment as Visit_number + best_shift
            phases = (g['Visit_number'].astype(int) + best_shift).to_numpy()
            logh(
                f"🔗 IRR {wid}: chosen shift={best_shift:+d}, "
                f"dist={best_dist:.4f} (raw phases={phases.tolist()})."
            )

            # --- Step 6: boundary fixes ---
            if best_shift == +1:
                # With shift +1:
                #   Visit 1 → Phase 2
                #   Visit 2 → Phase 3
                #   Visit 3 → Phase 4
                #   Visit 4 → Phase 5 (out of range)
                #
                # For Phase=5, we want to decide if visit 4 is better interpreted as Phase 3 or 4.
                if 3 in visits and 4 in visits:
                    i3 = visits.index(3)
                    i4 = visits.index(4)

                    # Subject's segment (visit3 → visit4) in standardized space
                    subj_pair = np.vstack([X[i3], X[i4]])  # shape (2,4)

                    # Compare this segment to the average Phase2→3 and Phase3→4 segments
                    d_23 = mean_euclidean(subj_pair, pair_templates.get((2,3), subj_pair * 0))
                    d_34 = mean_euclidean(subj_pair, pair_templates.get((3,4), subj_pair * 0))

                    # If subject's 3→4 is closer to template 2→3, interpret visit4 as Phase 3.
                    # Otherwise, interpret as Phase 4.
                    chosen = 3 if d_23 < d_34 else 4
                    phases[i4] = chosen
                    logh(
                        f"   ↪ boundary fix (+1): compare (v3,v4) → "
                        f"d23={d_23:.3f}, d34={d_34:.3f} → set Phase(visit4)={chosen}."
                    )
                else:
                    # If we don't have both visits 3 and 4, we cannot build a segment.
                    # Fallback: keep visit 4 as Phase 4.
                    if 4 in visits:
                        i4 = visits.index(4)
                        phases[i4] = 4
                        logh(
                            "   ↪ boundary fix (+1): (v3 or v4 missing) "
                            "→ fallback Phase(visit4)=4."
                        )

            elif best_shift == -1:
                # With shift -1:
                #   Visit 1 → Phase 0 (out of range)
                #   Visit 2 → Phase 1
                #   Visit 3 → Phase 2
                #   Visit 4 → Phase 3
                #
                # We force visit1 to at least Phase 1.
                if 1 in visits:
                    i1 = visits.index(1)
                    phases[i1] = 1
                    logh("   ↪ boundary fix (−1): set Phase(visit1)=1.")

            # Finally, clamp any residual value outside [1,4] into [1,4].
            phases = np.where(phases < 1, 1, np.where(phases > 4, 4, phases))

        # Save phases back (in the original non-standardized df_irr)
        g['Phase'] = phases.astype('int64')
        df_irr.loc[g.index, 'Phase'] = g['Phase'].values

else:
    logh("ℹ️ No IRREGULAR rows to reclassify.")

# ======================== Combine regulars + irregulars & metrics ========================
df_all = pd.concat([df_reg, df_irr], ignore_index=True)\
           .sort_values(['ID','Visit_number'])\
           .reset_index(drop=True)
df_all['Phase'] = df_all['Phase'].astype('Int64')

cm = save_confusion_matrix_png(df_all[['ID','Visit_number','Phase']], output_path)

diag = np.diag(cm.values)
row_sums = cm.sum(axis=1).values
diag_acc = float(np.nanmean(np.where(row_sums > 0, diag / row_sums, np.nan)))
avg_shift = float(np.nanmean(np.abs(df_all['Phase'] - df_all['Visit_number'])))

logh("=== Agreement Metrics ===")
logh(f"Diagonal accuracy: {diag_acc*100:.2f}%")
logh(f"Mean absolute shift: {avg_shift:.2f}")
logh("Counts per Phase row: " + ", ".join([f"Phase {i+1}={int(row_sums[i])}" for i in range(4)]))
logh("Confusion matrix counts (rows=Phase, cols=Visit_number):")
for i in range(4):
    logh(f"Phase {i+1}: " + " ".join([f"{int(cm.values[i,j])}" for j in range(4)]))

# ======================== Export ========================
Path(output_path).parent.mkdir(parents=True, exist_ok=True)
df_all[['Code','Phase']].to_csv(output_path, index=False)
logh(f"💾 Output CSV saved: {output_path}")

with open(log_path, 'w') as f:
    for line in log:
        f.write(line + "\n")

print(f"Done. Output: {output_path} | Log: {log_path}")
