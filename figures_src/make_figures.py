"""Figures for Jenkins and Graef, Genomic envirotyping (Theoretical and Applied Genetics).

Run from the repository root after run_all.R and R/16_figure_data.R:

    python figures_src/make_figures.py

Needs pandas, numpy, matplotlib, cartopy and pyshp. Natural Earth boundaries
(states and provinces, lakes, US counties) are downloaded by cartopy on first use.
Writes figures/Fig1.tif to Fig5.tif and FigS1.tif ... FigS3.tif (300 dpi, LZW),
with PDF and PNG copies.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf
import cartopy.io.shapereader as shpreader
import shapefile
import shapely.geometry as sg

ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
OUT = os.path.join(ROOT, "out")
FD = os.path.join(OUT, "figure_data")
FIGDIR = os.path.join(ROOT, "figures")
os.makedirs(FIGDIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42,
})
MGS = ["0", "I", "II", "III", "IV"]
MG4 = ["I", "II", "III", "IV"]
PAL = {1: "#0072B2", 2: "#E69F00", 3: "#009E73", 4: "#CC79A7"}
GEN, ENV, GEO, DARK = "#0072B2", "#E69F00", "#7f7f7f", "#3a3a3a"
KGHA = "kg ha$^{-1}$"
PROJ = ccrs.LambertConformal(central_longitude=-89, central_latitude=42)
PC = ccrs.PlateCarree()


def save(fig, name):
    for ext, kw in [("tif", {"pil_kwargs": {"compression": "tiff_lzw"}}), ("pdf", {}), ("png", {})]:
        fig.savefig(os.path.join(FIGDIR, f"{name}.{ext}"), dpi=300, **kw)
    plt.close(fig)


def letter(ax, s, x=-0.1, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")


def read_shapes(category, name, attrs=False):
    path = shpreader.natural_earth(resolution="10m" if name == "admin_2_counties" else "50m",
                                   category=category, name=name)
    base = str(path)[:-4]
    r = shapefile.Reader(shp=open(base + ".shp", "rb"), shx=open(base + ".shx", "rb"),
                         dbf=open(base + ".dbf", "rb"), encoding="latin1")
    geoms = [sg.shape(s.__geo_interface__) for s in r.shapes()]
    if not attrs:
        return geoms
    fields = [f[0] for f in r.fields[1:]]
    return geoms, [dict(zip(fields, rec)) for rec in r.records()]


def haversine(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 6371 * 2 * np.arcsin(np.sqrt(a))


# ---- shared inputs -------------------------------------------------------------
STATES = read_shapes("cultural", "admin_1_states_provinces_lakes")
LAKES = read_shapes("physical", "lakes")
cg, crec = read_shapes("cultural", "admin_2_counties", attrs=True)
COUNTY = {rec["CODE_LOCAL"]: g for g, rec in zip(cg, crec) if rec["ISO_A2"] == "US"}

types = pd.read_csv(os.path.join(OUT, "shared_lines", "Table_SL_Locations.csv"), dtype={"mg": str})
coords = pd.read_csv(os.path.join(FD, "loc_coords.csv"))
LOC = types.merge(coords, left_on="location", right_on="loc_key")
assert len(LOC) == len(types), "every profiled location needs coordinates"
cent = pd.read_csv(os.path.join(ROOT, "data", "county_centroids.csv"), dtype={"fips": str})
cent["fips"] = cent.fips.str.zfill(5)
latc = [c for c in cent.columns if c.lower().startswith("lat")][0]
lonc = [c for c in cent.columns if c.lower().startswith(("lon", "lng"))][0]
cent[latc] = cent[latc].astype(float)
cent[lonc] = cent[lonc].astype(float)


def base_map(ax, extent, land="#f7f7f7"):
    ax.set_extent(extent, crs=PC)
    ax.add_feature(cf.ShapelyFeature(STATES, PC, facecolor=land, edgecolor="none"), zorder=0)


def overlay(ax):
    ax.add_feature(cf.ShapelyFeature(STATES, PC, facecolor="none", edgecolor="#555555", linewidth=0.45), zorder=2)
    ax.add_feature(cf.ShapelyFeature(LAKES, PC, facecolor="#dfe8ef", edgecolor="none"), zorder=2)
    ax.spines["geo"].set_linewidth(0.5)


def points(ax, x, s=30):
    for t, g in x.groupby("profile_type"):
        ax.scatter(g.lon, g.lat, s=s, color=PAL[int(t)], edgecolor="black", linewidth=0.7, transform=PC, zorder=4)


def panel_label(ax, k, mg):
    ax.text(0.015, 0.98, "abcde"[k], transform=ax.transAxes, fontweight="bold", fontsize=10, va="top", zorder=6)
    ax.text(0.06, 0.978, f"MG {mg}", transform=ax.transAxes, fontsize=9, va="top", zorder=6)


def type_legend(ax, y, title="Environment type (numbered within each maturity group)"):
    h = [plt.Line2D([], [], marker="s", ls="", color=PAL[t], markersize=8, label=f"Type {t}") for t in [1, 2, 3, 4]]
    leg = ax.legend(handles=h, loc="upper left", bbox_to_anchor=(0.06, y), frameon=False, fontsize=8, ncol=2,
                    title=title, title_fontsize=8, alignment="left")
    ax.add_artist(leg)


def map_grid(figh):
    fig = plt.figure(figsize=(6.85, figh))
    gs = fig.add_gridspec(3, 2, left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.03, hspace=0.04)
    return fig, gs


# ---- Fig. 1: environment types at the testing locations on county climate ------
gdd = pd.read_csv(os.path.join(FD, "county_gdd.csv"), dtype={"fips": str}).set_index("fips").GDD_mean
gdd = gdd[gdd.index.isin(cent.fips)]
norm = mcolors.Normalize(*np.nanpercentile(gdd.values, [2, 98]))
gcmap = mcolors.LinearSegmentedColormap.from_list("gdd", ["#fbfbfb", "#7a7a7a"])
fig, gs = map_grid(5.6)
for k, mg in enumerate(MGS):
    ax = fig.add_subplot(gs[k // 2, k % 2], projection=PROJ)
    x = LOC[LOC.mg == mg]
    ext = [x.lon.min() - 2.2, x.lon.max() + 2.2, x.lat.min() - 1.4, x.lat.max() + 1.6]
    base_map(ax, ext, land="#eeeeee")
    inbox = cent[(cent[lonc].between(ext[0] - 2, ext[1] + 2)) & (cent[latc].between(ext[2] - 2, ext[3] + 2))]
    for key, f in inbox.assign(c=np.round(norm(gdd.reindex(inbox.fips).values), 2)).groupby("c"):
        col = gcmap(float(key))
        ax.add_geometries([COUNTY[q] for q in f.fips if q in COUNTY], PC, facecolor=col, edgecolor=col, linewidth=0.05, zorder=1)
    overlay(ax)
    points(ax, x)
    panel_label(ax, k, mg)
    if mg == "II":
        for nm, off in [("hooper", (-36, 12)), ("mead", (-30, -16))]:
            r = x[x.location == nm].iloc[0]
            ax.annotate(nm.title(), xy=(r.lon, r.lat), xycoords=PC._as_mpl_transform(ax), xytext=off,
                        textcoords="offset points", fontsize=7, zorder=7,
                        arrowprops=dict(arrowstyle="-", lw=0.5, color="#333333"))
ax = fig.add_subplot(gs[2, 1]); ax.axis("off")
type_legend(ax, 0.98)
cax = ax.inset_axes([0.08, 0.5, 0.62, 0.06])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=gcmap), cax=cax, orientation="horizontal")
cb.set_label("County growing degree days, May to September\n(mean of 1993 to 2024)", fontsize=8)
cb.ax.tick_params(labelsize=7)
h3 = [plt.Line2D([], [], marker="o", ls="", markerfacecolor="#bbbbbb", markeredgecolor="black",
                 markeredgewidth=0.7, markersize=6, label="Testing location, colored by its type\nfrom marker profiles")]
ax.legend(handles=h3, loc="upper left", bbox_to_anchor=(0.06, 0.16), frameon=False, fontsize=8)
save(fig, "Fig1")

# ---- Fig. 2: profile distance and response distance ------------------------------
pairs = pd.read_csv(os.path.join(FD, "pairs.csv"), dtype={"mg": str})
t1 = pd.read_csv(os.path.join(OUT, "Table_PhenoMantel_Final.csv"), dtype={"mg": str}).set_index("mg")
dj = pd.read_csv(os.path.join(OUT, "Table_DisjointMantel.csv"), dtype={"mg": str}).set_index("mg")
p80 = pd.read_csv(os.path.join(OUT, "power", "Table_PartialMantel_Power80.csv"), dtype={"mg": str}).set_index("mg")
for m in MG4:  # pairs reproduce Table 1
    x = pairs[pairs.mg == m]
    assert abs(np.corrcoef(x.gen, x.pheno)[0, 1] - t1.loc[m, "mantel_genomic"]) < 1e-6
fig = plt.figure(figsize=(6.85, 4.4))
gs = fig.add_gridspec(2, 4, height_ratios=[1, 0.95], left=0.08, right=0.99, top=0.91, bottom=0.08, wspace=0.12, hspace=0.75)
ylim = (pairs.pheno.min() - 0.05, pairs.pheno.max() + 0.05)
for k, m in enumerate(MG4):
    ax = fig.add_subplot(gs[0, k]); x = pairs[pairs.mg == m]
    ax.scatter(x.gen, x.pheno, s=5, color=DARK, alpha=0.6, lw=0)
    b = np.polyfit(x.gen, x.pheno, 1); xx = np.linspace(x.gen.min(), x.gen.max(), 20)
    ax.plot(xx, np.polyval(b, xx), color=GEN, lw=1)
    ax.set_title(f"MG {m} ({len(x)} location pairs)", loc="left", pad=2)
    ax.text(0.04, 0.96, f"r = {t1.loc[m, 'mantel_genomic']:.2f}", transform=ax.transAxes, va="top", fontsize=7)
    ax.set_xlabel("Genomic profile distance", fontsize=7); ax.set_ylim(*ylim)
    if k == 0:
        ax.set_ylabel("Phenotypic response distance", fontsize=7)
    else:
        ax.set_yticklabels([])
    letter(ax, "abcd"[k], x=-0.38 if k == 0 else -0.08, y=1.1)
ax = fig.add_subplot(gs[1, :])
rows = [("Genomic profile distance", "mantel_genomic", GEN, "o", GEN),
        ("Envirotyping distance", "mantel_enviro", ENV, "o", ENV),
        ("Geographic distance", "mantel_geo", GEO, "o", GEO),
        ("Envirotyping distance, given genomic profile distance", "partial_env_gen", ENV, "D", "white")]
xs = np.arange(4); w = 0.18
for j, (lab, col, c, mk, fc) in enumerate(rows):
    xp = xs + (j - 1.5) * w
    ax.scatter(xp, [t1.loc[m, col] for m in MG4], s=26, marker=mk, facecolor=fc, edgecolor=c, lw=1.2, zorder=3, label=lab)
    if j == 0:
        for q, m in zip(xp, MG4):
            ax.plot([q, q], [dj.loc[m, "q05"], dj.loc[m, "q95"]], color=c, lw=1.6, alpha=0.45, zorder=2, solid_capstyle="butt")
    if j == 3:
        for q, m in zip(xp, MG4):
            ax.plot([q - 0.07, q + 0.07], [p80.loc[m, "partial_r_first_80"]] * 2, color="#444444", lw=1.0, zorder=2)
ax.axhline(0, color="#999999", lw=0.5)
ax.set_xticks(xs); ax.set_xticklabels([f"MG {m}" for m in MG4]); ax.set_xlim(-0.55, 3.55); ax.set_ylim(-0.1, 0.72)
ax.set_ylabel("Mantel correlation with\nphenotypic response distance", fontsize=7)
h, l = ax.get_legend_handles_labels()
h += [plt.Line2D([], [], color=GEN, alpha=0.45, lw=1.6), plt.Line2D([], [], color="#444444", lw=1.0)]
l += ["Genomic, profiles and responses from separate halves of the lines\n(5th to 95th percentile of 30 splits)",
      "Partial correlation detectable with 80% power"]
ax.legend(h, l, loc="upper right", frameon=False, fontsize=6.5, ncol=2, bbox_to_anchor=(1.0, 1.36), columnspacing=1.2, handletextpad=0.4)
letter(ax, "e", x=-0.075, y=1.2)
save(fig, "Fig2")

# ---- Fig. 3: groupings from shared lines ------------------------------------------
SL = os.path.join(OUT, "shared_lines")
ari = pd.read_csv(os.path.join(SL, "Table_SL_ARI.csv"), dtype={"mg": str})
summ = pd.read_csv(os.path.join(SL, "Table_SL_Summary.csv"), dtype={"mg": str})
pc = pd.read_csv(os.path.join(SL, "Table_SL_PairCounts.csv"), dtype={"mg": str})
un = pd.read_csv(os.path.join(SL, "Table_SL_UnplacedLocations.csv"), dtype={"mg": str})
FOC, GREY, ORANGE = "#1f4e79", "#9a9a9a", "#d08c2c"
fig, axs = plt.subplots(2, 2, figsize=(6.85, 5.4))
ax = axs[0, 0]; x = np.arange(4)
v = ari[(ari.location_set == "validation") & (ari.profile_grouping == "published profile types")]
a_k = v[v.response_grouping == "response distance, k = profile k"].set_index("mg").loc[MG4]
a_s = v[v.response_grouping == "response distance, silhouette k"].set_index("mg").loc[MG4]
for i, m in enumerate(MG4):
    for dx, q in ((-0.08, a_k.chance_q95[m]), (0.08, a_s.chance_q95[m])):
        ax.add_patch(mpatches.Rectangle((i + dx - 0.07, -0.2), 0.14, q + 0.2, color="#e6e6e6", lw=0, zorder=0))
ax.axhline(0, color=GREY, lw=0.6, zorder=1)
ax.scatter(x - 0.08, a_k.ari, color=FOC, s=28, zorder=3, label="Same number of groups as profile types")
ax.scatter(x + 0.08, a_s.ari, facecolor="white", edgecolor=FOC, s=28, lw=1.1, zorder=3, label="Number of groups chosen by silhouette")
ax.add_patch(mpatches.Rectangle((0, 0), 0, 0, color="#e6e6e6", label="Up to the 95th percentile of chance agreement"))
ax.set_xticks(x, [f"MG {m}" for m in MG4]); ax.set_ylabel("Adjusted Rand index")
ax.set_ylim(-0.2, 0.62); ax.set_xlim(-0.5, 3.5)
ax.set_title("Profile types and groups formed on response distance", loc="left")
ax.legend(frameon=False, loc="upper left", fontsize=6.5, handletextpad=0.3)
ax = axs[0, 1]
s_ = summ[summ.location_set == "validation"].set_index("mg").loc[MG4]
for i, m in enumerate(MG4):
    ax.plot([i, i], [s_.response_dist_within_types[m], s_.response_dist_between_types[m]], color=GREY, lw=1, zorder=1)
ax.scatter(x, s_.response_dist_within_types, color=FOC, s=28, zorder=3, label="Pairs of the same profile type")
ax.scatter(x, s_.response_dist_between_types, color=ORANGE, s=28, zorder=3, marker="s", label="Pairs of different profile types")
ax.set_xticks(x, [f"MG {m}" for m in MG4]); ax.set_ylabel("Mean phenotypic response distance"); ax.set_ylim(0.5, 1.0)
ax.set_title("Response distance within and between profile types", loc="left")
ax.legend(frameon=False, loc="upper left", fontsize=6.5, handletextpad=0.3)
ax = axs[1, 0]
p = pc.set_index("mg").loc[MGS]
z = p.pairs_zero_shared; lo = p.pairs_lt8_shared - p.pairs_zero_shared; ok = p.location_pairs - p.pairs_lt8_shared
y = np.arange(5)[::-1]
ax.barh(y, ok, color="#d9d9d9", label="8 or more shared lines")
ax.barh(y, lo, left=ok, color="#6baed6", label="1 to 7 shared lines")
ax.barh(y, z, left=ok + lo, color="#08306b", label="No shared lines")
for yy, m in zip(y, MGS):
    ax.text(p.location_pairs[m] + 4, yy, f"{p.pairs_lt8_shared[m]} of {p.location_pairs[m]}", va="center", fontsize=6.5)
ax.set_yticks(y, [f"MG {m}" for m in MGS]); ax.set_xlabel("Pairs of profiled locations"); ax.set_xlim(0, 480)
ax.set_title("Location pairs by number of shared lines", loc="left")
ax.legend(frameon=False, loc="lower right", fontsize=6.5, handletextpad=0.3)
ax = axs[1, 1]
g = un.groupby("mg").agg(n=("agrees", "size"), agree=("agrees", "sum"), chance=("chance_rate_mean", "first")).loc[MG4]
ax.bar(x, g.agree, color=FOC, width=0.55, label="Same type as the profile")
ax.bar(x, g.n - g.agree, bottom=g.agree, color="#c6dbef", width=0.55, label="Different type")
ax.scatter(x, g.n * g.chance, marker="_", s=260, color=ORANGE, lw=2, zorder=4, label="Expected by chance")
ax.set_xticks(x, [f"MG {m}" for m in MG4]); ax.set_ylabel("Locations placed through shared lines"); ax.set_ylim(0, 11.5)
ax.set_title("Locations with too few shared lines for direct grouping", loc="left")
ax.legend(frameon=False, loc="upper left", fontsize=6.5, handletextpad=0.3)
for a_, l_ in zip(axs.flat, "abcd"):
    a_.text(-0.14, 1.02, l_, transform=a_.transAxes, fontsize=10, fontweight="bold", va="bottom")
fig.tight_layout()
save(fig, "Fig3")

# ---- Fig. 5: projected environment types by county --------------------------------
RADIUS = 200


def alpha_of(c):
    return 0.12 + 0.83 * np.clip((c - 0.4) / 0.6, 0, 1)


fig, gs = map_grid(5.6)
for k, mg in enumerate(MGS):
    ax = fig.add_subplot(gs[k // 2, k % 2], projection=PROJ)
    pr = pd.read_csv(os.path.join(OUT, f"County_Predictions_v2_all_87_MG{mg}.csv"), dtype={"fips": str, "predicted": str})
    pr["fips"] = pr.fips.str.zfill(5)
    pr = pr.merge(cent[["fips", latc, lonc]], on="fips")
    x = LOC[LOC.mg == mg]
    dist = np.min(haversine(pr[latc].values[:, None], pr[lonc].values[:, None], x.lat.values[None, :], x.lon.values[None, :]), axis=1)
    pr = pr[dist <= RADIUS]
    lo_ = np.r_[pr[lonc], x.lon]; la_ = np.r_[pr[latc], x.lat]
    base_map(ax, [lo_.min() - 0.6, lo_.max() + 0.6, la_.min() - 0.5, la_.max() + 0.9])
    for (t, ab), f in pr.assign(ab=np.round(alpha_of(pr.confidence) * 10) / 10).groupby(["predicted", "ab"]):
        col = mcolors.to_rgba(PAL[int(t)], alpha=ab)
        ax.add_geometries([COUNTY[q] for q in f.fips if q in COUNTY], PC, facecolor=col, edgecolor=col, linewidth=0.05, zorder=1)
    overlay(ax)
    points(ax, x)
    panel_label(ax, k, mg)
ax = fig.add_subplot(gs[2, 1]); ax.axis("off")
type_legend(ax, 0.98)
h2 = [mpatches.Patch(facecolor=mcolors.to_rgba(PAL[1], alpha_of(c)), edgecolor="none", label=f"{c:.1f}") for c in [0.4, 0.6, 0.8, 1.0]]
l2 = ax.legend(handles=h2, loc="upper left", bbox_to_anchor=(0.06, 0.62), frameon=False, fontsize=8, ncol=4,
               title="Predicted probability of the county's type", title_fontsize=8, alignment="left", handlelength=1.6, columnspacing=1.0)
ax.add_artist(l2)
h3 = [plt.Line2D([], [], marker="o", ls="", markerfacecolor="#bbbbbb", markeredgecolor="black", markeredgewidth=0.7,
                 markersize=6, label="Testing location, colored by its type\nfrom marker profiles")]
ax.legend(handles=h3, loc="upper left", bbox_to_anchor=(0.06, 0.32), frameon=False, fontsize=8)
save(fig, "Fig5")

# ---- Fig. 4: yield differences between types --------------------------------------
yn = pd.read_csv(os.path.join(OUT, "Table_YieldDivergence_Null.csv"), dtype={"mg": str}).set_index("mg")
yp = pd.read_csv(os.path.join(OUT, "Table_YieldPenalty_byYear.csv"), dtype={"mg": str})
ypen = pd.read_csv(os.path.join(OUT, "Table_YieldPenalty.csv"), dtype={"mg": str})
ypen = ypen[ypen.variant == "pen_primary"].set_index("mg")
fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.85, 2.7), gridspec_kw=dict(width_ratios=[1, 1.25]))
for i, m in enumerate(MG4):
    r = yn.loc[m]
    a1.plot([i, i], [r.null_mean_kgha, r.null_p95_kgha], color=GREY, lw=5, solid_capstyle="butt", alpha=0.6)
    a1.scatter(i, r.null_mean_kgha, marker="_", s=120, color="#6a6a6a", lw=1.5)
    a1.scatter(i, r.observed_kgha, s=34, color=GEN, zorder=3)
    a1.text(i + 0.12, r.observed_kgha, f"P = {r.p:.3f}", fontsize=6.5, va="center")
a1.set_xticks(range(4)); a1.set_xticklabels([f"MG {m}" for m in MG4]); a1.set_xlim(-0.5, 3.7); a1.set_ylim(150, 400)
a1.set_ylabel(f"Mean difference in a line's adjusted\nyield between types ({KGHA})", fontsize=7)
a1.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=GEN, label="Environment types"),
                   plt.Line2D([], [], color=GREY, lw=5, alpha=0.6, label="Random grouping of locations\n(mean to 95th percentile)")],
          frameon=False, fontsize=6.5, loc="lower right")
pp = yp[yp.mg.isin(MG4)].copy(); pp["kg"] = pp.penalty_bu * 67.25
for i, m in enumerate(MG4):
    v_ = pp[pp.mg == m].kg.values
    a2.scatter(i + np.random.default_rng(i).uniform(-0.18, 0.18, len(v_)), v_, s=6, color=DARK, alpha=0.55, lw=0)
    a2.plot([i - 0.28, i + 0.28], [ypen.loc[m, "mean_penalty_kgha"]] * 2, color=GEN, lw=1.6)
    a2.text(i, -120, f"{int(ypen.loc[m, 'n_years'])} years", ha="center", fontsize=6.5, color="#444444")
mx = pp.loc[pp.kg.idxmax()]
a2.annotate(f"MG {mx.mg}, {int(mx.Year)}", xy=(MG4.index(mx.mg), mx.kg), xytext=(MG4.index(mx.mg) + 0.45, mx.kg - 80), fontsize=6.5,
            arrowprops=dict(arrowstyle="-", lw=0.5))
a2.set_xticks(range(4)); a2.set_xticklabels([f"MG {m}" for m in MG4]); a2.set_ylim(-170, 1700)
a2.set_ylabel(f"Yield lost by choosing one best\nline for all types ({KGHA})", fontsize=7)
a2.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=DARK, alpha=0.55, markersize=3, label="One type in one year"),
                   plt.Line2D([], [], color=GEN, lw=1.6, label="Mean of the yearly averages")], frameon=False, fontsize=6.5, loc="upper left")
letter(a1, "a", x=-0.22, y=1.0); letter(a2, "b", x=-0.2, y=1.0)
fig.subplots_adjust(left=0.1, right=0.99, top=0.93, bottom=0.12, wspace=0.32)
save(fig, "Fig4")

# ---- Fig. S1: permutation test of cluster separation ------------------------------
dr = pd.read_csv(os.path.join(OUT, "permutation_null", "PermNull_classification_draws.csv"), dtype={"mg": str})
pn = pd.read_csv(os.path.join(OUT, "permutation_null", "Table_PermNull_classification.csv"), dtype={"mg": str}).set_index("mg")
val = [c for c in dr.columns if c not in ("mg", "perm", "mode", "k", "null_k")][-1]
fig, ax = plt.subplots(figsize=(3.3, 2.6))
parts = ax.violinplot([dr[dr.mg == m][val].values for m in MGS], positions=range(5), widths=0.7, showextrema=False)
for b in parts["bodies"]:
    b.set_facecolor("#cccccc"); b.set_edgecolor("none"); b.set_alpha(1)
ax.scatter(range(5), [pn.loc[m, "observed"] for m in MGS], color="black", s=18, zorder=3)
ax.set_xticks(range(5)); ax.set_xticklabels([f"MG {m}" for m in MGS]); ax.set_ylim(0, 0.42)
ax.set_ylabel("Mean silhouette width")
ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", color="black", markersize=4, label="Observed"),
                   mpatches.Patch(color="#cccccc", label="1,000 permutations")], frameon=False, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.12))
fig.tight_layout(); save(fig, "FigS1")

# ---- Fig. S3: covariate importance for the county model ---------------------------
NAMES = {
    "GDD_mean": "Growing degree days", "GDD_sd": "Growing degree days, SD between years",
    "T2M": "Air temperature (monthly mean)", "T2M_MAX": "Maximum air temperature", "T2M_MIN": "Minimum air temperature",
    "T2MDEW": "Dew point temperature", "T2MDEW_sd": "Dew point temperature, SD", "Tmean_mean": "Mean temperature",
    "Tmean_early_mean": "Mean temperature, May to June", "Tmean_late_mean": "Mean temperature, July to August",
    "Tmax_mean": "Mean daily maximum temperature", "Tmax_july_90th": "July maximum temperature, 90th percentile",
    "VPD_mean": "Vapor pressure deficit", "VPD_sd": "Vapor pressure deficit, SD", "VPD_90th": "Vapor pressure deficit, 90th percentile",
    "VPD_early_mean": "Vapor pressure deficit, May to June", "VPD_late_mean": "Vapor pressure deficit, July to August",
    "VPD_late_early_ratio": "Vapor pressure deficit, late to early ratio", "RH2M": "Relative humidity",
    "WS2M": "Wind speed", "WS2M_sd": "Wind speed, SD", "PS": "Surface pressure", "PS_sd": "Surface pressure, SD",
    "PRECTOTCORR_sd": "Precipitation, SD", "Precip_late_mean": "Precipitation, July to August",
    "Solar_early_mean": "Solar radiation, May to June", "Solar_late_mean": "Solar radiation, July to August",
    "awc_0_30": "Available water capacity, 0 to 30 cm", "awc_30_60": "Available water capacity, 30 to 60 cm",
    "awc_total": "Available water capacity, 0 to 100 cm", "awc_deep_frac": "Share of water capacity below 30 cm",
    "clay_0_30": "Clay, 0 to 30 cm", "ksat_0_30": "Saturated hydraulic conductivity, 0 to 30 cm",
    "om_30_60": "Organic matter, 30 to 60 cm", "om_60_100": "Organic matter, 60 to 100 cm", "cec_0_30": "Cation exchange capacity, 0 to 30 cm",
    "silt_60_100": "Silt, 60 to 100 cm", "elevation": "Elevation", "corn_soy_ratio": "Corn to soybean area ratio",
    "irrig_fraction": "Irrigated share of cropland", "tile_drain_fraction": "Tile-drained share of cropland",
    "avg_farm_size": "Average farm size", "Solar_mean": "Solar radiation",
}
fig, axs = plt.subplots(3, 2, figsize=(6.85, 7.0))
for k, (ax, mg) in enumerate(zip(axs.flat, MGS)):
    im = pd.read_csv(os.path.join(OUT, f"RF_Importance_v2_MG{mg}.csv")).sort_values("MDA", ascending=False).head(10)[::-1]
    ax.barh(range(len(im)), im.MDA, color=GEN, height=0.65)
    ax.set_yticks(range(len(im))); ax.set_yticklabels([NAMES.get(f, f) for f in im.feature], fontsize=6)
    ax.set_title(f"MG {mg}", loc="left", pad=2); ax.tick_params(axis="x", labelsize=6)
    ax.text(-0.02, 1.06, "abcde"[k], transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right")
axs.flat[5].axis("off")
fig.supxlabel("Mean decrease in accuracy when the covariate is permuted", fontsize=7)
fig.subplots_adjust(left=0.3, right=0.98, top=0.96, bottom=0.07, wspace=1.05, hspace=0.38); save(fig, "FigS3")

# ---- Fig. S2: two-dimensional scaling of profile distance --------------------------
mds = pd.read_csv(os.path.join(FD, "mds.csv"), dtype={"mg": str})
fig, axs = plt.subplots(2, 3, figsize=(6.85, 4.4))
for k, (ax, mg) in enumerate(zip(axs.flat, MGS)):
    x = mds[mds.mg == mg]
    for t, g in x.groupby("type"):
        ax.scatter(g.MDS1, g.MDS2, s=14, color=PAL[int(t)], edgecolor="white", lw=0.3)
    ax.set_title(f"MG {mg} ({len(x)} locations, {int(x.k.iloc[0])} types)", loc="left", pad=2)
    ax.set_xticks([]); ax.set_yticks([]); ax.margins(0.08)
    ax.text(-0.02, 1.06, "abcde"[k], transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right")
    if k in (0, 3):
        ax.set_ylabel("Scaling axis 2")
    if k >= 3:
        ax.set_xlabel("Scaling axis 1")
ax = axs.flat[5]; ax.axis("off")
type_legend(ax, 0.8, title="Environment type")
fig.tight_layout(); save(fig, "FigS2")
print("figures written to", FIGDIR)
