"""utils.py — MAP-Elites JSON loading and phenotype utilities.

Pure Python/NumPy — no Congrego package dependency.
Used by explore_mapelites.ipynb.
"""

import json
import os
import re

import numpy as np
import pandas as pd


# ── JSON loading ───────────────────────────────────────────────────────────────

def _parse_descriptor_key(s):
    """Parse '(10, 1250)' → (10, 1250)."""
    m = re.match(r'\((\d+),\s*(\d+)\)', s)
    return (int(m.group(1)), int(m.group(2)))


def load_mapelites_json(filepath):
    """
    Returns (archive, pc_meta, max_fitness, metadata).

    - archive: dict[(n_slots, price)] = dict with keys:
        fitness (float), phenotypes (list[int]), n_circuits (int), generator (raw slots)
    - pc_meta: raw part-catalogue dict (string-keyed from JSON)
        key fields: "inputs" (list[str]), "parts", "activatable_promoters", etc.
    - max_fitness: list[float]
    - metadata: raw dict with "algorithm" and "run" sub-dicts
    """
    with open(filepath) as f:
        d = json.load(f)
    archive = {_parse_descriptor_key(k): v for k, v in d["archive"].items()}
    # TODO: remove "max_fitness_per_generation" fallback once all old archive JSONs are retired
    max_fitness = d.get("max_fitness_per_iteration") or d["max_fitness_per_generation"]
    return archive, d["part_catalogue"], max_fitness, d["metadata"]


# ── Phenotype utilities ────────────────────────────────────────────────────────

def n_phenotypes(k):
    """2**(2**k)"""
    return 2 ** (2 ** k)


def pheno_labels(k):
    """
    Mirror of Julia pheno_labels(k):
    k=1 → ['zero', 'identity', 'not', 'one']                          (hardcoded)
    k=2 → ['false','and','a_not_b','a','b_not_a','b','xor','or',
            'nor','xnor','not_b','a_or_not_b','not_a','not_a_or_b','nand','true']
    k≥3 → binary strings '00000000'…'11111111' (2^(2^k) entries)
    """
    if k == 1:
        return ['zero', 'identity', 'not', 'one']
    if k == 2:
        return [
            'false', 'and', 'a_not_b', 'a', 'b_not_a', 'b', 'xor', 'or',
            'nor', 'xnor', 'not_b', 'a_or_not_b', 'not_a', 'not_a_or_b', 'nand', 'true',
        ]
    n = n_phenotypes(k)
    width = 2 ** k
    return [format(i, f'0{width}b') for i in range(n)]


# ── Archive helpers ────────────────────────────────────────────────────────────

def archive_matrix(archive, f):
    """
    Build 2-D matrix from archive for heatmaps.

    Returns (ns, ps, mat) where:
    - ns  = sorted unique n_slots values
    - ps  = sorted unique price values
    - mat = np.ndarray shape (len(ns), len(ps)), NaN for missing cells
    """
    ns = sorted(set(k[0] for k in archive))
    ps = sorted(set(k[1] for k in archive))
    ni = {v: i for i, v in enumerate(ns)}
    pi = {v: i for i, v in enumerate(ps)}
    mat = np.full((len(ns), len(ps)), np.nan)
    for k, cell in archive.items():
        mat[ni[k[0]], pi[k[1]]] = f(cell)
    return ns, ps, mat


def pheno_entropy(cell):
    """Shannon entropy of normalised phenotype distribution."""
    p = np.array(cell['phenotypes'], dtype=float) / cell['n_circuits']
    return -float(np.sum([x * np.log(x) for x in p if x > 0]))


def best_cell(archive):
    """Return (key, cell) for the cell with highest fitness."""
    return max(archive.items(), key=lambda kv: kv[1]['fitness'])


# ── Useful extras ──────────────────────────────────────────────────────────────

def summarize_cell(key, cell, pc_meta):
    """
    Return a multi-line string summary of an archive cell:
    n_slots, price, fitness, n_circuits, coverage (% of phenotypes present).
    """
    k = len(pc_meta['inputs'])
    n_ph = n_phenotypes(k)
    coverage = sum(1 for x in cell['phenotypes'] if x > 0)
    lines = [
        f"n_slots  : {key[0]}",
        f"price    : {key[1]} GBP",
        f"fitness  : {cell['fitness']:.4f}",
        f"circuits : {cell['n_circuits']}",
        f"coverage : {coverage}/{n_ph} ({100 * coverage / n_ph:.1f}%)",
    ]
    return '\n'.join(lines)


def reachability(archive, k):
    """
    Return array of length n_phenotypes(k): for each function, count of archive
    cells that have ≥1 circuit achieving it.
    """
    n_ph = n_phenotypes(k)
    return np.array([
        sum(1 for cell in archive.values() if cell['phenotypes'][i] > 0)
        for i in range(n_ph)
    ])


def presence_matrix(archive, k):
    """
    Return (sorted_keys, mat) where mat is shape (n_phenotypes(k), n_cells),
    sorted by fitness ascending. Binary float: 1.0 if phenotype present, 0.0 otherwise.
    """
    sorted_items = sorted(archive.items(), key=lambda kv: kv[1]['fitness'])
    sorted_keys = [kv[0] for kv in sorted_items]
    n_ph = n_phenotypes(k)
    mat = np.array([
        [1.0 if cell['phenotypes'][i] > 0 else 0.0 for _, cell in sorted_items]
        for i in range(n_ph)
    ])
    return sorted_keys, mat


# ── Sweep loader (generic — caller supplies files + data_dir) ─────────────────

def load_sweep(files, data_dir):
    """Load every sweep JSON listed in `files` and build three structures.

    Generic across sweeps — caller supplies the manifest and directory so the
    same function works for any future sweep without modification.

    Parameters
    ----------
    files : iterable of tuples
        Each tuple: (filename, run_idx, metric, w_diversity, c_parsimony).
    data_dir : str
        Directory containing the JSON files.

    Returns
    -------
    df_runs : pd.DataFrame       one row per run (metadata + best-cell summary)
    df_archive : pd.DataFrame    one row per archive cell across all runs
    convergence : dict[int, np.ndarray]   run_idx → max-fitness-per-iteration
    """
    runs_rows = []
    archive_rows = []
    convergence = {}

    for filename, run_idx, metric, w_div, c_pars in files:
        filepath = os.path.join(data_dir, filename)
        with open(filepath) as f:
            data = json.load(f)

        archive = data["archive"]
        convergence[run_idx] = np.array(data["max_fitness_per_iteration"])

        best_fitness = -np.inf
        best_components = None
        for key_str, cell in archive.items():
            # Descriptor key looks like "(6, 1500)" — parse it
            n_slots, price = [int(x.strip("() ")) for x in key_str.split(",")]
            fc = cell["fitness_components"]
            phenotypes = cell["phenotypes"]
            n_phenotypes_covered = sum(1 for c in phenotypes if c > 0)

            archive_rows.append({
                "run_idx": run_idx,
                "metric": metric,
                "w_diversity": w_div,
                "c_parsimony": c_pars,
                "n_slots": n_slots,
                "price": price,
                "fitness": cell["fitness"],
                "entropy": fc["entropy"],
                "coverage": fc["coverage"],
                "diversity": fc["diversity"],
                "norm_mean_part": fc["norm_mean_part"],
                "n_circuits": cell["n_circuits"],
                "n_phenotypes": n_phenotypes_covered,
                "phenotypes": phenotypes,
            })

            if cell["fitness"] > best_fitness:
                best_fitness = cell["fitness"]
                best_components = fc

        runs_rows.append({
            "run_idx": run_idx,
            "metric": metric,
            "w_diversity": w_div,
            "c_parsimony": c_pars,
            "n_cells": len(archive),
            "best_fitness": best_fitness,
            "best_entropy": best_components["entropy"],
            "best_coverage": best_components["coverage"],
            "best_diversity": best_components["diversity"],
            "best_norm_mean_part": best_components["norm_mean_part"],
        })

    df_runs = pd.DataFrame(runs_rows)
    df_archive = pd.DataFrame(archive_rows)
    return df_runs, df_archive, convergence


# ── Archive heatmap primitive (generic — caller supplies cmap) ────────────────

def plot_archive_heatmap(ax, df, *, df_archive, cmap,
                         component="fitness", vmin=0.0, vmax=1.0,
                         title=None, show_xlabel=True, show_ylabel=True,
                         show_xticklabels=True, show_yticklabels=True):
    """One panel of a descriptor-space heatmap.

    Parameters
    ----------
    ax : matplotlib Axes
    df : pd.DataFrame or None
        Filtered cells for this panel (one run, or an aggregation).
        If None or empty, draws a 'No data' placeholder.
    df_archive : pd.DataFrame
        Full archive; used only to derive consistent (n_slots, price) axes so
        panels across a grid share the same descriptor ordering.
    cmap : matplotlib Colormap
        Caller-supplied colormap — styling is not the utility's concern.
    component : str
        Column in `df` to colour by.
    vmin, vmax : float
        Colour scale bounds.
    title : str, optional
    show_xlabel, show_ylabel : bool
        Whether to draw the axis *title* text (e.g. "Price (GBP)").
    show_xticklabels, show_yticklabels : bool
        Whether to draw tick *labels*. Tick marks are always drawn so
        adjacent panels remain visually aligned; only the labels are hidden.
        Set False on panels that share an axis with a neighbour to avoid
        collisions in tight grids (mirrors matplotlib's sharex/sharey).

    Returns
    -------
    AxesImage or None
        The image object (for attaching a colorbar). None if no data.
    """
    # Derive consistent descriptor axes from the full archive
    n_slots_vals = sorted(df_archive["n_slots"].unique())
    price_vals = sorted(df_archive["price"].unique())
    slot_idx = {v: i for i, v in enumerate(n_slots_vals)}
    price_idx = {v: i for i, v in enumerate(price_vals)}

    if df is None or df.empty:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center", fontsize=9, color="#888")
        ax.set_facecolor("#f5f5f5")
        if title:
            ax.set_title(title, fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        return None

    grid = np.full((len(n_slots_vals), len(price_vals)), np.nan)
    cov_x, cov_y = [], []
    for _, row in df.iterrows():
        si, pi = slot_idx[row["n_slots"]], price_idx[row["price"]]
        grid[si, pi] = row[component]
        if row["n_phenotypes"] == 16:
            cov_x.append(pi)
            cov_y.append(si)

    im = ax.imshow(grid, aspect="auto", origin="lower", cmap=cmap,
                   vmin=vmin, vmax=vmax)
    ax.scatter(cov_x, cov_y, marker="o", s=16, facecolors="white",
               edgecolors="#333", linewidths=0.4, alpha=0.9,
               zorder=3, label="16/16 coverage")

    # Tick marks always shown so shared-axis panels align visually; labels
    # are optional so tight grids can hide redundant text.
    ax.set_xticks(range(0, len(price_vals), 2))
    if show_xticklabels:
        ax.set_xticklabels([price_vals[i] for i in range(0, len(price_vals), 2)],
                           rotation=45, ha="right", fontsize=7)
    else:
        ax.set_xticklabels([])
    ax.set_yticks(range(len(n_slots_vals)))
    if show_yticklabels:
        ax.set_yticklabels(n_slots_vals, fontsize=8)
    else:
        ax.set_yticklabels([])

    if show_xlabel:
        ax.set_xlabel("Price (GBP)")
    if show_ylabel:
        ax.set_ylabel("Number of slots")
    if title:
        ax.set_title(title, fontsize=10)

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    return im


# ── Component-stack primitive (generic — one panel per filtered slice) ────────

def plot_components_stack(ax, df, *, colors,
                          component_cols=("norm_mean_part", "entropy",
                                          "coverage", "diversity"),
                          component_labels=None,
                          ylim=(0, 3.5),
                          mark_full_coverage=True,
                          title=None,
                          show_xlabel=True, show_ylabel=True,
                          show_xticklabels=True, show_yticklabels=True,
                          show_legend=False):
    """One panel: stacked-bar chart of fitness components vs price.

    Unlike `plot_archive_heatmap`, this primitive shows the *decomposition*
    of unweighted fitness — each bar is the sum of four components, and the
    stack order reveals which term dominates for each archive cell.

    Parameters
    ----------
    ax : matplotlib Axes
    df : pd.DataFrame or None
        Pre-filtered rows (one metric × one c_p × one w_d × one n_slots).
        Required columns: `price`, `n_phenotypes`, and every key in `colors`.
    colors : dict[str, str]
        Maps component name → matplotlib colour. Iteration order defines the
        stack order (bottom → top). Python 3.7+ dicts preserve insertion order.
    component_cols : tuple[str, ...]
        Which fitness-component columns of `df` to stack. Defaults to the
        four tracked across the sweep.
    component_labels : dict[str, str] or None
        Optional display-name overrides for the legend. Maps df column name
        → label text; columns not present fall back to the column name.
    ylim : tuple[float, float]
        Y-axis limits. 3.5 leaves headroom for the full-coverage annotation
        above the tallest bar.
    mark_full_coverage : bool
        Mark bars where n_phenotypes == 16 with a white dot (matching the
        heatmap coverage markers).
    show_{xlabel,ylabel,xticklabels,yticklabels} : bool
        Matches plot_archive_heatmap's shared-axis convention — outer panels
        draw labels, inner panels hide them.
    show_legend : bool
        Axis-level legend. Caller turns on for one panel only.

    Returns
    -------
    list[matplotlib.container.BarContainer]
        The bar containers in stack order, or [] if `df` was empty.
    """
    if df is None or df.empty:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center", fontsize=9, color="#888")
        ax.set_facecolor("#f5f5f5")
        if title:
            ax.set_title(title, fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        return []

    df_sorted = df.sort_values("price").reset_index(drop=True)
    prices = df_sorted["price"].to_numpy()
    x = np.arange(len(prices))

    _labels = component_labels or {}

    bottom = np.zeros(len(prices))
    bars = []
    for col in component_cols:
        vals = df_sorted[col].to_numpy(dtype=float)
        bc = ax.bar(x, vals, bottom=bottom, color=colors[col],
                    label=_labels.get(col, col), width=0.7, edgecolor="none")
        bars.append(bc)
        bottom = bottom + vals

    if mark_full_coverage:
        full = df_sorted["n_phenotypes"].to_numpy() == 16
        if np.any(full):
            ax.scatter(x[full], bottom[full] + 0.08, marker="o", s=16,
                       facecolors="white", edgecolors="#333",
                       linewidths=0.4, alpha=0.9, zorder=3,
                       label="16/16 coverage")

    ax.set_ylim(*ylim)
    ax.set_xticks(x)
    if show_xticklabels:
        ax.set_xticklabels([str(p) for p in prices],
                           rotation=60, ha="right", fontsize=7)
    else:
        ax.set_xticklabels([])
    if not show_yticklabels:
        ax.set_yticklabels([])

    if show_xlabel:
        ax.set_xlabel("Price (GBP)")
    if show_ylabel:
        ax.set_ylabel("Fitness")
    if title:
        ax.set_title(title, fontsize=10)
    if show_legend:
        ax.legend(loc="upper left", fontsize=6, framealpha=0.85)

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    return bars


# ── Convergence analysis (per-cell history + plateau detection) ───────────────

def load_archive_history(files, data_dir, n_slots=6):
    """Long-form per-cell improvement events for archive cells of a given n_slots.

    Reads each sweep JSON's ``archive_history`` — the per-cell log of
    ``{iter, fitness, entropy, coverage, diversity, norm_mean_part}`` written
    whenever a cell's elite is replaced — and returns the events for cells whose
    descriptor has ``n_slots`` slots (the niche of interest).

    Parameters
    ----------
    files : iterable of tuples
        Same manifest as ``load_sweep``:
        (filename, run_idx, metric, w_diversity, c_parsimony).
    data_dir : str
    n_slots : int
        Keep only cells with this slot count.

    Returns
    -------
    pd.DataFrame
        One row per (cell, improvement event):
        [run_idx, metric, w_diversity, c_parsimony, price, iter,
         entropy, coverage, diversity].
    """
    rows = []
    for filename, run_idx, metric, w_div, c_pars in files:
        with open(os.path.join(data_dir, filename)) as fh:
            data = json.load(fh)
        for key_str, events in data["archive_history"].items():
            ns, price = _parse_descriptor_key(key_str)
            if ns != n_slots:
                continue
            for e in events:
                rows.append({
                    "run_idx": run_idx, "metric": metric,
                    "w_diversity": w_div, "c_parsimony": c_pars,
                    "price": price, "iter": e["iter"],
                    "entropy": e["entropy"], "coverage": e["coverage"],
                    "diversity": e["diversity"],
                })
    return pd.DataFrame(rows)


def plateau_iter(iters, values, *, tol=0.02):
    """Earliest iteration after which ``values`` stays within a relative ``tol``
    band of its final value for the entire remainder of the series.

    Parameters
    ----------
    iters, values : 1-D array-like, same length, ``iters`` ascending.
    tol : float
        Relative half-width of the band around ``values[-1]``.

    Returns
    -------
    (t_star, converged) : (number, bool)
        ``t_star`` is the iteration from which the series never again leaves the
        band. ``converged`` is False when the series is still outside the band at
        the final point (it never settled); then ``t_star`` is the last iteration
        — a lower bound on the true convergence point.
    """
    iters = np.asarray(iters)
    values = np.asarray(values, dtype=float)
    final = values[-1]
    band = tol * abs(final) if final != 0 else tol
    outside = np.abs(values - final) > band
    if not outside.any():
        return iters[0], True
    last_outside = int(np.flatnonzero(outside)[-1])
    if last_outside >= len(values) - 1:
        return iters[-1], False
    return iters[last_outside + 1], True


def elite_value_at(events, component, at):
    """Per-cell forward-filled `component` value as of iteration `at`.

    Parameters
    ----------
    events : pd.DataFrame
        Rows from ``load_archive_history`` (typically one run). Must contain
        ``price``, ``iter`` and the named ``component`` column.
    component : str
    at : int
        Query iteration.

    Returns
    -------
    pd.Series indexed by ``price`` — each cell's elite ``component`` value at
    iteration ``at`` (the last improvement with ``iter <= at``; 0.0 if the cell had
    no improvement by then).
    """
    out = {}
    for price, cell in events.groupby("price"):
        prior = cell[cell["iter"] <= at]
        out[price] = (float(prior.sort_values("iter")[component].iloc[-1])
                      if len(prior) else 0.0)
    return pd.Series(out, name=component)
