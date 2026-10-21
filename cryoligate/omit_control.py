"""
Omit control for CryoLigATE predictions.

CryoLigATE is trained against density simulated from atomic models and can learn priors about what
ligand density looks like. The omit control estimates how much of the predicted ligand density is
supported by the experimental data: the density around the putative ligand is replaced by noise
matching the solvent of the same volume, the network is run again (E_omit), and the result is compared
with the standard prediction (E). Density that appears in both runs reflects the network's priors;
E - E_omit shows the density attributable to the experimental data.

The ligand Q-score (Pintilie et al., Nat. Methods 2020) is computed here in NumPy, at the deposited
ligand position and without re-fitting, following the ChimeraX "qscore" implementation (shells from
0.1 to 2.0 Å in 0.1 Å steps, 8 points per shell, reference Gaussian with sigma 0.6 Å). Only the ligand's
own atoms are used to assign shell points, as in the evaluation reported in the manuscript.
"""
import numpy as np
import scipy.ndimage

OMIT_RADIUS = 2.5          # Å around ligand atoms replaced by noise
BACKGROUND_DISTANCE = 4.0  # Å from any atom: voxels used to estimate the solvent background
SIGNAL_THRESHOLD = 0.1     # Q(E) - Q(E_omit) below this: density likely driven by learned priors
PRESENT_THRESHOLD = 0.2    # Q(E) below this: no ligand-like density predicted at the putative position

# Q-score parameters (ChimeraX qscore defaults)
REF_SIGMA = 0.6
MAX_RADIUS = 2.0
STEP = 0.1
POINTS_PER_SHELL = 8
TEST_POINTS = 128


def atom_occupancy(coords, phys_origin, grid_size, voxel):
    """Boolean grid with True at voxels containing an atom (as used to build the training masks)."""
    grid = np.zeros((grid_size,) * 3, dtype=bool)
    if len(coords):
        v = np.round((np.asarray(coords) - phys_origin) / voxel).astype(int)
        inside = np.all((v >= 0) & (v < grid_size), axis=1)
        v = v[inside]
        grid[v[:, 0], v[:, 1], v[:, 2]] = True
    return grid


def make_omit_input(density, ligand_coords, other_coords, phys_origin, voxel, seed=0):
    """Replace density within OMIT_RADIUS of the ligand atoms with noise matching the solvent background.

    `density` is indexed (x, y, z) and already normalised as the network input.
    """
    n = density.shape[0]
    ligand = atom_occupancy(ligand_coords, phys_origin, n, voxel)
    atoms = ligand | atom_occupancy(other_coords, phys_origin, n, voxel)
    near_ligand = scipy.ndimage.distance_transform_edt(~ligand) * voxel <= OMIT_RADIUS
    background = scipy.ndimage.distance_transform_edt(~atoms) * voxel > BACKGROUND_DISTANCE
    if background.sum() < 100:   # crowded box: fall back to everything outside the omitted region
        background = ~near_ligand
    bg = density[background]
    rng = np.random.default_rng(seed)
    omitted = density.copy()
    omitted[near_ligand] = rng.normal(bg.mean(), bg.std(), size=int(near_ligand.sum()))
    return omitted


def _fibonacci_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)


def _cluster_representatives(points, k, iterations=5):
    """k-means on the candidate shell points; return the point nearest each cluster centre (as ChimeraX)."""
    # deterministic, well-spread seeds (farthest-point sampling)
    seeds = [0]
    d = np.linalg.norm(points - points[0], axis=1)
    while len(seeds) < k:
        nxt = int(np.argmax(d))
        seeds.append(nxt)
        d = np.minimum(d, np.linalg.norm(points - points[nxt], axis=1))
    centres = points[seeds].copy()
    for _ in range(iterations):
        labels = np.argmin(np.linalg.norm(points[:, None, :] - centres[None], axis=2), axis=1)
        for c in range(k):
            if np.any(labels == c):
                centres[c] = points[labels == c].mean(axis=0)
    nearest = np.argmin(np.linalg.norm(points[:, None, :] - centres[None], axis=2), axis=0)
    return points[np.unique(nearest)]


def ligand_q_score(volume, phys_origin, voxel, ligand_coords):
    """Mean Q-score of the ligand atoms in `volume` (indexed x, y, z) at their given positions."""
    coords = np.asarray(ligand_coords, dtype=float)
    radii = np.arange(STEP, MAX_RADIUS + STEP / 2, STEP)
    sphere = _fibonacci_sphere(TEST_POINTS)
    sphere_small = _fibonacci_sphere(POINTS_PER_SHELL)
    upper = np.array(volume.shape) - 1

    def owned(pts, centre, others):
        """Points closer to this atom than to any other ligand atom."""
        if not len(others):
            return pts
        d_self = np.linalg.norm(pts - centre, axis=1)
        d_other = np.min(np.linalg.norm(pts[:, None, :] - others[None, :, :], axis=2), axis=1)
        return pts[d_self <= d_other]

    def sample(points):
        idx = (points - phys_origin) / voxel
        return scipy.ndimage.map_coordinates(volume, idx.T, order=1, mode="nearest")

    scores = []
    for a, centre in enumerate(coords):
        idx = (centre - phys_origin) / voxel
        if np.any(idx < 0) or np.any(idx > upper):
            continue                                    # atom outside the box, as ChimeraX excludes it
        others = np.delete(coords, a, axis=0)
        measured = [np.repeat(sample(centre[None]), POINTS_PER_SHELL)]
        reference = [np.ones(POINTS_PER_SHELL)]
        for r in radii:
            # small radii: try 8 evenly spread points first (all must belong to this atom)
            pts = centre + sphere_small * r
            if r >= 0.7 or len(owned(pts, centre, others)) < POINTS_PER_SHELL:
                pts = owned(centre + sphere * r, centre, others)
                if len(pts) == 0:
                    continue
                if len(pts) > POINTS_PER_SHELL:
                    pts = _cluster_representatives(pts, POINTS_PER_SHELL)
            measured.append(sample(pts))
            reference.append(np.full(len(pts), np.exp(-0.5 * (r / REF_SIGMA) ** 2)))
        u, v = np.concatenate(measured), np.concatenate(reference)
        u, v = u - u.mean(), v - v.mean()
        denom = np.sqrt((u * u).sum() * (v * v).sum())
        scores.append(float((u * v).sum() / denom) if denom > 0 else 0.0)
    return float(np.mean(scores)) if scores else float("nan")


def report(identifier, q_original, q_pred, q_omit):
    signal = q_pred - q_omit
    no_density = q_pred < PRESENT_THRESHOLD
    supported = (not no_density) and signal >= SIGNAL_THRESHOLD
    lines = [
        f"CryoLigATE omit control: {identifier}",
        "Ligand Q-score at the input ligand position (no re-fitting):",
        f"  original map             {q_original:6.2f}",
        f"  CryoLigATE (E)           {q_pred:6.2f}",
        f"  CryoLigATE omit (E_omit) {q_omit:6.2f}",
        f"  signal-supported Q(E) - Q(E_omit) = {signal:+.2f}",
        "",
        (f"Interpretation: no ligand-like density was predicted at this position (Q(E) < {PRESENT_THRESHOLD}). "
         "CryoLigATE does not support the putative ligand here."
         if no_density else
         "Interpretation: the predicted ligand density depends on the experimental data at this site. "
         "Inspect difference_*.mrc (E - E_omit) for the supported density."
         if supported else
         f"Interpretation: CAUTION. Q(E) - Q(E_omit) is below {SIGNAL_THRESHOLD}. The predicted ligand density at "
         "this site is likely driven by the network's learned priors and should not be taken as evidence "
         "for the ligand. Rely on the original map."),
        "",
        f"For reference, on the CryoLigATE test set (n = 647) the median Q(E) - Q(E_omit) was 0.34, "
        f"and 14% of complexes were below {SIGNAL_THRESHOLD}.",
    ]
    return "\n".join(lines), signal, supported
