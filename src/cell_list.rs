//! Cell list spatial indexing for O(N) pairwise neighbor search
//!
//! Partitions 3D space into uniform cubic cells. For a given cutoff
//! distance, only atoms in the same or adjacent cells need to be checked,
//! reducing pairwise search from O(N^2) to O(N) for uniform distributions.
//!
//! Author: Moroya Sakamoto

/// Configuration for the cell list.
#[derive(Debug, Clone, Copy)]
pub struct CellListConfig {
    /// Cutoff distance (Angstroms). Pairs beyond this are ignored.
    pub cutoff: f64,
    /// Minimum coordinate on each axis (auto-computed if None).
    pub origin: Option<[f64; 3]>,
}

impl Default for CellListConfig {
    fn default() -> Self {
        Self {
            cutoff: 8.0,
            origin: None,
        }
    }
}

/// Cell list for efficient spatial neighbor queries.
#[derive(Debug, Clone)]
pub struct CellList {
    /// Cell size (= cutoff distance).
    cell_size: f64,
    /// Number of cells in each dimension.
    dims: [usize; 3],
    /// Origin (minimum corner) of the grid.
    origin: [f64; 3],
    /// For each cell, list of atom indices.
    cells: Vec<Vec<usize>>,
    /// Cell assignment for each atom.
    assignments: Vec<usize>,
}

impl CellList {
    /// Build a cell list from atom positions with given cutoff.
    pub fn build(positions: &[[f64; 3]], config: &CellListConfig) -> Self {
        if positions.is_empty() {
            return Self {
                cell_size: config.cutoff,
                dims: [0, 0, 0],
                origin: [0.0; 3],
                cells: Vec::new(),
                assignments: Vec::new(),
            };
        }

        let cell_size = config.cutoff;

        // Compute bounding box
        let mut min = [f64::MAX; 3];
        let mut max = [f64::MIN; 3];
        for p in positions {
            for i in 0..3 {
                if p[i] < min[i] { min[i] = p[i]; }
                if p[i] > max[i] { max[i] = p[i]; }
            }
        }

        let origin = config.origin.unwrap_or([
            min[0] - cell_size,
            min[1] - cell_size,
            min[2] - cell_size,
        ]);

        // Grid dimensions (at least 1 cell per axis)
        let inv_cell = 1.0 / cell_size;
        let dims = [
            ((max[0] - origin[0]) * inv_cell).ceil() as usize + 2,
            ((max[1] - origin[1]) * inv_cell).ceil() as usize + 2,
            ((max[2] - origin[2]) * inv_cell).ceil() as usize + 2,
        ];

        let total_cells = dims[0] * dims[1] * dims[2];
        let mut cells = vec![Vec::new(); total_cells];
        let mut assignments = Vec::with_capacity(positions.len());

        for (idx, p) in positions.iter().enumerate() {
            let cx = ((p[0] - origin[0]) * inv_cell) as usize;
            let cy = ((p[1] - origin[1]) * inv_cell) as usize;
            let cz = ((p[2] - origin[2]) * inv_cell) as usize;
            let cx = cx.min(dims[0] - 1);
            let cy = cy.min(dims[1] - 1);
            let cz = cz.min(dims[2] - 1);
            let cell_idx = cx * dims[1] * dims[2] + cy * dims[2] + cz;
            cells[cell_idx].push(idx);
            assignments.push(cell_idx);
        }

        Self {
            cell_size,
            dims,
            origin,
            cells,
            assignments,
        }
    }

    /// Find all pairs within cutoff distance. Returns `(i, j, dist_sq)` where `i < j`.
    pub fn find_pairs(&self, positions: &[[f64; 3]]) -> Vec<(usize, usize, f64)> {
        let cutoff_sq = self.cell_size * self.cell_size;
        let mut pairs = Vec::new();

        if self.dims[0] == 0 {
            return pairs;
        }

        let [nx, ny, nz] = self.dims;

        for cx in 0..nx {
            for cy in 0..ny {
                for cz in 0..nz {
                    let cell_idx = cx * ny * nz + cy * nz + cz;

                    // Check same cell (self-pairs)
                    let atoms = &self.cells[cell_idx];
                    for a in 0..atoms.len() {
                        for b in (a + 1)..atoms.len() {
                            let i = atoms[a];
                            let j = atoms[b];
                            let dsq = dist_sq(&positions[i], &positions[j]);
                            if dsq <= cutoff_sq {
                                let (lo, hi) = if i < j { (i, j) } else { (j, i) };
                                pairs.push((lo, hi, dsq));
                            }
                        }
                    }

                    // Check 13 forward neighbors (to avoid double counting)
                    for &(dx, dy, dz) in &FORWARD_NEIGHBORS {
                        let nx2 = cx as isize + dx;
                        let ny2 = cy as isize + dy;
                        let nz2 = cz as isize + dz;
                        if nx2 < 0 || ny2 < 0 || nz2 < 0 {
                            continue;
                        }
                        let (nx2, ny2, nz2) = (nx2 as usize, ny2 as usize, nz2 as usize);
                        if nx2 >= nx || ny2 >= ny || nz2 >= nz {
                            continue;
                        }
                        let neighbor_idx = nx2 * ny * nz + ny2 * nz + nz2;
                        let neighbors = &self.cells[neighbor_idx];

                        for &i in atoms {
                            for &j in neighbors {
                                let dsq = dist_sq(&positions[i], &positions[j]);
                                if dsq <= cutoff_sq {
                                    let (lo, hi) = if i < j { (i, j) } else { (j, i) };
                                    pairs.push((lo, hi, dsq));
                                }
                            }
                        }
                    }
                }
            }
        }

        pairs
    }

    /// Find all neighbors of a specific atom within cutoff.
    /// Returns `(neighbor_index, dist_sq)`.
    pub fn neighbors_of(&self, atom_idx: usize, positions: &[[f64; 3]]) -> Vec<(usize, f64)> {
        let cutoff_sq = self.cell_size * self.cell_size;
        let mut result = Vec::new();

        if self.dims[0] == 0 || atom_idx >= self.assignments.len() {
            return result;
        }

        let cell_idx = self.assignments[atom_idx];
        let [_nx, ny, nz] = self.dims;
        let cx = cell_idx / (ny * nz);
        let rem = cell_idx % (ny * nz);
        let cy = rem / nz;
        let cz = rem % nz;

        // Check all 27 cells (self + 26 neighbors)
        for dx in -1isize..=1 {
            for dy in -1isize..=1 {
                for dz in -1isize..=1 {
                    let nx2 = cx as isize + dx;
                    let ny2 = cy as isize + dy;
                    let nz2 = cz as isize + dz;
                    if nx2 < 0 || ny2 < 0 || nz2 < 0 {
                        continue;
                    }
                    let (nx2, ny2, nz2) = (nx2 as usize, ny2 as usize, nz2 as usize);
                    if nx2 >= self.dims[0] || ny2 >= self.dims[1] || nz2 >= self.dims[2] {
                        continue;
                    }
                    let neighbor_cell = nx2 * ny * nz + ny2 * nz + nz2;
                    for &j in &self.cells[neighbor_cell] {
                        if j == atom_idx {
                            continue;
                        }
                        let dsq = dist_sq(&positions[atom_idx], &positions[j]);
                        if dsq <= cutoff_sq {
                            result.push((j, dsq));
                        }
                    }
                }
            }
        }

        result
    }

    /// Number of cells in the grid.
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Grid dimensions.
    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }

    /// Grid origin (minimum corner).
    pub fn origin(&self) -> [f64; 3] {
        self.origin
    }
}

/// 13 forward neighbors for half-shell enumeration.
const FORWARD_NEIGHBORS: [(isize, isize, isize); 13] = [
    (1, 0, 0), (0, 1, 0), (0, 0, 1),
    (1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1),
    (0, 1, 1), (0, 1, -1),
    (1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1),
];

#[inline]
fn dist_sq(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

// ── Tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_positions() {
        let cl = CellList::build(&[], &CellListConfig::default());
        assert_eq!(cl.cell_count(), 0);
        assert!(cl.find_pairs(&[]).is_empty());
    }

    #[test]
    fn single_atom_no_pairs() {
        let pos = [[0.0, 0.0, 0.0]];
        let cl = CellList::build(&pos, &CellListConfig::default());
        assert!(cl.find_pairs(&pos).is_empty());
    }

    #[test]
    fn two_atoms_within_cutoff() {
        let pos = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let cfg = CellListConfig { cutoff: 5.0, origin: None };
        let cl = CellList::build(&pos, &cfg);
        let pairs = cl.find_pairs(&pos);
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].0, 0);
        assert_eq!(pairs[0].1, 1);
        assert!((pairs[0].2 - 9.0).abs() < 1e-10); // dist_sq = 9
    }

    #[test]
    fn two_atoms_beyond_cutoff() {
        let pos = [[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]];
        let cfg = CellListConfig { cutoff: 5.0, origin: None };
        let cl = CellList::build(&pos, &cfg);
        let pairs = cl.find_pairs(&pos);
        assert!(pairs.is_empty());
    }

    #[test]
    fn matches_brute_force() {
        // Grid of 27 atoms at integer positions
        let mut pos = Vec::new();
        for x in 0..3 {
            for y in 0..3 {
                for z in 0..3 {
                    pos.push([x as f64 * 2.0, y as f64 * 2.0, z as f64 * 2.0]);
                }
            }
        }
        let cutoff = 3.0;
        let cfg = CellListConfig { cutoff, origin: None };
        let cl = CellList::build(&pos, &cfg);
        let cell_pairs = cl.find_pairs(&pos);

        // Brute force
        let cutoff_sq = cutoff * cutoff;
        let mut brute_pairs = Vec::new();
        for i in 0..pos.len() {
            for j in (i + 1)..pos.len() {
                let dsq = dist_sq(&pos[i], &pos[j]);
                if dsq <= cutoff_sq {
                    brute_pairs.push((i, j));
                }
            }
        }

        // Cell list should find all brute-force pairs
        let cell_set: std::collections::HashSet<(usize, usize)> =
            cell_pairs.iter().map(|&(i, j, _)| (i, j)).collect();
        for &(i, j) in &brute_pairs {
            assert!(cell_set.contains(&(i, j)), "Missing pair ({}, {})", i, j);
        }
        assert_eq!(cell_set.len(), brute_pairs.len());
    }

    #[test]
    fn neighbors_of_specific_atom() {
        let pos = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [100.0, 0.0, 0.0],
        ];
        let cfg = CellListConfig { cutoff: 5.0, origin: None };
        let cl = CellList::build(&pos, &cfg);
        let nbrs = cl.neighbors_of(0, &pos);
        let ids: Vec<usize> = nbrs.iter().map(|&(i, _)| i).collect();
        assert!(ids.contains(&1));
        assert!(!ids.contains(&2));
    }

    #[test]
    fn neighbors_of_empty() {
        let cl = CellList::build(&[], &CellListConfig::default());
        assert!(cl.neighbors_of(0, &[]).is_empty());
    }

    #[test]
    fn cutoff_boundary_included() {
        // Two atoms exactly at cutoff distance
        let pos = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let cfg = CellListConfig { cutoff: 5.0, origin: None };
        let cl = CellList::build(&pos, &cfg);
        let pairs = cl.find_pairs(&pos);
        assert_eq!(pairs.len(), 1);
    }

    #[test]
    fn large_system_no_panic() {
        let mut pos = Vec::new();
        for i in 0..100 {
            pos.push([i as f64 * 0.5, 0.0, 0.0]);
        }
        let cfg = CellListConfig { cutoff: 2.0, origin: None };
        let cl = CellList::build(&pos, &cfg);
        let pairs = cl.find_pairs(&pos);
        assert!(!pairs.is_empty());
    }

    #[test]
    fn dims_reasonable() {
        let pos = [[0.0, 0.0, 0.0], [10.0, 10.0, 10.0]];
        let cfg = CellListConfig { cutoff: 5.0, origin: None };
        let cl = CellList::build(&pos, &cfg);
        let dims = cl.dims();
        assert!(dims[0] >= 2);
        assert!(dims[1] >= 2);
        assert!(dims[2] >= 2);
    }

    #[test]
    fn custom_origin() {
        let pos = [[5.0, 5.0, 5.0], [7.0, 5.0, 5.0]];
        let cfg = CellListConfig {
            cutoff: 5.0,
            origin: Some([0.0, 0.0, 0.0]),
        };
        let cl = CellList::build(&pos, &cfg);
        let origin = cl.origin();
        assert_eq!(origin, [0.0, 0.0, 0.0]);
        // Should still find the pair
        let pairs = cl.find_pairs(&pos);
        assert_eq!(pairs.len(), 1);
    }

    #[test]
    fn neighbors_of_out_of_bounds_index() {
        let pos = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let cfg = CellListConfig::default();
        let cl = CellList::build(&pos, &cfg);
        // atom_idx=999 is out of bounds → should return empty, not panic
        let nbrs = cl.neighbors_of(999, &pos);
        assert!(nbrs.is_empty());
    }

    #[test]
    fn collinear_atoms_all_within_cutoff() {
        // 5 atoms in a line, each 1.0 apart, cutoff=5.0
        let pos: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let cfg = CellListConfig { cutoff: 5.0, origin: None };
        let cl = CellList::build(&pos, &cfg);
        let pairs = cl.find_pairs(&pos);
        // All 5*4/2 = 10 pairs should be within cutoff (max dist = 4.0 < 5.0)
        assert_eq!(pairs.len(), 10, "Expected 10 pairs, got {}", pairs.len());
    }

    #[test]
    fn symmetry_of_neighbors() {
        // If atom j is a neighbor of atom i, then atom i is a neighbor of atom j
        let pos = [
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [6.0, 0.0, 0.0],
        ];
        let cfg = CellListConfig { cutoff: 4.0, origin: None };
        let cl = CellList::build(&pos, &cfg);

        let nbrs_of_0: Vec<usize> = cl.neighbors_of(0, &pos).iter().map(|&(i, _)| i).collect();
        let nbrs_of_1: Vec<usize> = cl.neighbors_of(1, &pos).iter().map(|&(i, _)| i).collect();

        if nbrs_of_0.contains(&1) {
            assert!(nbrs_of_1.contains(&0), "Neighbor symmetry violated: 0 sees 1 but 1 doesn't see 0");
        }
    }

    #[test]
    fn cell_count_positive_for_nonempty() {
        let pos = [[1.0, 2.0, 3.0]];
        let cl = CellList::build(&pos, &CellListConfig::default());
        assert!(cl.cell_count() > 0);
    }
}
