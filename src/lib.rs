//! # Step-Back Topology
//!
//! The Step-Back Operator: β₁ = E - V + C
//! Topological data analysis via simplicial complexes, Betti numbers, and fishing simulations.

// ── graph ───────────────────────────────────────────────────────────────────

/// An undirected graph represented as adjacency sets.
#[derive(Debug, Clone)]
pub struct Graph {
    vertices: usize,
    edges: Vec<(usize, usize)>,
}

impl Graph {
    pub fn new(vertices: usize) -> Self {
        Self {
            vertices,
            edges: Vec::new(),
        }
    }

    pub fn add_edge(&mut self, a: usize, b: usize) {
        if a < self.vertices && b < self.vertices && a != b {
            let edge = if a < b { (a, b) } else { (b, a) };
            if !self.edges.contains(&edge) {
                self.edges.push(edge);
            }
        }
    }

    pub fn vertex_count(&self) -> usize {
        self.vertices
    }

    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Count connected components using union-find.
    pub fn connected_components(&self) -> usize {
        let mut parent: Vec<usize> = (0..self.vertices).collect();

        fn find(parent: &mut [usize], x: usize) -> usize {
            if parent[x] != x {
                parent[x] = find(parent, parent[x]);
            }
            parent[x]
        }

        for &(a, b) in &self.edges {
            let ra = find(&mut parent, a);
            let rb = find(&mut parent, b);
            if ra != rb {
                parent[ra] = rb;
            }
        }

        let mut roots = std::collections::HashSet::new();
        for i in 0..self.vertices {
            roots.insert(find(&mut parent, i));
        }
        roots.len()
    }

    /// Build an adjacency list representation.
    pub fn adjacency(&self) -> Vec<Vec<usize>> {
        let mut adj = vec![vec![]; self.vertices];
        for &(a, b) in &self.edges {
            adj[a].push(b);
            adj[b].push(a);
        }
        adj
    }
}

// ── betti ───────────────────────────────────────────────────────────────────

/// A simplicial complex with simplices up to dimension 2.
#[derive(Debug, Clone)]
pub struct SimplicialComplex {
    /// Vertices (0-simplices).
    pub vertices: usize,
    /// Edges (1-simplices) as pairs.
    pub edges: Vec<(usize, usize)>,
    /// Triangles (2-simplices) as triples (sorted).
    pub triangles: Vec<(usize, usize, usize)>,
}

impl SimplicialComplex {
    pub fn new(vertices: usize) -> Self {
        Self {
            vertices,
            edges: Vec::new(),
            triangles: Vec::new(),
        }
    }

    pub fn add_edge(&mut self, a: usize, b: usize) {
        let edge = if a < b { (a, b) } else { (b, a) };
        if !self.edges.contains(&edge) {
            self.edges.push(edge);
        }
    }

    pub fn add_triangle(&mut self, a: usize, b: usize, c: usize) {
        let mut tri = [a, b, c];
        tri.sort();
        let tri = (tri[0], tri[1], tri[2]);
        if !self.triangles.contains(&tri) {
            self.triangles.push(tri);
            // Ensure edges exist
            self.add_edge(tri.0, tri.1);
            self.add_edge(tri.0, tri.2);
            self.add_edge(tri.1, tri.2);
        }
    }

    /// Compute Betti numbers β₀, β₁, β₂.
    ///
    /// β₀ = connected components of the 1-skeleton
    /// β₁ = E - V + C + T  (E=edges, V=vertices, C=components, T=triangles, in the cycle space formula)
    ///     Actually: β₁ = E - V + C for the 1-skeleton, then subtract triangles that fill cycles.
    ///     With triangles: β₁ = rank H₁ = cycles - filled cycles.
    ///     Simplified: β₁ = E - V + C (edges that form independent cycles)
    /// β₂ = triangles that are not boundaries of tetrahedra (here, just count independent voids)
    pub fn betti_numbers(&self) -> (usize, usize, usize) {
        // β₀: connected components
        let g = {
            let mut g = Graph::new(self.vertices);
            for &(a, b) in &self.edges {
                g.add_edge(a, b);
            }
            g
        };
        let c = g.connected_components();
        let v = self.vertices;
        let e = self.edges.len();
        let t = self.triangles.len();

        // β₁ = E - V + C - T (each filled triangle reduces β₁ by 1)
        // For a pure 1-skeleton: β₁ = E - V + C
        let beta1 = (e as i64 + c as i64 - v as i64 - t as i64).max(0) as usize;

        // β₂: for our 2-complex, it's 0 unless we have enclosed voids
        // Simplified: no 3-simplices, so β₂ = 0 in basic cases
        let beta2 = 0;

        (c, beta1, beta2)
    }
}

// ── fishing ─────────────────────────────────────────────────────────────────

/// A hook in the fishing simulation — a boolean probe on a simplicial complex.
#[derive(Debug, Clone)]
pub struct Hook {
    /// Which vertex the hook is attached to.
    pub vertex: usize,
    /// Whether the hook "catches" (true = occupied).
    pub caught: bool,
}

/// Simulate longline fishing as topological data analysis.
///
/// Deploy `hooks` across vertices of a complex. The topology of the
/// catch pattern reveals the Betti numbers of the underlying complex.
pub fn fish(complex: &SimplicialComplex, hooks: &[Hook]) -> FishingResult {
    let caught: Vec<usize> = hooks.iter().filter(|h| h.caught).map(|h| h.vertex).collect();
    let (beta0, beta1, beta2) = complex.betti_numbers();

    FishingResult {
        hooks_deployed: hooks.len(),
        hooks_caught: caught.len(),
        catch_vertices: caught,
        beta0,
        beta1,
        beta2,
    }
}

/// Result of a fishing simulation.
#[derive(Debug, Clone)]
pub struct FishingResult {
    pub hooks_deployed: usize,
    pub hooks_caught: usize,
    pub catch_vertices: Vec<usize>,
    pub beta0: usize,
    pub beta1: usize,
    pub beta2: usize,
}

// ── profile ─────────────────────────────────────────────────────────────────

/// An emergent distribution detected from boolean hook data.
#[derive(Debug, Clone, PartialEq)]
pub enum Distribution {
    /// Uniform catch across vertices.
    Uniform,
    /// Clustered in one region.
    Clustered { center: usize },
    /// No catches.
    Empty,
    /// Sparse scattered pattern.
    Sparse,
}

/// Detect the emergent distribution from a fishing result.
pub fn detect_distribution(result: &FishingResult) -> Distribution {
    if result.hooks_caught == 0 {
        return Distribution::Empty;
    }
    if result.hooks_deployed == 0 {
        return Distribution::Empty;
    }
    let ratio = result.hooks_caught as f64 / result.hooks_deployed as f64;
    if ratio >= 0.8 {
        return Distribution::Uniform;
    }
    if ratio < 0.2 {
        return Distribution::Sparse;
    }
    // Check clustering: are most catches near one vertex?
    let vertices = &result.catch_vertices;
    if !vertices.is_empty() {
        let mean = vertices.iter().sum::<usize>() as f64 / vertices.len() as f64;
        let variance = vertices.iter().map(|&v| {
            let diff = v as f64 - mean;
            diff * diff
        }).sum::<f64>() / vertices.len() as f64;
        if variance < 0.5 {
            // If most hooks caught, it's uniform, not clustered
            if ratio > 0.5 {
                return Distribution::Uniform;
            }
            return Distribution::Clustered { center: vertices[0] };
        }
    }
    Distribution::Sparse
}

// ── tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_graph_vertex_edge_count() {
        let mut g = Graph::new(4);
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        assert_eq!(g.vertex_count(), 4);
        assert_eq!(g.edge_count(), 2);
    }

    #[test]
    fn test_graph_connected_components_single() {
        let mut g = Graph::new(4);
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        g.add_edge(2, 3);
        assert_eq!(g.connected_components(), 1);
    }

    #[test]
    fn test_graph_connected_components_multiple() {
        let mut g = Graph::new(6);
        g.add_edge(0, 1);
        g.add_edge(2, 3);
        // 4 and 5 isolated
        assert_eq!(g.connected_components(), 4);
    }

    #[test]
    fn test_betti_cycle_graph() {
        // Cycle graph: 4 vertices, 4 edges → β₁ = E - V + C = 4 - 4 + 1 = 1
        let mut cx = SimplicialComplex::new(4);
        cx.add_edge(0, 1);
        cx.add_edge(1, 2);
        cx.add_edge(2, 3);
        cx.add_edge(3, 0);
        let (b0, b1, b2) = cx.betti_numbers();
        assert_eq!(b0, 1);
        assert_eq!(b1, 1);
        assert_eq!(b2, 0);
    }

    #[test]
    fn test_betti_tree() {
        // Tree: 4 vertices, 3 edges → β₁ = 3 - 4 + 1 = 0
        let mut tree = SimplicialComplex::new(4);
        tree.add_edge(0, 1);
        tree.add_edge(1, 2);
        tree.add_edge(1, 3);
        let (b0, b1, _b2) = tree.betti_numbers();
        assert_eq!(b0, 1);
        assert_eq!(b1, 0);
    }

    #[test]
    fn test_betti_two_cycles() {
        // Figure-8: 6 vertices, 7 edges → β₁ = 7 - 6 + 1 = 2
        let mut fig = SimplicialComplex::new(6);
        // Cycle 1: 0-1-2-0
        fig.add_edge(0, 1);
        fig.add_edge(1, 2);
        fig.add_edge(2, 0);
        // Cycle 2: 0-3-4-5-0
        fig.add_edge(0, 3);
        fig.add_edge(3, 4);
        fig.add_edge(4, 5);
        fig.add_edge(5, 0);
        let (_, b1, _) = fig.betti_numbers();
        assert_eq!(b1, 2);
    }

    #[test]
    fn test_betti_triangle_fills_cycle() {
        // Triangle 0-1-2 with all 3 edges + face → β₁ = 3 - 3 + 1 - 1 = 0
        let mut cx = SimplicialComplex::new(3);
        cx.add_triangle(0, 1, 2);
        let (b0, b1, b2) = cx.betti_numbers();
        assert_eq!(b0, 1);
        assert_eq!(b1, 0);
        assert_eq!(b2, 0);
    }

    #[test]
    fn test_fishing_simulation() {
        let mut cx = SimplicialComplex::new(4);
        cx.add_edge(0, 1);
        cx.add_edge(1, 2);
        cx.add_edge(2, 3);
        cx.add_edge(3, 0);
        let hooks = vec![
            Hook { vertex: 0, caught: true },
            Hook { vertex: 1, caught: false },
            Hook { vertex: 2, caught: true },
            Hook { vertex: 3, caught: false },
        ];
        let result = fish(&cx, &hooks);
        assert_eq!(result.hooks_deployed, 4);
        assert_eq!(result.hooks_caught, 2);
        assert_eq!(result.beta1, 1);
    }

    #[test]
    fn test_fishing_empty_catch() {
        let cx = SimplicialComplex::new(3);
        let hooks = vec![
            Hook { vertex: 0, caught: false },
            Hook { vertex: 1, caught: false },
        ];
        let result = fish(&cx, &hooks);
        assert_eq!(result.hooks_caught, 0);
    }

    #[test]
    fn test_distribution_empty() {
        let result = FishingResult {
            hooks_deployed: 4,
            hooks_caught: 0,
            catch_vertices: vec![],
            beta0: 1, beta1: 0, beta2: 0,
        };
        assert_eq!(detect_distribution(&result), Distribution::Empty);
    }

    #[test]
    fn test_distribution_uniform() {
        let result = FishingResult {
            hooks_deployed: 5,
            hooks_caught: 4,
            catch_vertices: vec![0, 1, 2, 3],
            beta0: 1, beta1: 0, beta2: 0,
        };
        assert_eq!(detect_distribution(&result), Distribution::Uniform);
    }

    #[test]
    fn test_distribution_sparse() {
        let result = FishingResult {
            hooks_deployed: 10,
            hooks_caught: 1,
            catch_vertices: vec![5],
            beta0: 1, beta1: 0, beta2: 0,
        };
        assert_eq!(detect_distribution(&result), Distribution::Sparse);
    }

    #[test]
    fn test_distribution_clustered() {
        let result = FishingResult {
            hooks_deployed: 6,
            hooks_caught: 3,
            catch_vertices: vec![2, 2, 3], // tight cluster
            beta0: 1, beta1: 0, beta2: 0,
        };
        match detect_distribution(&result) {
            Distribution::Clustered { center } => assert_eq!(center, 2),
            other => panic!("expected Clustered, got {:?}", other),
        }
    }
}
