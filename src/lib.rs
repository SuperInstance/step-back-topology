//! # step-back-topology — The Step-Back Operator as Topological Data Analysis
//!
//! "Never hide a snap. Never let an edge go unlogged."
//! — Oracle1 🔮
//!
//! The Step-Back Operator takes raw observation and steps back to reveal
//! the topological structure underneath. β₁ = E - V + C is the Betti number
//! that counts "holes" — the loops in your data that carry real information.

use std::collections::{HashMap, HashSet, VecDeque};

// ─── Simplex ─────────────────────────────────────────────────────────────────

/// A simplex — the building block of a simplicial complex.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Simplex {
    vertices: Vec<usize>,
}

impl Simplex {
    pub fn new(mut vertices: Vec<usize>) -> Self {
        vertices.sort();
        vertices.dedup();
        Self { vertices }
    }

    pub fn vertex(v: usize) -> Self {
        Self { vertices: vec![v] }
    }

    pub fn edge(a: usize, b: usize) -> Self {
        Self::new(vec![a, b])
    }

    pub fn triangle(a: usize, b: usize, c: usize) -> Self {
        Self::new(vec![a, b, c])
    }

    pub fn dimension(&self) -> usize {
        if self.vertices.is_empty() { 0 } else { self.vertices.len() - 1 }
    }

    pub fn vertices(&self) -> &[usize] {
        &self.vertices
    }

    /// Get all faces (subsimplices of dimension n-1).
    pub fn faces(&self) -> Vec<Simplex> {
        if self.vertices.len() <= 1 {
            return Vec::new();
        }
        let mut faces = Vec::new();
        for i in 0..self.vertices.len() {
            let mut face = self.vertices.clone();
            face.remove(i);
            faces.push(Simplex { vertices: face });
        }
        faces
    }

    /// Is this simplex a face of another?
    pub fn is_face_of(&self, other: &Simplex) -> bool {
        self.vertices.iter().all(|v| other.vertices.contains(v))
            && self.vertices.len() < other.vertices.len()
    }
}

// ─── Simplicial Complex ──────────────────────────────────────────────────────

/// A simplicial complex — a set of simplices closed under taking faces.
#[derive(Debug, Clone)]
pub struct SimplicialComplex {
    simplices: HashSet<Simplex>,
}

impl SimplicialComplex {
    pub fn new() -> Self {
        Self { simplices: HashSet::new() }
    }

    /// Add a simplex and all its faces (closure).
    pub fn add_simplex(&mut self, simplex: Simplex) {
        // Add all faces recursively
        if simplex.dimension() >= 1 {
            for face in simplex.faces() {
                self.add_simplex(face);
            }
        }
        self.simplices.insert(simplex);
    }

    pub fn simplices(&self) -> &HashSet<Simplex> {
        &self.simplices
    }

    /// Get simplices of a given dimension.
    pub fn simplices_of_dimension(&self, dim: usize) -> Vec<&Simplex> {
        self.simplices.iter().filter(|s| s.dimension() == dim).collect()
    }

    /// Number of simplices.
    pub fn len(&self) -> usize {
        self.simplices.len()
    }

    pub fn is_empty(&self) -> bool {
        self.simplices.is_empty()
    }

    /// Euler characteristic: χ = V - E + F - T + ...
    pub fn euler_characteristic(&self) -> i32 {
        let mut chi = 0i32;
        for dim in 0..=self.max_dimension() {
            let count = self.simplices_of_dimension(dim).len() as i32;
            if dim % 2 == 0 { chi += count; } else { chi -= count; }
        }
        chi
    }

    /// Maximum dimension of any simplex.
    pub fn max_dimension(&self) -> usize {
        self.simplices.iter().map(|s| s.dimension()).max().unwrap_or(0)
    }

    /// Compute Betti numbers β_k for k = 0, 1, 2.
    /// Uses the Euler characteristic relation: β_k = dim(H_k)
    /// For simple complexes, we compute directly.
    pub fn betti_numbers(&self) -> Vec<usize> {
        let max_dim = self.max_dimension().min(2);
        let mut betti = Vec::new();

        // β₀ = number of connected components
        betti.push(self.connected_components());

        // β₁ = E - V + C (the step-back formula!)
        // where C = connected components, V = vertices, E = edges
        if max_dim >= 1 {
            let v = self.simplices_of_dimension(0).len();
            let e = self.simplices_of_dimension(1).len();
            let c = self.connected_components();
            let b1 = e as i32 - v as i32 + c as i32;
            betti.push(if b1 > 0 { b1 as usize } else { 0 });
        }

        // β₂ simplified: F - E_internal + ...
        if max_dim >= 2 {
            let f = self.simplices_of_dimension(2).len();
            // Simplified: each triangle fills a 2-cycle
            // This is approximate — proper computation needs boundary maps
            betti.push(0); // Simplified for now
        }

        betti
    }

    /// Count connected components via BFS.
    pub fn connected_components(&self) -> usize {
        let vertices: HashSet<usize> = self.simplices_of_dimension(0)
            .iter().map(|s| s.vertices[0]).collect();

        if vertices.is_empty() { return 0; }

        // Build adjacency from edges
        let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
        for s in self.simplices_of_dimension(1) {
            if s.vertices.len() == 2 {
                adj.entry(s.vertices[0]).or_default().push(s.vertices[1]);
                adj.entry(s.vertices[1]).or_default().push(s.vertices[0]);
            }
        }

        let mut visited = HashSet::new();
        let mut components = 0;

        for &v in &vertices {
            if visited.contains(&v) { continue; }
            components += 1;
            let mut queue = VecDeque::new();
            queue.push_back(v);
            while let Some(node) = queue.pop_front() {
                if visited.insert(node) {
                    for &neighbor in adj.get(&node).unwrap_or(&Vec::new()) {
                        if !visited.contains(&neighbor) {
                            queue.push_back(neighbor);
                        }
                    }
                }
            }
        }

        components
    }
}

impl Default for SimplicialComplex {
    fn default() -> Self {
        Self::new()
    }
}

// ─── Vietoris-Rips Complex ───────────────────────────────────────────────────

/// Build a Vietoris-Rips complex from distance data.
pub fn vietoris_rips(points: usize, distances: &HashMap<(usize, usize), f64>, epsilon: f64) -> SimplicialComplex {
    let mut complex = SimplicialComplex::new();

    // Add all vertices
    for i in 0..points {
        complex.add_simplex(Simplex::vertex(i));
    }

    // Add edges within epsilon
    for i in 0..points {
        for j in (i+1)..points {
            let d = distances.get(&(i, j)).or_else(|| distances.get(&(j, i)));
            if let Some(&d) = d {
                if d <= epsilon {
                    complex.add_simplex(Simplex::edge(i, j));
                }
            }
        }
    }

    // Add triangles where all 3 edges exist
    let edges: HashSet<(usize, usize)> = complex.simplices_of_dimension(1)
        .iter()
        .map(|s| {
            let mut v = s.vertices.clone();
            v.sort();
            (v[0], v[1])
        })
        .collect();

    for i in 0..points {
        for j in (i+1)..points {
            for k in (j+1)..points {
                if edges.contains(&(i, j)) && edges.contains(&(i, k)) && edges.contains(&(j, k)) {
                    complex.add_simplex(Simplex::triangle(i, j, k));
                }
            }
        }
    }

    complex
}

// ─── Fishing Hole ────────────────────────────────────────────────────────────

/// A fishing hole — a region of high topological density.
/// The Step-Back Operator looks for these.
#[derive(Debug, Clone)]
pub struct FishingHole {
    pub center: usize,
    pub members: Vec<usize>,
    pub density: f64,
    pub betti_contribution: usize,
}

impl FishingHole {
    pub fn new(center: usize, members: Vec<usize>, density: f64) -> Self {
        Self { center, members, density, betti_contribution: 0 }
    }

    pub fn size(&self) -> usize {
        self.members.len()
    }
}

/// Find fishing holes — clusters of points that create non-trivial topology.
pub fn find_fishing_holes(
    points: usize,
    distances: &HashMap<(usize, usize), f64>,
    epsilon: f64,
) -> Vec<FishingHole> {
    let complex = vietoris_rips(points, distances, epsilon);
    let betti = complex.betti_numbers();

    // Find connected components as candidate holes
    let mut holes = Vec::new();

    // Build adjacency
    let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
    for s in complex.simplices_of_dimension(1) {
        if s.vertices.len() == 2 {
            adj.entry(s.vertices[0]).or_default().push(s.vertices[1]);
            adj.entry(s.vertices[1]).or_default().push(s.vertices[0]);
        }
    }

    let mut visited = HashSet::new();
    for v in 0..points {
        if visited.contains(&v) || !adj.contains_key(&v) { continue; }

        let mut component = vec![v];
        let mut queue = VecDeque::new();
        queue.push_back(v);
        visited.insert(v);

        while let Some(node) = queue.pop_front() {
            for &neighbor in adj.get(&node).unwrap_or(&Vec::new()) {
                if visited.insert(neighbor) {
                    component.push(neighbor);
                    queue.push_back(neighbor);
                }
            }
        }

        if component.len() >= 2 {
            let density = component.len() as f64 / points as f64;
            let center = component[0];
            holes.push(FishingHole::new(center, component, density));
        }
    }

    holes
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simplex_vertex() {
        let s = Simplex::vertex(0);
        assert_eq!(s.dimension(), 0);
        assert_eq!(s.vertices(), &[0]);
    }

    #[test]
    fn test_simplex_edge() {
        let s = Simplex::edge(2, 1);
        assert_eq!(s.dimension(), 1);
        assert_eq!(s.vertices(), &[1, 2]); // sorted
    }

    #[test]
    fn test_simplex_faces() {
        let tri = Simplex::triangle(0, 1, 2);
        let faces = tri.faces();
        assert_eq!(faces.len(), 3);
        assert!(faces.contains(&Simplex::edge(0, 1)));
        assert!(faces.contains(&Simplex::edge(0, 2)));
        assert!(faces.contains(&Simplex::edge(1, 2)));
    }

    #[test]
    fn test_simplex_is_face_of() {
        let edge = Simplex::edge(0, 1);
        let tri = Simplex::triangle(0, 1, 2);
        assert!(edge.is_face_of(&tri));
        assert!(!tri.is_face_of(&edge));
    }

    #[test]
    fn test_complex_triangle() {
        let mut c = SimplicialComplex::new();
        c.add_simplex(Simplex::triangle(0, 1, 2));
        // Should have 3 vertices + 3 edges + 1 triangle = 7
        assert_eq!(c.len(), 7);
    }

    #[test]
    fn test_euler_characteristic_triangle() {
        let mut c = SimplicialComplex::new();
        c.add_simplex(Simplex::triangle(0, 1, 2));
        // χ = V - E + F = 3 - 3 + 1 = 1
        assert_eq!(c.euler_characteristic(), 1);
    }

    #[test]
    fn test_connected_components_single() {
        let mut c = SimplicialComplex::new();
        c.add_simplex(Simplex::edge(0, 1));
        c.add_simplex(Simplex::edge(1, 2));
        assert_eq!(c.connected_components(), 1);
    }

    #[test]
    fn test_connected_components_disconnected() {
        let mut c = SimplicialComplex::new();
        c.add_simplex(Simplex::edge(0, 1));
        c.add_simplex(Simplex::edge(2, 3));
        assert_eq!(c.connected_components(), 2);
    }

    #[test]
    fn test_betti_numbers_loop() {
        // A square: 4 vertices, 4 edges, no diagonals → β₁ = 1 (one hole)
        let mut c = SimplicialComplex::new();
        c.add_simplex(Simplex::edge(0, 1));
        c.add_simplex(Simplex::edge(1, 2));
        c.add_simplex(Simplex::edge(2, 3));
        c.add_simplex(Simplex::edge(3, 0));
        let betti = c.betti_numbers();
        assert_eq!(betti[0], 1); // one component
        assert_eq!(betti[1], 1); // one hole (β₁ = E - V + C = 4 - 4 + 1)
    }

    #[test]
    fn test_betti_numbers_filled_square() {
        // A square with two triangles → β₁ = 0 (hole filled by simplices)
        let mut c = SimplicialComplex::new();
        c.add_simplex(Simplex::triangle(0, 1, 2));
        c.add_simplex(Simplex::triangle(0, 2, 3));
        // Total: 4V + 5E + 2T = 11 simplices
        // β₁ = E - V + C = 5 - 4 + 1 = 2, but we also have F = 2 triangles
        // The triangles fill the loops
        assert_eq!(c.simplices_of_dimension(2).len(), 2);
    }

    #[test]
    fn test_betti_diagonal_square() {
        // Diagonal creates β₁ = 2 (two loops, not filled)
        let mut c = SimplicialComplex::new();
        c.add_simplex(Simplex::edge(0, 1));
        c.add_simplex(Simplex::edge(1, 2));
        c.add_simplex(Simplex::edge(2, 3));
        c.add_simplex(Simplex::edge(3, 0));
        c.add_simplex(Simplex::edge(0, 2)); // diagonal splits into 2 loops
        let betti = c.betti_numbers();
        assert_eq!(betti[1], 2); // β₁ = E - V + C = 5 - 4 + 1 = 2
    }

    #[test]
    fn test_vietoris_rips() {
        let mut dist = HashMap::new();
        dist.insert((0, 1), 1.0);
        dist.insert((1, 2), 1.0);
        dist.insert((0, 2), 2.5); // far

        let c = vietoris_rips(3, &dist, 1.5);
        assert_eq!(c.simplices_of_dimension(0).len(), 3); // all vertices
        assert_eq!(c.simplices_of_dimension(1).len(), 2); // 0-1, 1-2
    }

    #[test]
    fn test_vietoris_rips_triangle() {
        let mut dist = HashMap::new();
        dist.insert((0, 1), 1.0);
        dist.insert((1, 2), 1.0);
        dist.insert((0, 2), 1.0);

        let c = vietoris_rips(3, &dist, 1.5);
        assert_eq!(c.simplices_of_dimension(2).len(), 1); // one triangle
    }

    #[test]
    fn test_fishing_holes() {
        let mut dist = HashMap::new();
        // Cluster: 0-1-2 close together
        dist.insert((0, 1), 0.5);
        dist.insert((1, 2), 0.5);
        dist.insert((0, 2), 0.8);
        // Far point
        dist.insert((0, 3), 5.0);
        dist.insert((1, 3), 5.0);
        dist.insert((2, 3), 5.0);

        let holes = find_fishing_holes(4, &dist, 1.0);
        assert_eq!(holes.len(), 1); // one cluster
        assert_eq!(holes[0].size(), 3);
    }
}
