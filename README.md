# Step-Back Topology

[![crates.io](https://img.shields.io/crates/v/step-back-topology.svg)](https://crates.io/crates/step-back-topology)
[![docs.rs](https://docs.rs/step-back-topology/badge.svg)](https://docs.rs/step-back-topology)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

> **The Step-Back Operator as full topological data analysis — simplicial complexes, Betti numbers, Vietoris-Rips complexes, and fishing hole detection.**

---

*"Never hide a snap. Never let an edge go unlogged."*
— Oracle1 🔮

## The Problem

The `step-back-operator` crate gives you the formula β₁ = E − V + C. But what if you need the full machinery of topological data analysis? Simplicial complexes, Vietoris-Rips complexes, connected components via BFS, and the ability to find "fishing holes" — clusters of points creating non-trivial topology?

## Why This Exists

Step-Back Topology extends the core Betti number computation into a full TDA toolkit:
- **Simplices** (vertices, edges, triangles) as first-class objects
- **Simplicial Complexes** with automatic face closure
- **Euler Characteristic** computation
- **Betti Numbers** β₀, β₁, β₂ from the simplicial structure
- **Vietoris-Rips Complex** construction from distance data
- **Fishing Hole Detection** — finding regions of high topological density

## Architecture

```
  Points + Distances ──→ Vietoris-Rips ──→ Simplicial Complex
                                                    │
                                          ┌─────────┼─────────┐
                                          │         │         │
                                     Euler χ    Betti βk   Fishing
                                                   │       Holes
                                              β₀ = C
                                         β₁ = E - V + C
                                         β₂ (simplified)
```

## Installation

```toml
[dependencies]
step-back-topology = "0.1"
```

## API Reference

### `Simplex`

Building blocks of a simplicial complex:

```rust
use step_back_topology::Simplex;

let v = Simplex::vertex(0);           // 0-dimensional
let e = Simplex::edge(0, 1);          // 1-dimensional
let t = Simplex::triangle(0, 1, 2);   // 2-dimensional

assert_eq!(t.faces().len(), 3); // 3 edges
assert!(e.is_face_of(&t));
```

### `SimplicialComplex`

A set of simplices closed under taking faces:

```rust
use step_back_topology::SimplicialComplex;

let mut c = SimplicialComplex::new();
c.add_simplex(Simplex::triangle(0, 1, 2));
// Automatically adds 3 vertices + 3 edges = 7 total
assert_eq!(c.len(), 7);
assert_eq!(c.euler_characteristic(), 1); // χ = V - E + F = 3 - 3 + 1
```

### Betti Numbers

```rust
// Square with no diagonal: β₁ = 1 (one hole)
let mut c = SimplicialComplex::new();
c.add_simplex(Simplex::edge(0, 1));
c.add_simplex(Simplex::edge(1, 2));
c.add_simplex(Simplex::edge(2, 3));
c.add_simplex(Simplex::edge(3, 0));

let betti = c.betti_numbers();
assert_eq!(betti[0], 1); // one component
assert_eq!(betti[1], 1); // one hole (β₁ = 4 - 4 + 1)
```

### `vietoris_rips`

Build a complex from distance data:

```rust
use step_back_topology::vietoris_rips;
use std::collections::HashMap;

let mut dist = HashMap::new();
dist.insert((0, 1), 1.0);
dist.insert((1, 2), 1.0);
dist.insert((0, 2), 1.0);

let complex = vietoris_rips(3, &dist, 1.5);
assert_eq!(complex.simplices_of_dimension(2).len(), 1); // one triangle
```

### `find_fishing_holes`

Find clusters creating non-trivial topology:

```rust
use step_back_topology::find_fishing_holes;

let mut dist = HashMap::new();
dist.insert((0, 1), 0.5);
dist.insert((1, 2), 0.5);
dist.insert((0, 2), 0.8);
dist.insert((0, 3), 5.0); // far point

let holes = find_fishing_holes(4, &dist, 1.0);
assert_eq!(holes.len(), 1); // one cluster of 3 points
```

## Usage Examples

### Example 1: Analyze a Social Network

```rust
use step_back_topology::*;
use std::collections::HashMap;

let mut dist = HashMap::new();
// Close relationships
dist.insert((0, 1), 0.3);
dist.insert((1, 2), 0.4);
dist.insert((2, 3), 0.3);
dist.insert((3, 0), 0.5);
// Distant member
dist.insert((0, 4), 2.0);

let complex = vietoris_rips(5, &dist, 1.0);
let betti = complex.betti_numbers();
println!("Components: {}, Loops: {}", betti[0], betti[1]);
```

### Example 2: Euler Characteristic

```rust
use step_back_topology::*;

let mut c = SimplicialComplex::new();
c.add_simplex(Simplex::triangle(0, 1, 2)); // χ = 1 (solid triangle)
assert_eq!(c.euler_characteristic(), 1);
```

### Example 3: Connected Components

```rust
use step_back_topology::*;

let mut c = SimplicialComplex::new();
c.add_simplex(Simplex::edge(0, 1));
c.add_simplex(Simplex::edge(2, 3));
assert_eq!(c.connected_components(), 2); // two separate groups
```

## Mathematical Background

**Betti Numbers** βk count the number of k-dimensional "holes":
- β₀ = connected components
- β₁ = independent loops (β₁ = E − V + C, the step-back formula!)
- β₂ = enclosed cavities

**Euler Characteristic**: χ = Σᵢ₌₀ⁿ (-1)ⁱ |simplices_of_dim(i)|

**Vietoris-Rips Complex**: Given points with pairwise distances, connect all points within distance ε, and fill in all simplices whose edges all exist.

## Performance

| Operation | Complexity |
|-----------|-----------|
| Add simplex (with closure) | O(2^k × k) |
| Euler characteristic | O(n × max_dim) |
| Connected components (BFS) | O(V + E) |
| Betti numbers | O(V + E) |
| Vietoris-Rips | O(n³) worst case |

## License

Licensed under the [MIT License](LICENSE).

## Contributing

1. Fork the repository
2. Create a feature branch
3. Write tests
4. Push and open a Pull Request
