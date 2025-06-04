use std::fs;

pub const EPS: f64 = 1e-6;

#[derive(Debug, Clone, Copy)]
pub struct Config {
    pub gamma: usize,
    pub dev: f64,
    pub tlimit: usize,
    pub frac_st: bool,
    pub frac_cap: bool,
    pub robust_capacity_cuts: bool,
}
pub struct RBRPData {
    pub name: String,
    pub num_vertices: usize,
    pub capacity: usize,
    pub d: Vec<isize>,
    pub p: Vec<isize>,
    pub q: Vec<isize>,
    pub distances: Vec<Vec<usize>>,
}

pub fn read_from_file(file_path: &str) -> RBRPData {
    let file = fs::read_to_string(file_path).unwrap();

    let mut contents = file.lines();
    
    assert!(contents.next().unwrap() == "#name");
    let name = contents.next().unwrap().to_string();

    contents.next(); // Skip the next line
    assert!(contents.next().unwrap() == "#numVertices");
    let num_vertices: usize = contents.next().unwrap().parse().unwrap();

    contents.next(); // Skip the next line
    assert!(contents.next().unwrap() == "#capacity");
    let capacity: usize = contents.next().unwrap().parse().unwrap();

    contents.next(); // Skip the next line
    assert!(contents.next().unwrap() == "#d/p/q");
    let mut d = Vec::new();
    let mut p = Vec::new();
    let mut q = Vec::new();

    for _ in 0..num_vertices {
        let line = contents.next().unwrap();
        let parts: Vec<isize> = line.split_ascii_whitespace()
            .map(|s| {s.parse().unwrap()})
            .collect();
        d.push(parts[1]);
        p.push(parts[2]);
        q.push(parts[3]);
    }

    contents.next(); // Skip the next line
    assert!(contents.next().unwrap() == "#distances");
    let mut distances = vec![vec![0; num_vertices]; num_vertices];
    contents.next(); // Skip the next line
    for i in 0..num_vertices {
        let line = contents.next().unwrap();
        let parts: Vec<usize> = line.split_whitespace()
            .map(|s| s.parse().unwrap())
            .collect();
        for j in 0..num_vertices {
            distances[i][j] = parts[j+1];
        }
    }

    RBRPData {
        name,
        num_vertices: num_vertices,
        capacity,
        d,
        p,
        q,
        distances,
    }
}


#[derive(Debug, Clone)]
struct Arc {
    to: usize,
    cost: f64,
    flow: f64
}

#[derive(Debug, Clone)]
pub struct Graph {
    vertices: Vec<usize>,
    adjacent: Vec<Vec<Arc>>,
    pub max_flow: f64,
    pub source_cut: Vec<usize>,
    pub target_cut: Vec<usize>,
}

impl Graph {
    pub fn new(data: Vec<(usize, usize, f64)>, max_v: usize) -> Self {
        let mut graph = Graph {
            vertices: (0..max_v).collect(),
            adjacent: vec![Vec::new(); max_v],
            max_flow: 0.0,
            source_cut: Vec::new(),
            target_cut: Vec::new(),
        };
        for (from, to, cost) in data {
            graph.adjacent[from].push(Arc { to, cost, flow: 0.0 });
        }

        graph
    }

    pub fn min_cut(&mut self, source: usize, target: usize) -> f64 {
        // Edmonds Karp algorithm
        let mut flow = 0.0;
        loop {
            let mut q = Vec::new();
            q.push(source);
            let mut pred = vec![None; self.vertices.len()];
            while q.len() > 0 && pred[target].is_none() {
                let c = q.remove(0);
                for arc in &self.adjacent[c] {
                    if pred[arc.to].is_none() && arc.to != source && arc.flow < arc.cost {
                        pred[arc.to] = Some(c);
                        q.push(arc.to);
                    }
                }
            }

            if !(pred[target].is_none()) {
                let mut df = f64::INFINITY;
                let mut curr = target;
                loop {
                    // println!("Current node: {}, Previous node: {}", curr, prev);
                    if pred[curr].is_none() {
                        break;
                    }
                    let prev = pred[curr].unwrap();
                    let arc = self.adjacent[prev].iter()
                        .find(|a| a.to == curr)
                        .unwrap();
                    df = df.min(arc.cost - arc.flow);
                    curr = prev;
                }

                let mut curr = target;
                loop {
                    if pred[curr].is_none() {
                        break;
                    }
                    let prev = pred[curr].unwrap();
                    self.adjacent[prev].iter_mut()
                        .find(|a| a.to == curr)
                        .unwrap()
                        .flow += df;
                    curr = prev;
                }
                flow += df;
            }
            else {
                break;
            }
        }
        self.max_flow = flow;
        // Find the source cut
        let mut q = Vec::new();
        q.push(source);
        let mut visited = vec![false; self.vertices.len()];
        while q.len() > 0 {
            let c = q.remove(0);
            visited[c] = true;
            for arc in &self.adjacent[c] {
                if visited[arc.to] == false && arc.flow < arc.cost {
                    q.push(arc.to);
                }
            }
        }

        for i in 0..self.vertices.len() {
            if visited[i] {
                self.source_cut.push(i);
            } else {
                self.target_cut.push(i);
            }
        }

        self.max_flow
    }

    pub fn reset(&mut self) {
        for arcs in &mut self.adjacent {
            for arc in arcs {
                arc.flow = 0.0;
            }
        }
        self.max_flow = 0.0;
        self.source_cut.clear();
        self.target_cut.clear();
    }

    pub fn print(&self) {
        println!("Graph with {} vertices", self.vertices.len());
        for (i, arcs) in self.adjacent.iter().enumerate() {
            println!("Vertex {}: {:?}", i, arcs);
        }
    }
}

// Some basic testing
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_from_file() {
        let data = read_from_file("instances/01-Bari30_RBRP.txt");
        assert_eq!(data.name, "01-Bari30_RBRP");
        assert_eq!(data.num_vertices, 13);
        assert_eq!(data.capacity, 30);
        assert_eq!(data.d.len(), 13);
        assert_eq!(data.p.len(), 13);
        assert_eq!(data.q.len(), 13);
        assert_eq!(data.distances.len(), 13);
    }

    #[test]
    fn test_graph_min_cut_1() {
        let data = vec![
            (0, 1, 15.0),
            (0, 2, 5.0),
            (1, 2, 15.0),
            (1, 3, 10.0),
            (2, 3, 15.0),
            (3, 4, 14.0),
        ];
        let mut graph = Graph::new(data, 5);
        let flow = graph.min_cut(0, 4);
        assert_eq!(flow, 14.0);
        assert_eq!(graph.source_cut, vec![0, 1, 2, 3]);
        assert_eq!(graph.target_cut, vec![4]);
    }

    #[test]
    fn test_graph_min_cut_2() {
        let data = vec![
            (0, 1, 15.0),
            (0, 2, 5.0),
            (1, 2, 10.0),
            (1, 3, 10.0),
            (2, 3, 20.0),
            (3, 4, 30.0),
        ];
        let mut graph = Graph::new(data, 5);
        let flow = graph.min_cut(0, 4);
        assert_eq!(flow, 20.0);
        assert_eq!(graph.source_cut, vec![0]);
        assert_eq!(graph.target_cut, vec![1, 2, 3, 4]);
    }
}