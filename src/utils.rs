use std::fs;

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
}