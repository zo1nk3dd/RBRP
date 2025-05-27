use std::ops::Index;

use grb::{attribute::{ConstrDoubleAttr::Slack, ModelDoubleAttr::ObjVal, VarDoubleAttr::{UB, X}}, callback::*, parameter::IntParam::{LazyConstraints, Threads}, prelude::*};

use crate::utils::*;

pub struct CompactModel {
    model: Model,
    gamma: usize,
    deviation: Option<f64>,
    data: RBRPData,

    // Variables
    x: Vec<Vec<Var>>,
    _u: Vec<Var>,
    _theta_p: Vec<Vec<Var>>,
    _theta_m: Vec<Vec<Var>>,
    theta_d: Vec<Var>,

    clique_cuts: Vec<Constr>,
}

impl CompactModel {
    pub fn new(filename: &str, gamma: usize, deviation: Option<f64>) -> Self {
        // Read data from file
        let data = read_from_file(filename);
        let mut model = Model::new("RBRP Compact Model").unwrap();

        // Create variables
        let x: Vec<Vec<Var>> = (0..data.num_vertices)
            .map(|i| (0..data.num_vertices).map(|j| add_binvar!(model, name: &format!("x_{}_{}", i, j)).unwrap()).collect())
            .collect::<Vec<_>>();

        for i in 0..data.num_vertices {
            model.set_obj_attr(UB, &x[i][i], 0.0).unwrap();
        }

        let _u: Vec<Var> = (0..data.num_vertices).map(|i| add_ctsvar!(model, name: &format!("u_{}", i), bounds: 0.0..data.num_vertices).unwrap()).collect();

        let theta_p: Vec<Vec<Var>> = (1..data.num_vertices).map(|i| (0..gamma+1).map(|g| add_ctsvar!(model, name: &format!("theta_p_{}_{}", i, g), bounds: 0.0..data.capacity).unwrap()).collect()).collect();
        let theta_m: Vec<Vec<Var>> = (1..data.num_vertices).map(|i| (0..gamma+1).map(|g| add_ctsvar!(model, name: &format!("theta_p_{}_{}", i, g), bounds: 0.0..data.capacity).unwrap()).collect()).collect();

        let theta_d: Vec<Var> = (1..data.num_vertices).map(|i| add_ctsvar!(model, name: &format!("theta_d_{}", i), bounds: 0.0..data.capacity).unwrap()).collect();


        let obj_expr = (0..data.num_vertices)
            .map(|i| (0..data.num_vertices)
                .map(|j| data.distances[i][j] as f64 * x[i][j])
                .collect::<Vec<Expr>>())
            .flatten()
            .grb_sum();

        model.set_objective(obj_expr, Minimize).unwrap();

        // Flow out constraints
        for i in 1..data.num_vertices {
            let cover_expr = (0..data.num_vertices)
                .map(|j| x[i][j])
                .collect::<Vec<Var>>()
                .grb_sum();
            model.add_constr(&format!("Cover_{}", i), c!(cover_expr == 1.0)).unwrap();
        }

        // Flow in constraints
        for j in 1..data.num_vertices {
            let cover_expr = (0..data.num_vertices)
                .map(|i| x[i][j])
                .collect::<Vec<Var>>()
                .grb_sum();
            model.add_constr(&format!("Cover_{}", j), c!(cover_expr == 1.0)).unwrap();
        }

        // Leaving the depot
        let leaving_depot_expr = (1..data.num_vertices)
            .map(|j| x[0][j])
            .collect::<Vec<Var>>()
            .grb_sum();
        model.add_constr("Leaving_depot", c!(leaving_depot_expr >= 1.0)).unwrap();

        // Depot flow conservation
        let leaving_depot_expr = (1..data.num_vertices)
            .map(|j| x[0][j])
            .collect::<Vec<Var>>()
            .grb_sum();
        let arrival_depot_expr = (1..data.num_vertices)
            .map(|i| x[i][0])
            .collect::<Vec<Var>>()
            .grb_sum();
        model.add_constr("Depot_flow_conservation", c!(leaving_depot_expr - arrival_depot_expr == 0.0)).unwrap();

        // MTZ constraints
        for i in 0..data.num_vertices {
            for j in 1..data.num_vertices {
                if i != j {
                    model.add_constr(&format!("MTZ_{}_{}", i, j), c!(_u[j] >= _u[i] + 1 - ((data.num_vertices) as f64) * (1 - x[i][j]))).unwrap();
                }
            }
        }

        // let subtours = vec![(2, 1), (1, 2)];

        // for (i, j) in subtours {
        //     model.add_constr(&format!("MTZ_{}_{}", i, j), c!(u[j] >= u[i] + 1 - ((data.num_vertices) as f64) * (1 - x[i][j]))).unwrap();
        // }

        // Negative case loads
        for i in 1..data.num_vertices {
            for j in 1..data.num_vertices {
                for g in 0..gamma+1 {
                    // TODO: Double check d is what i want it to do
                    model.add_constr(&format!("Neg_Norm_{}_{}_{}", i, j, g), c!(theta_m[j-1][g] <= theta_m[i-1][g] + data.d[j] * x[i][j] + (data.capacity as f64) * (1.0 - x[i][j]))).unwrap();
                }
            }
        }

        for j in 1..data.num_vertices {
            model.add_constr(&format!("Fix_theta_d_m_{}", j), c!(theta_m[j-1][0] <= theta_d[j-1] + data.d[j] as f64 * x[0][j] + data.capacity as f64 * (1.0 - x[0][j]))).unwrap();
        }

        // Positive case loads
        for i in 1..data.num_vertices {
            for j in 1..data.num_vertices {
                for g in 0..gamma+1 {
                    // TODO: Double check d is what i want it to do
                    model.add_constr(&format!("Pos_Norm_{}_{}_{}", i, j, g), c!(theta_p[j-1][g] >= theta_p[i-1][g] + data.d[j] * x[i][j] - (data.capacity as f64) * (1.0 - x[i][j]))).unwrap();
                }
            }
        }
        
        for j in 1..data.num_vertices {
            model.add_constr(&format!("Fix_theta_d_p_{}", j), c!(theta_p[j-1][0] >= theta_d[j-1] + data.d[j] as f64 * x[0][j] - (data.capacity as f64) * (1.0 - x[0][j]))).unwrap();
        }


        let pc = if deviation.is_none() {
            0.0
        } else {
            deviation.unwrap()
        };
        for i in 1..data.num_vertices {
            for j in 1..data.num_vertices {
                for g in 1..gamma+1 {
                    // TODO: Double check d is what i want it to do
                    model.add_constr(&format!("Pos_WC_{}_{}_{}", i, j, g), c!(theta_p[j-1][g] >= theta_p[i-1][g-1] + (data.d[j] as f64 + dev(data.d[j], pc)) * x[i][j] - (data.capacity as f64) * (1.0 - x[i][j]))).unwrap();
                }
            }
        }

        for i in 1..data.num_vertices {
            for j in 1..data.num_vertices {
                for g in 1..gamma+1 {
                    // TODO: Double check d is what i want it to do
                    model.add_constr(&format!("Neg_WC_{}_{}_{}", i, j, g), c!(theta_m[j-1][g] <= theta_m[i-1][g-1] + (data.d[j] as f64 - dev(data.d[j], pc)) * x[i][j] + (data.capacity as f64) * (1.0 - x[i][j]))).unwrap();
                }
            }
        }

        for j in 1..data.num_vertices {
            for g in 1..gamma+1 {
                model.add_constr(&format!("Fix_theta_d_p_{}_{}", j, g), c!(theta_p[j-1][g] >= theta_d[j-1] + (data.d[j] as f64 + dev(data.d[j], pc)) * x[0][j] - data.capacity as f64 * (1.0 - x[0][j]))).unwrap();
            }
        }

        for j in 1..data.num_vertices {
            for g in 1..gamma+1 {
                model.add_constr(&format!("Fix_theta_d_m_{}_{}", j, g), c!(theta_m[j-1][g] <= theta_d[j-1] + (data.d[j] as f64 - dev(data.d[j], pc)) * x[0][j] + data.capacity as f64 * (1.0 - x[0][j]))).unwrap();
            }
        }

        model.update().unwrap();

        CompactModel { 
            model,
            gamma,
            deviation,
            data,

            // Variables
            x,
            _u,
            _theta_p: theta_p,
            _theta_m: theta_m,
            theta_d,
            clique_cuts: vec![],
        }
    }

    fn add_valid_inequalities(&mut self) {
        for i in 1..self.data.num_vertices {
            for j in 1..self.data.num_vertices {
                if i != j {
                    let mut invalid = vec![];
                    let demand = self.data.d[i] + self.data.d[j];

                    for k in 1..self.data.num_vertices {
                        if k != i && k != j {
                            if (demand + self.data.d[k]).abs() as usize > self.data.capacity {
                                invalid.push(k);
                            }
                        }
                    }

                    if invalid.len() > 0 {
                        let lhs = self.x[i][j] + invalid.iter().map(|&k| self.x[j][k]).grb_sum();
                        let cut = self.model.add_constr(&format!("Clique_{}_{}_out", i, j), c!(lhs <= 1.0)).unwrap();
                        self.clique_cuts.push(cut);

                        let lhs = self.x[j][i] + invalid.iter().map(|&k| self.x[k][i]).grb_sum();
                        let cut = self.model.add_constr(&format!("Clique_{}_{}_in", i, j), c!(lhs <= 1.0)).unwrap();
                        self.clique_cuts.push(cut);
                    }
                }
            }
        }
        println!("Added {} clique cuts", self.clique_cuts.len());
    }

    pub fn solve(&mut self) {
        self.add_valid_inequalities();
        self.model.optimize().unwrap();
        self.print_solution();
    }

    fn print_solution(&self) {
        let sol = self.model.get_attr(ObjVal).unwrap();
        let mut arcs = vec![];
        for i in 0..self.data.num_vertices {
            for j in 0..self.data.num_vertices {
                if self.model.get_obj_attr(X, &self.x[i][j]).unwrap() > 0.5 {
                    arcs.push((i, j));
                }
            }
        }

        while arcs.len() > 0 {
            let mut curr = 0;
            let mut idx = 0;
            let mut gamma = 0;
            let mut load = 0.0;
            let mut deviation = 0.0;
            let mut is_wc = false;
            while idx < arcs.len() {
                if arcs[idx].0 == curr {
                    let next = arcs[idx].1;
                    if curr == 0 {
                        let var = self.theta_d[next-1];
                        let val = self.model.get_obj_attr(X, &var).unwrap();
                        load = val;
                        print!("0 (Initial load: {}) ->\n", load);
                    }
                    let demand = self.data.d[next];
                    load += demand as f64;
                    if next == 0 {
                        print!("Depot\n\n");
                        deviation = 0.0;
                    } else {
                        if gamma < self.gamma {
                            // let wc_minus_c = self.model.get_constr_by_name(&format!("Neg_WC_{}_{}_{}", curr, arcs[idx].1, gamma)).unwrap().unwrap();
                            let norm_minus_c = if curr == 0 {
                                if gamma == 0 {
                                        self.model.get_constr_by_name(&format!("Fix_theta_d_m_{}", next)).unwrap()
                                    } else {
                                        self.model.get_constr_by_name(&format!("Fix_theta_d_m_{}_{}", next, gamma)).unwrap()
                                    }
                                } else {
                                    self.model.get_constr_by_name(&format!("Neg_Norm_{}_{}_{}", curr, next, gamma)).unwrap()
                                };
                            if norm_minus_c.is_none() {
                                println!("Constraint not found");
                                println!("{:?}", format!("Neg_Norm_{}_{}_{}", curr, next, gamma));
                            }
                            let slack_norm = self.model.get_obj_attr(Slack, &norm_minus_c.unwrap()).unwrap();
                            if slack_norm > 1e-6 {
                                gamma += 1;
                                is_wc = true;
                                deviation += dev(self.data.d[curr], self.deviation.unwrap());
                            }
                        }
                        
                        print!("\tWC: {}\n\tDemand: {}\n{} (Load: {} - {}) ->\n", is_wc, demand, next, load - deviation, load + deviation);
                    }
                    is_wc = false;
                    curr = next;
                    arcs.swap_remove(idx);
                    idx = 0;
                }
                else {
                    idx += 1;
                }
            }
        }
        println!("\n\nObjective value: {}\n", sol);
    }
}

fn dev(d: isize, p: f64) -> f64 {
    ((d as f64).abs() * p).ceil()
}

pub struct BranchCutModel {
    model: Model,
    gamma: usize,
    deviation: Option<f64>,
    data: RBRPData,

    // Variables
    x: Vec<Vec<Var>>,
    
    _rounded_capacity_cuts: Vec<Constr>,
    clique_cuts: Vec<Constr>,
    _infeasibility_cuts: Vec<Constr>,
}

impl BranchCutModel {
    pub fn new(filename: &str, gamma: usize, deviation: Option<f64>) -> Self {
        // Read data from file
        let data = read_from_file(filename);
        let mut model = Model::new("RBRP Compact Model").unwrap();

        // Create variables
        let x: Vec<Vec<Var>> = (0..data.num_vertices)
            .map(|i| (0..data.num_vertices).map(|j| add_binvar!(model, name: &format!("x_{}_{}", i, j)).unwrap()).collect())
            .collect::<Vec<_>>();

        for i in 0..data.num_vertices {
            model.set_obj_attr(UB, &x[i][i], 0.0).unwrap();
        }

        let obj_expr = (0..data.num_vertices)
            .map(|i| (0..data.num_vertices)
                .map(|j| data.distances[i][j] as f64 * x[i][j])
                .collect::<Vec<Expr>>())
            .flatten()
            .grb_sum();

        model.set_objective(obj_expr, Minimize).unwrap();

        // Flow out constraints
        for i in 1..data.num_vertices {
            let cover_expr = (0..data.num_vertices)
                .map(|j| x[i][j])
                .collect::<Vec<Var>>()
                .grb_sum();
            model.add_constr(&format!("Cover_{}", i), c!(cover_expr == 1.0)).unwrap();
        }

        // Flow in constraints
        for j in 1..data.num_vertices {
            let cover_expr = (0..data.num_vertices)
                .map(|i| x[i][j])
                .collect::<Vec<Var>>()
                .grb_sum();
            model.add_constr(&format!("Cover_{}", j), c!(cover_expr == 1.0)).unwrap();
        }

        // Leaving the depot
        let leaving_depot_expr = (1..data.num_vertices)
            .map(|j| x[0][j])
            .collect::<Vec<Var>>()
            .grb_sum();
        model.add_constr("Leaving_depot", c!(leaving_depot_expr >= 1.0)).unwrap();

        // Depot flow conservation
        let leaving_depot_expr = (1..data.num_vertices)
            .map(|j| x[0][j])
            .collect::<Vec<Var>>()
            .grb_sum();
        let arrival_depot_expr = (1..data.num_vertices)
            .map(|i| x[i][0])
            .collect::<Vec<Var>>()
            .grb_sum();
        model.add_constr("Depot_flow_conservation", c!(leaving_depot_expr - arrival_depot_expr == 0.0)).unwrap();

        model.update().unwrap();

        let mut bcutmodel = BranchCutModel { 
            model,
            gamma,
            deviation,
            data,

            // Variables
            x,
            
            _rounded_capacity_cuts: vec![],
            clique_cuts: vec![],
            _infeasibility_cuts: vec![],
        };

        bcutmodel.add_valid_inequalities();

        bcutmodel
    }

    fn add_valid_inequalities(&mut self) {
        for i in 1..self.data.num_vertices {
            for j in 1..self.data.num_vertices {
                if i != j {
                    let mut invalid = vec![];
                    let demand = self.data.d[i] + self.data.d[j];

                    for k in 1..self.data.num_vertices {
                        if k != i && k != j {
                            if (demand + self.data.d[k]).abs() as usize > self.data.capacity {
                                invalid.push(k);
                            }
                        }
                    }

                    if invalid.len() > 0 {
                        let lhs = self.x[i][j] + invalid.iter().map(|&k| self.x[j][k]).grb_sum();
                        let cut = self.model.add_constr(&format!("Clique_{}_{}_out", i, j), c!(lhs <= 1.0)).unwrap();
                        self.clique_cuts.push(cut);

                        let lhs = self.x[j][i] + invalid.iter().map(|&k| self.x[k][i]).grb_sum();
                        let cut = self.model.add_constr(&format!("Clique_{}_{}_in", i, j), c!(lhs <= 1.0)).unwrap();
                        self.clique_cuts.push(cut);
                    }
                }
            }
        }
        println!("Added {} clique cuts", self.clique_cuts.len());
    }

    pub fn solve(&mut self) {
        let mut cb = SeperationCallback {
            var_data: self.x.iter().enumerate().flat_map(|(i, row)| {
                row.iter().enumerate().map(|(j, &x_var)| VarData { x: x_var, i, j }).collect::<Vec<_>>()
            }).collect(),
            d: self.data.d.clone(),
            p: self.deviation,
            gamma: self.gamma,
            capacity: self.data.capacity,
            num_vertices: self.data.num_vertices,
        };
        self.model.set_param(LazyConstraints, 1).unwrap();
        self.model.optimize_with_callback(&mut cb).unwrap();
        self.print_solution();
    }

    fn print_solution(&self) {
        let sol = self.model.get_attr(ObjVal).unwrap();
        let mut arcs = vec![];
        for i in 0..self.data.num_vertices {
            for j in 0..self.data.num_vertices {
                if self.model.get_obj_attr(X, &self.x[i][j]).unwrap() > 0.5 {
                    arcs.push((i, j));
                }
            }
        }
        let mut routes = vec![];
        let mut route = vec![];
        let mut curr = 0;
        while arcs.len() > 0 {
            let mut to_remove: isize = -1;
            for (idx, arc) in arcs.iter().enumerate() {
                if arc.0 == curr {
                    route.push(curr);
                    print!(" {} ->", curr);
                    let next = arc.1;
                    if next == 0 {
                        println!(" Depot");
                        routes.push(route);
                        route = vec![];
                    }
                    curr = next;
                    to_remove = idx as isize;
                    break;
                }
            }
            arcs.remove(to_remove as usize);
        }
        for route in routes.iter() {
            println!("Route: {:?}", route);
            self.verify_route_is_robust(route.clone());
        }
        println!("\n\nObjective value: {}\n", sol);
    }

    fn verify_route_is_robust(&self, route: Vec<usize>) {
        let mut theta_m = vec![vec![0.0; self.gamma+1]; route.len()];
        let mut theta_p = vec![vec![0.0; self.gamma+1]; route.len()];

        for j in 0..route.len() {
            for g in 0..self.gamma + 1 {
                if j == 0 {
                    theta_m[j][g] = 0.0;
                }
                else if g == 0 {
                    theta_m[j][g] = theta_m[j-1][g] + self.data.d[route[j]] as f64;
                }
                else {
                    theta_m[j][g] = (theta_m[j-1][g] + self.data.d[route[j]] as f64)
                        .min(theta_m[j-1][g] + self.data.d[route[j]] as f64 - dev(self.data.d[route[j]], self.deviation.unwrap_or(0.0)));
                }
            }
        }

        let mut theta_d = -(theta_m.iter().flatten().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap());
        if theta_d < 0.0 {
            theta_d = 0.0;
        }

        for j in 0..route.len() {
            for g in 0..self.gamma + 1 {
                if 0.0 <= theta_m[j][g] + theta_d && theta_m[j][g] + theta_d <= self.data.capacity as f64 {
                    // Valid
                } else {
                    println!("Route {:?} is not robust at index {} with gamma {}", route, j, g);
                    println!("theta_m[{}][{}] = {}", j, g, theta_m[j][g]);
                    println!("theta_d = {}", theta_d);
                    return;
                }
            }
        }

        for j in 0..route.len() {
            for g in 0..self.gamma + 1 {
                if j == 0 {
                    theta_p[j][g] = theta_d;
                }
                else if g == 0 {
                    theta_p[j][g] = theta_p[j-1][g] + self.data.d[route[j]] as f64;
                }
                else {
                    theta_p[j][g] = (theta_p[j-1][g] + self.data.d[route[j]] as f64)
                        .min(theta_p[j-1][g] + self.data.d[route[j]] as f64 - dev(self.data.d[route[j]], self.deviation.unwrap_or(0.0)));
                }
            }
        }

        for j in 0..route.len() {
            for g in 0..self.gamma + 1 {
                if 0.0 <= theta_p[j][g] && theta_p[j][g] <= self.data.capacity as f64 {
                    // Valid
                } else {
                    println!("Route {:?} is not robust at index {} with gamma {}", route, j, g);
                    println!("theta_p[{}][{}] = {}", j, g, theta_p[j][g]);
                    println!("theta_d = {}", theta_d);
                    return;
                }
            }
        }
    }
}

#[derive(Debug, Clone)]
struct VarData {
    x: Var,
    i: usize,
    j: usize,
}

struct SeperationCallback {
    var_data: Vec<VarData>,
    d: Vec<isize>,
    p: Option<f64>,
    gamma: usize,
    capacity: usize,
    num_vertices: usize,
}

impl SeperationCallback {
    fn seperate_routes(&self, ctx: &MIPSolCtx) -> Vec<Vec<&VarData>> {
        let vars = self.var_data.iter().map(|v| v.x);
        let sol = ctx.get_solution(vars).unwrap();
        let mut arcs: Vec<&VarData> = self.var_data.iter().enumerate()
            .filter(|(idx, _)| sol[*idx] > 0.5)
            .map(|(_, vd)| vd)
            .collect::<Vec<_>>();
        arcs.sort_by(|a, b| a.i.cmp(&b.i));
        let mut routes = vec![];
        let mut curr_idx = 0;
        let mut start = arcs[curr_idx].i;
        let mut route = vec![];

        loop {
            if arcs.len() == 0 {
                assert!(route.len() == 0, "Semi route found when no arcs left");
                break;
            }
            let vd = arcs.remove(curr_idx);
            if vd.j == start {
                route.push(vd);
                // println!("Route found: {:?}", route.iter().map(|vd| (vd.i, vd.j)).collect::<Vec<_>>());
                routes.push(route);
                route = vec![];
                curr_idx = arcs.iter().position(|vd_nxt| vd_nxt.i == 0).unwrap_or(0);
                if arcs.len() == 0 {
                    break;
                }
                start = arcs[curr_idx].i;
            }
            else {
                for (idx, vd_nxt) in arcs.iter().enumerate() {
                    if vd_nxt.i == vd.j {
                        route.push(vd);
                        curr_idx  = idx;
                        break;
                    }
                }
            }
        }
        routes
    }

    fn add_capacity_cut(&self, ctx: &MIPSolCtx, in_nodes: &Vec<usize>, wc_nodes: &Vec<usize>) {
        let mut rhs: usize = 1;
        rhs = rhs.max(
            ((in_nodes.iter().map(|&i| self.d[i] as f64).sum::<f64>()
            + wc_nodes.iter().map(|&i| dev(self.d[i], self.p.unwrap_or(0.0))).sum::<f64>())
            / (self.capacity as f64)).abs().ceil() as usize);
        rhs = rhs.max(
            ((in_nodes.iter().map(|&i| self.d[i] as f64).sum::<f64>()
            - wc_nodes.iter().map(|&i| dev(self.d[i], self.p.unwrap_or(0.0))).sum::<f64>())
            / (self.capacity as f64)).abs().ceil() as usize);

        let mut vars: Vec<Var> = vec![];
        let out_nodes: Vec<usize> = (0..self.num_vertices).filter(|i| !in_nodes.contains(i)).collect();
        
        for out_node in out_nodes.iter() {
            for in_node in in_nodes.iter() {
                let var_data = &self.var_data[out_node * self.num_vertices + in_node];
                assert!(var_data.j == *in_node && var_data.i == *out_node, "VarData indices do not match");
                vars.push(var_data.x);
            }
        }

        let lhs = vars.iter().grb_sum();
        ctx.add_lazy(c!(lhs >= rhs as f64)).unwrap();
    }

    fn add_tournament_cut(&self, ctx: &MIPSolCtx, route: Vec<&VarData>) {
        let k = route.len();

        let lhs = route[0].x +
            (1..k-1).map(|i| {
                (i+1..k).map(|j| {
                    let var_data = &self.var_data[route[i].i * self.num_vertices + route[j].j];
                    assert!(var_data.i == route[i].i && var_data.j == route[j].j, "VarData indices do not match");
                    var_data.x
                }).grb_sum()
            }).grb_sum();

        let rhs = k-1;
        ctx.add_lazy(c!(lhs <= rhs as f64)).unwrap();
        println!("Added tournament cut for route {:?}", route.iter().map(|vd| vd.i).collect::<Vec<_>>());
    }   
}

impl Callback for SeperationCallback {
    fn callback(&mut self, w: Where) -> CbResult {
        match w {
            Where::MIPSol(ctx) => {
                let routes = self.seperate_routes(&ctx);
                let mut added_cuts = false;
                for route in routes.iter() {
                    if !route.iter().any(|vd| vd.i == 0) {
                        // Sub-tour cuts
                        // let lhs = route.iter().map(|vd| vd.x).grb_sum();
                        // let mut rhs = 1.0;
                        // let tot_demand = route.iter().map(|vd| self.d[vd.i] as f64).sum::<f64>();

                        // if self.p.is_some() {
                        //     let mut devs: Vec<f64> = route.iter().map(|vd| dev(self.d[vd.i], self.p.unwrap())).collect();
                        //     devs.sort_by(|a, b| b.partial_cmp(a).unwrap());
                        //     let mut dev_sum = 0.0;
                        //     for i in 0..self.gamma {
                        //         if i < devs.len() {
                        //             dev_sum += devs[i];
                        //         } else {
                        //             break;
                        //         }
                        //     }

                        //     let up = ((tot_demand + dev_sum) / (self.capacity as f64)).abs().ceil();
                        //     let down = ((tot_demand - dev_sum) / (self.capacity as f64)).abs().ceil();

                        //     if up > rhs {
                        //         rhs = up;
                        //     }
                        //     if down > rhs {
                        //         rhs = down;
                        //     }
                        // } else {
                        //     if (tot_demand / (self.capacity as f64)).abs() > rhs {
                        //         rhs = (tot_demand / (self.capacity as f64)).abs().ceil();
                        //     }
                        // }
                        // ctx.add_lazy(c!(lhs >= rhs)).unwrap();
                        let nodes = route.iter().map(|vd| vd.i).collect::<Vec<_>>();
                        self.add_capacity_cut(&ctx, &nodes, &vec![]);
                        // println!("Added sub-tour cut for nodes {:?}", nodes);
                        added_cuts = true;
                    }
                }
                if added_cuts == false {
                    for route in routes.iter() {
                        let mut theta_minus: f64 = 0.0;
                        let mut theta_plus: f64 = 0.0;
                        let mut theta_minus_min: f64 = 0.0;
                        let mut theta_plus_max: f64 = 0.0;
                        let mut l = vec![];
                        let mut first_infeasible: isize = -1;
                        let mut sep_cap_added_cuts = false;

                        for k in 1..route.len() {
                            let i = route[k].i;
                            theta_minus += self.d[i] as f64;
                            theta_plus += self.d[i] as f64;

                            if self.gamma > 0 {
                                if l.len() < self.gamma {
                                    theta_minus -= dev(self.d[i], self.p.unwrap());
                                    theta_plus += dev(self.d[i], self.p.unwrap());
                                    l.push(i);
                                }
                                else {
                                    let j = l.iter().enumerate().fold(0, |min_idx, (idx, &wc)| {
                                        if dev(self.d[wc], self.p.unwrap()) < dev(self.d[l[min_idx]], self.p.unwrap()) {
                                            idx
                                        } else {
                                            min_idx
                                        }
                                    });
                                    let dev_i = dev(self.d[i], self.p.unwrap());
                                    let dev_j = dev(self.d[l[j]], self.p.unwrap());
                                    if dev_i > dev_j {
                                        theta_minus += dev_j - dev_i;
                                        theta_plus += dev_i - dev_j;
                                        l[j] = i;
                                    }
                                }
                            }
                            theta_minus_min = theta_minus_min.min(theta_minus);
                            theta_plus_max = theta_plus_max.max(theta_plus);

                            if theta_plus_max - theta_minus_min > self.capacity as f64 {
                                if first_infeasible == -1 {
                                    // println!("Tournament required");
                                    first_infeasible = k as isize;
                                }
                            }

                            if theta_minus.abs() > self.capacity as f64 || theta_plus.abs() > self.capacity as f64 {
                                // Add capacity cut
                                let nodes = route[1..k+1].iter().map(|vd| vd.i).collect();
                                self.add_capacity_cut(&ctx, &nodes, &l);
                                // println!("Added capacity cut for nodes {:?}", nodes);
                                sep_cap_added_cuts = true;
                            }
                        }
                        if sep_cap_added_cuts == false && first_infeasible > 0 {
                            // add tournament cut
                            self.add_tournament_cut(&ctx, route[..(first_infeasible as usize + 1)].iter().map(|vd| *vd).collect());
                            // println!("Added tournament cut");
                        }
                    }

                }
            },
            _ => {},
        }
        Ok(())
    }
}