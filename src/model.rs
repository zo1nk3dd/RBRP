use grb::{attribute::{ConstrDoubleAttr::Slack, ModelDoubleAttr::ObjVal, VarDoubleAttr::{UB, X}}, callback::*, constr::IneqExpr, parameter::{DoubleParam::TimeLimit, IntParam::{BranchDir, LazyConstraints}}, prelude::*};

use crate::utils::*;

pub struct CompactModel {
    pub model: Model,
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
    config: Config,
}

impl CompactModel {
    pub fn new(filename: &str, config: Config) -> Self {
        // Read data from file
        let data = read_from_file(filename);
        let mut model = Model::new("RBRP Compact Model").unwrap();
        let gamma = config.gamma;
        let deviation = if config.dev > 0.0 {
            Some(config.dev)
        } else {
            None
        };

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
            config,
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
        self.model.set_param(TimeLimit, self.config.tlimit as f64).unwrap();
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
    pub model: Model,
    gamma: usize,
    deviation: f64,
    data: RBRPData,
    // Variables
    x: Vec<Vec<Var>>,
    config: Config,
}

impl BranchCutModel {
    pub fn new(filename: &str, config: Config) -> Self {
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

        // let vehicle_obj_expr = (1..data.num_vertices)
        //     .map(|i| {
        //         model.set_obj_attr(BranchPriority, &x[0][i], 1).unwrap();
        //         x[0][i]
        //     })
        //     .collect::<Vec<Var>>()
        //     .grb_sum();

        // model.add_constr("vehicle", c!(vehicle_obj_expr == 6.0)).unwrap();

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
            gamma: config.gamma,
            deviation: config.dev,
            data,
            x,
            config,
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
                            if (demand + self.data.d[k]).abs() as f64 > self.data.capacity as f64 {
                                invalid.push(k);
                            }
                        }
                    }

                    if invalid.len() > 0 {
                        let lhs = self.x[i][j] + invalid.iter().map(|&k| self.x[j][k]).grb_sum();
                        self.model.add_constr(&format!("Clique_{}_{}_out", i, j), c!(lhs <= 1.0)).unwrap();

                        let lhs = self.x[j][i] + invalid.iter().map(|&k| self.x[k][i]).grb_sum();
                        self.model.add_constr(&format!("Clique_{}_{}_in", i, j), c!(lhs <= 1.0)).unwrap();
                    }
                }
            }
        }
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
            stats: CallbackStats {
                subtour: 0,
                capacity: 0,
                tournament: 0,
                frac_st: 0,
                frac_cap: 0,
            },
            config: self.config,
        };
        self.model.set_param(LazyConstraints, 1).unwrap();
        self.model.set_param(BranchDir, 1).unwrap(); // todo: is this better?
        self.model.set_param(TimeLimit, self.config.tlimit as f64).unwrap();
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
                        .min(theta_m[j-1][g-1] + self.data.d[route[j]] as f64 - dev(self.data.d[route[j]], self.deviation));
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
                        .min(theta_p[j-1][g-1] + self.data.d[route[j]] as f64 - dev(self.data.d[route[j]], self.deviation));
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

struct CallbackStats {
    subtour: usize,
    capacity: usize,
    tournament: usize,
    frac_st: usize,
    frac_cap: usize,
}

struct SeperationCallback {
    var_data: Vec<VarData>,
    d: Vec<isize>,
    p: f64,
    gamma: usize,
    capacity: usize,
    num_vertices: usize,
    stats: CallbackStats,
    config: Config,
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

    fn create_capacity_cut(&self, in_nodes: &Vec<usize>, wc_nodes: &Vec<usize>) -> IneqExpr {
        let wc_dev = wc_nodes.iter().map(|&i| dev(self.d[i], self.p)).sum::<f64>();
        let normal_demand = in_nodes.iter().map(|&i| self.d[i] as f64).sum::<f64>();
        
        let mut rhs: usize = (wc_dev * 2.0 / self.capacity as f64).ceil() as usize; // could be 0 if wc_nodes is empty
        rhs = rhs.max(1);
        rhs = rhs.max(
            ((normal_demand + wc_dev) / (self.capacity as f64)).abs().ceil() as usize);
        rhs = rhs.max(
            ((normal_demand - wc_dev) / (self.capacity as f64)).abs().ceil() as usize);

        // if rhs == 1 {
        //     if 2.0 * wc_dev > self.capacity as f64 {
        //         // println!("Warning: Capacity cut rhs is 1, but wc_dev is too large: {}", wc_dev);
        //         rhs = 2;
        //     }
        // }
        
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
        c!(lhs >= rhs as f64)
    }

    fn add_capacity_cut(&self, ctx: &MIPSolCtx, in_nodes: &Vec<usize>, wc_nodes: &Vec<usize>) {
        ctx.add_lazy(self.create_capacity_cut(in_nodes, wc_nodes)).unwrap();
        // println!("Added capacity cut for nodes {:?}, wc: {:?} with rhs: {:?}", in_nodes, wc_nodes, rhs);
    }

    fn add_tournament_cut(&self, ctx: &MIPSolCtx, route: Vec<usize>) {
        let k = route.len();

        let lhs = self.var_data[route[0] * self.num_vertices + route[1]].x + 
            (1..k-1).map(|i| {
                // println!("");
                (i+1..k).map(|j| {
                    let var_data = &self.var_data[route[i] * self.num_vertices + route[j]];
                    assert!(var_data.i == route[i] && var_data.j == route[j], "VarData indices do not match");
                    // print!("{}, {}: ", var_data.i, var_data.j);
                    var_data.x
                }).grb_sum()
            }).grb_sum();

        // let lhs_perm = (0..k).map(|i| {
        //     (i+1..k).map(|j| {
        //         &self.var_data[route[i] * self.num_vertices + route[j]].x
        //     }).grb_sum()
        // }).grb_sum() + (1..k).map(|i| {
        //     (i+1..k).map(|j| {
        //         &self.var_data[route[j] * self.num_vertices + route[i]].x
        //     }).grb_sum()
        // }).grb_sum();

        let rhs = k-2;
        ctx.add_lazy(c!(lhs <= rhs as f64)).unwrap();
    }   

    fn seperate_frac_sep_st(&self, ctx: &MIPNodeCtx) -> bool {
        let soln = ctx.get_solution(self.var_data.iter().map(|vd| vd.x).collect::<Vec<_>>());

        if soln.is_err() {
            return false;
        }

        let soln = soln.unwrap();

        let mut data = vec![];
        for (idx, vd) in self.var_data.iter().enumerate() {
            if soln[idx] > EPS {
                let from = vd.i;
                let to = vd.j;
                let cost = soln[idx];
                data.push((from, to, cost));
            }
        }

        let mut graph = Graph::new(data, self.num_vertices);
        let mut added_cuts = false;
        for target in 1..self.num_vertices {
            let flow = graph.min_cut(0, target);
            if flow < 1.0 - EPS {
                added_cuts = true;
                ctx.add_lazy(self.create_capacity_cut(&graph.target_cut, &Vec::new())).unwrap();
                return added_cuts;
            }
            graph.reset();
        }
        added_cuts
    }

    fn seperate_frac_sep_cap(&self, ctx: &MIPNodeCtx) -> bool {
        let mut added_cuts = false;
        let soln = ctx.get_solution(self.var_data.iter().map(|vd| vd.x).collect::<Vec<_>>());

        if soln.is_err() {
            return false;
        }

        let soln = soln.unwrap();
        let n = self.num_vertices;

        let mut data = vec![];
        for (idx, vd) in self.var_data.iter().enumerate() {
            if soln[idx] > EPS {
                let from = vd.i;
                let to = vd.j;
                let cost = soln[idx];
                data.push((from, to, cost));
            }
        }

        let mut required_flow = 0.0;
        for station in 1..n {
            if self.d[station] > 0 {
                data.push((n, station, self.d[station] as f64 / self.capacity as f64));
            }
            if self.d[station] < 0 {
                required_flow += -self.d[station] as f64 / self.capacity as f64;
                data.push((station, n+1, -self.d[station] as f64 / self.capacity as f64));
            }
        }

        let mut graph = Graph::new(data, n + 2);
        let min_cut = graph.min_cut(n, n+1);

        if min_cut < required_flow - EPS {
            let in_nodes: Vec<usize> = graph.target_cut.iter().filter_map(|&i| {if i < n {Some(i)} else {None}}).collect();
            let wc_nodes: Vec<usize> = Vec::new();
            if in_nodes.len() != self.num_vertices {
                let cut = self.create_capacity_cut(&in_nodes, &wc_nodes);
                ctx.add_lazy(cut).unwrap();
                added_cuts = true;
                return added_cuts;
            }
        }

        let mut data = vec![];
        
        for (idx, vd) in self.var_data.iter().enumerate() {
            if soln[idx] > EPS {
                let from = vd.i;
                let to = vd.j;
                let cost = soln[idx];
                data.push((from, to, cost));
            }
        }

        let mut required_flow = 0.0;
        for station in 1..n {
            if self.d[station] > 0 {
                required_flow += self.d[station] as f64 / self.capacity as f64;
                data.push((station, n+1, self.d[station] as f64 / self.capacity as f64));
            }
            if self.d[station] < 0 {
                data.push((n, station, -self.d[station] as f64 / self.capacity as f64));
            }
        }

        let mut graph = Graph::new(data, n + 2);
        let min_cut = graph.min_cut(n, n+1);

        if min_cut < required_flow - EPS {
            let in_nodes: Vec<usize> = graph.target_cut.iter().filter_map(|&i| {if i < n {Some(i)} else {None}}).collect();
            let wc_nodes: Vec<usize> = Vec::new();
            if in_nodes.len() != self.num_vertices {
                let cut = self.create_capacity_cut(&in_nodes, &wc_nodes);
                ctx.add_lazy(cut).unwrap();
                added_cuts = true;
                return added_cuts;
            }
        }

        added_cuts
    }
}

impl Callback for SeperationCallback {
    fn callback(&mut self, w: Where) -> CbResult {
        match w {
            Where::MIPSol(ctx) => {
                let routes = self.seperate_routes(&ctx);
                let mut subtour_cuts = 0;
                let mut cap_cuts = 0;
                let mut tour_cuts = 0;

                for route in routes.iter() {
                    if !route.iter().any(|vd| vd.i == 0) {
                        let nodes = route.iter().map(|vd| vd.i).collect::<Vec<_>>();
                        self.add_capacity_cut(&ctx, &nodes, &vec![]);
                        subtour_cuts += 1;
                    }
                }
                if subtour_cuts == 0 {
                    for route in routes.iter() {
                        let mut theta_minus: f64 = 0.0;
                        let mut theta_plus: f64 = 0.0;
                        let mut theta_minus_min: f64 = 0.0;
                        let mut theta_plus_max: f64 = 0.0;
                        let mut l = vec![];
                        let mut first_infeasible: isize = -1;

                        for k in 1..route.len() {
                            let i = route[k].i;
                            theta_minus += self.d[i] as f64;
                            theta_plus += self.d[i] as f64;

                            if self.gamma > 0 {
                                if l.len() < self.gamma {
                                    theta_minus -= dev(self.d[i], self.p);
                                    theta_plus += dev(self.d[i], self.p);
                                    l.push(i);
                                }
                                else {
                                    // Find the argmin of l
                                    let mut argmin = 0;
                                    let mut min = dev(self.d[l[argmin]], self.p);
                                    for (idx, station) in l.iter().enumerate() {
                                        if dev(self.d[*station], self.p) < min {
                                            argmin = idx;
                                            min = dev(self.d[*station], self.p);
                                        }
                                    }
                                    let dev_i = dev(self.d[i], self.p);
                                    let dev_min = dev(self.d[l[argmin]], self.p);
                                    if dev_i > dev_min {
                                        theta_minus += dev_min - dev_i;
                                        theta_plus += dev_i - dev_min;
                                        l[argmin] = i;
                                    }
                                }
                            }
                            theta_minus_min = theta_minus_min.min(theta_minus);
                            theta_plus_max = theta_plus_max.max(theta_plus);

                            if theta_plus_max - theta_minus_min > self.capacity as f64 {
                                if first_infeasible == -1 {
                                    first_infeasible = k as isize;
                                }
                            }

                            if theta_minus.abs() > self.capacity as f64 || theta_plus.abs() > self.capacity as f64 || (if self.config.robust_capacity_cuts { (theta_plus - theta_minus).abs() > self.capacity as f64 } else { false }) {
                                // Add capacity cut
                                let nodes = route[1..k+1].iter().map(|vd| vd.i).collect();
                                self.add_capacity_cut(&ctx, &nodes, &l);
                                cap_cuts += 1;
                                // break 'in_route;
                            }
                        }
                        if cap_cuts == 0 && first_infeasible > 0 {
                            // add tournament cut
                            self.add_tournament_cut(&ctx, route[..(first_infeasible as usize + 1)].iter().map(|vd| vd.i).collect());
                            tour_cuts += 1;
                        }
                    }
                }
                self.stats.subtour += subtour_cuts;
                self.stats.capacity += cap_cuts;
                self.stats.tournament += tour_cuts;
            },
            Where::MIPNode(ctx) => {
                let result = if self.config.frac_st { self.seperate_frac_sep_st(&ctx) } else { false };
                if result {
                    self.stats.frac_st += 1;
                } else {
                    let new_result = if self.config.frac_cap { self.seperate_frac_sep_cap(&ctx) } else { false };
                    self.stats.frac_cap += if new_result { 1 } else { 0 };
                }
            }
            _ => {},
        }
        Ok(())
    }
}