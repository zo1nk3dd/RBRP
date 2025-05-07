use grb::{attribute::{ModelDoubleAttr::ObjVal, VarDoubleAttr::{UB, X}}, prelude::*};

use crate::utils::*;

pub struct CompactModel {
    model: Model,
    gamma: usize,
    deviation: Option<f64>,
    data: RBRPData,

    // Variables
    x: Vec<Vec<Var>>,
    u: Vec<Var>,
    theta_p: Vec<Vec<Var>>,
    theta_m: Vec<Vec<Var>>,
    theta_d: Vec<Var>,
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

        let u: Vec<Var> = (0..data.num_vertices).map(|i| add_ctsvar!(model, name: &format!("u_{}", i), bounds: 0.0..data.num_vertices).unwrap()).collect();

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
                    model.add_constr(&format!("MTZ_{}_{}", i, j), c!(u[j] >= u[i] + 1 - ((data.num_vertices) as f64) * (1 - x[i][j]))).unwrap();
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
                model.add_constr(&format!("Fix_theta_d_p_{}_{}", j, g), c!(theta_m[j-1][g] <= theta_d[j-1] + (data.d[j] as f64 - dev(data.d[j], pc)) * x[0][j] + data.capacity as f64 * (1.0 - x[0][j]))).unwrap();
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
            u,
            theta_p,
            theta_m,
            theta_d,
        }
    }

    pub fn solve(&mut self) {
        self.model.optimize().unwrap();
        self.print_solution();
    }

    fn print_solution(&self) {
        let sol = self.model.get_attr(ObjVal).unwrap();
        println!("\n\nObjective value: {}\n", sol);
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
            while idx < arcs.len() {
                if arcs[idx].0 == curr {
                    print!(" {} ({}-{}) -> ", curr);
                    curr = arcs[idx].1;
                    arcs.swap_remove(idx);
                    idx = 0;
                    if curr == 0 {
                        print!(" 0\n");
                    }
                }
                else {
                    idx += 1;
                }
            }
        }
    }
}

fn dev(d: isize, p: f64) -> f64 {
    ((d as f64).abs() * p).ceil()
}