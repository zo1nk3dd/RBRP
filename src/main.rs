use rbrp::model::*;
use rbrp::utils::Config;
use std::{env, fs::File, io::Write};
use chrono::Local;

use grb::{attribute::ModelDoubleAttr::{MIPGap, ObjBound, ObjVal, Runtime}, Status};

const FILENAMES: [&str; 65]  = [
    "instances/01-Bari30_RBRP.txt",
    "instances/02-Bari20_RBRP.txt",
    "instances/03-Bari10_RBRP.txt",
    "instances/04-ReggioEmilia30_RBRP.txt",
    "instances/05-ReggioEmilia20_RBRP.txt",
    "instances/06-ReggioEmilia10_RBRP.txt",
    "instances/07-Bergamo30_RBRP.txt",
    "instances/08-Bergamo20_RBRP.txt",
    "instances/09-Bergamo12_RBRP.txt",
    "instances/10-Parma30_RBRP.txt",
    "instances/11-Parma20_RBRP.txt",
    "instances/12-Parma10_RBRP.txt",
    "instances/13-Treviso30_RBRP.txt",
    "instances/14-Treviso20_RBRP.txt",
    "instances/15-Treviso10_RBRP.txt",
    "instances/16-LaSpezia30_RBRP.txt",
    "instances/17-LaSpezia20_RBRP.txt",
    "instances/18-LaSpezia10_RBRP.txt",
    "instances/19-BuenosAires30_RBRP.txt",
    "instances/20-BuenosAires20_RBRP.txt",
    "instances/21-Ottawa30_RBRP.txt",
    "instances/22-Ottawa20_RBRP.txt",
    "instances/23-Ottawa10_RBRP.txt",
    "instances/24-SanAntonio30_RBRP.txt",
    "instances/25-SanAntonio20_RBRP.txt",
    "instances/26-SanAntonio10_RBRP.txt",
    "instances/27-Brescia30_RBRP.txt",
    "instances/28-Brescia20_RBRP.txt",
    "instances/29-Brescia11_RBRP.txt",
    "instances/30-Roma30_RBRP.txt",
    "instances/31-Roma20_RBRP.txt",
    "instances/32-Roma18_RBRP.txt",
    "instances/33-Madison30_RBRP.txt",
    "instances/34-Madison20_RBRP.txt",
    "instances/35-Madison10_RBRP.txt",
    "instances/36-Guadalajara30_RBRP.txt",
    "instances/37-Guadalajara20_RBRP.txt",
    "instances/38-Guadalajara11_RBRP.txt",
    "instances/39-Dublin30_RBRP.txt",
    "instances/40-Dublin20_RBRP.txt",
    "instances/41-Dublin11_RBRP.txt",
    "instances/42-Denver30_RBRP.txt",
    "instances/43-Denver20_RBRP.txt",
    "instances/44-Denver10_RBRP.txt",
    "instances/45-RioDeJaneiro30_RBRP.txt",
    "instances/46-RioDeJaneiro20_RBRP.txt",
    "instances/47-RioDeJaneiro10_RBRP.txt",
    "instances/48-Boston30_RBRP.txt",
    "instances/49-Boston20_RBRP.txt",
    "instances/50-Boston16_RBRP.txt",
    "instances/51-Torino30_RBRP.txt",
    "instances/52-Torino20_RBRP.txt",
    "instances/53-Torino10_RBRP.txt",
    "instances/54-Toronto30_RBRP.txt",
    "instances/55-Toronto20_RBRP.txt",
    "instances/56-Toronto12_RBRP.txt",
    "instances/57-Miami30_RBRP.txt",
    "instances/58-Miami20_RBRP.txt",
    "instances/59-Miami10_RBRP.txt",
    "instances/60-CiudadDeMexico_30_RBRP.txt",
    "instances/61-CiudadDeMexico_20_RBRP.txt",
    "instances/62-CiudadDeMexico_17_RBRP.txt",
    "instances/63-Minneapolis30_RBRP.txt",
    "instances/64-Minneapolis20_RBRP.txt",
    "instances/65-Minneapolis10_RBRP.txt",
];

const TLIMIT: usize = 900;

fn benchmark() {
    let instances= vec![1, 3, 5, 8, 11, 14, 19];
    let gammas = vec![0, 1, 5, 10];
    let devs = vec![0.0, 0.1, 0.25, 0.5];
    let tlimit = TLIMIT;

    let folder = "Results";
    let datetime = Local::now().format("Final").to_string();
    std::fs::create_dir_all(format!("{}/{}", folder, datetime)).expect("Failed to create directory");

    let frac_st = false;
    let frac_cap = false;
    let robust_capacity_cuts = false;
    let name = "Compact";

    std::fs::create_dir_all(format!("{}/{}/{}", folder, datetime, name)).expect("Failed to create directory");

    for instance in instances.iter() {
        println!("Benchmarking instance {}:", instance);
        let mut writer = File::create(format!("{}/{}/{}/{}.txt", folder, datetime, name, instance)).expect("Failed to create file");
        for gamma in gammas.iter() {
            for dev in devs.iter() {
                if (*dev == 0.0 && *gamma > 0) || (*dev > 0.0 && *gamma == 0) {
                    continue;
                }
                let config = Config {
                    gamma: *gamma,
                    dev: *dev,
                    tlimit,
                    frac_st,
                    frac_cap,
                    robust_capacity_cuts,
                };

                let mut model = CompactModel::new(FILENAMES[instance-1], config);
                model.solve();

                let status = model.model.status()
                    .expect("Failed to get model status");

                match status {
                    Status::Optimal => {
                        let obj_val = model.model.get_attr(ObjVal)
                            .expect("Failed to get objective value");
                        let time = model.model.get_attr(Runtime).expect("Failed to get runtime");
                        let buf = format!(
                            "OPTIMAL, Gamma: {}, Dev: {}, LB: {:.2}, Gap: 0, Time: {:.2}s\n",
                            gamma, dev, obj_val, time
                        );
                        writer.write_all(buf.as_bytes()).expect("Failed to write to file");
                        },
                    Status::Infeasible => {
                        let buf = format!(
                            "INFEASIBLE, Gamma: {}, Dev: {}, LB: -, Gap: 0, Time: 0s\n",
                            gamma, dev,
                        );
                        writer.write_all(buf.as_bytes()).expect("Failed to write to file");
                    },
                    Status::TimeLimit => {
                        let lb = model.model.get_attr(ObjBound)
                            .expect("Failed to get lower bound");
                        let gap = model.model.get_attr(MIPGap)
                            .expect("Failed to get MIP gap");
                        let buf = format!(
                            "TIMELIMIT, Gamma: {}, Dev: {}, LB: {:.2}, Gap: {:.2}, Time: {}s\n",
                            gamma, dev, lb, gap, TLIMIT
                        );
                        writer.write_all(buf.as_bytes()).expect("Failed to write to file");
                    },
                    _ => println!("Instance {} with gamma {} and dev {} has status: {:?}", instance, gamma, dev, status),
                }
                writer.flush().expect("Failed to flush writer");
            }
        }
    }
}

fn main() {
    let msg = "cargo run <instance_number> <gamma> <deviation>  OR  cargo run benchmark";
    let args: Vec<String> = env::args().collect();
    if args.len() > 1 && args[1] == "benchmark" {
        benchmark();
        return;
    }
    let instance = if args.len() > 1 {
        args[1].parse::<usize>().unwrap()
    } else {
        println!("Usage: {}", msg);
        panic!();
    };
    assert!(instance > 0 && instance <= 66, "Instance number must be between 1 and 66.");

    let gamma = if args.len() > 2 {
        args[2].parse::<usize>().unwrap()
    } else {
        println!("Usage: {}", msg);
        panic!();
    };
    
    let dev = if args.len() > 3 {
        args[3].parse::<f64>().unwrap()
    } else {
        0.0
    };

    if dev == 0.0 && gamma > 0 {
        panic!("Deviation is 0.0 and budget is non-zero.");
    }

    let tlimit = TLIMIT;
    let frac_st = false;
    let frac_cap = false;
    let robust_capacity_cuts = true;

    let config = Config {
        gamma,
        dev,
        tlimit,
        frac_st,
        frac_cap,
        robust_capacity_cuts,
    };

    let mut model = BranchCutModel::new(FILENAMES[instance-1], config);
    model.solve();
}
