use rbrp::model::*;
use std::{env};

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
fn main() {
    let args: Vec<String> = env::args().collect();
    let instance = if args.len() > 1 {
        args[1].parse::<usize>().unwrap()
    } else {
        panic!("Please provide an instance number as a command line argument.");
    };
    assert!(instance > 0 && instance <= 66, "Instance number must be between 1 and 66.");

    let gamma = if args.len() > 2 {
        args[2].parse::<usize>().unwrap()
    } else {
        panic!("Please provide a gamma value as a command line argument.");
    };
    
    let dev = if args.len() > 3 {
        Some(args[3].parse::<f64>().unwrap())
    } else {
        None
    };

    if dev.is_none() && gamma > 0 {
        panic!("Deviation is none and budget is non-zero.");
    }

    let mut model = BranchCutModel::new(FILENAMES[instance-1], gamma, dev);
    model.solve();
}
