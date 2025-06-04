# RBRP

In order to run this program, you will need to install rust. Follow the installation guide [here](https://doc.rust-lang.org/book/ch01-01-installation.html).

Once you have installed rust, the program is run using the following command.
```
cargo run -- <instance> <gamma> <deviation>
```
The instance should just be the integer at the beginning of the name (1-65). Results are reported for gamma values of 1, 5 and 10, but any value could be used. Results are reported for deviations of 0.1, 0.25 and 0.5, but any float between 0 and 1 could be used. The program might complain if you ask for the deterministic case (gamma = 0) and provide a non-zero deviation, or vice versa. 

The main.rs file contains some config information for running the program. If you wish, altering the following block (lines 187-190) in the main function will change the behaviour of the algorithm. 

```
let tlimit = TLIMIT;
let frac_st = true;
let frac_cap = true;
let robust_capacity_cuts = true;
```

`frac_st` is whether the fractional subtour cuts are generated. `frac_cap` is whether the fractional rounded capacity cuts are generated. `robust_capacity_cuts` is whether the new idea of removing all permutations of infeasible routes (when appropriate) is used. 

If you wish to use the compact formulation, change the `BranchCutModel` to `CompactModel` on line 201.