use crate::probs::Construction;

mod probs;
mod utils;
pub type U = u128;

fn main() {
    let c = Construction::default();
    for i in 0..c.deg_stab*c.deg_bit {
    // for i in 49..51 {
        // Defo should not be 0.000000478088549214902 for g = 1
        println!("{}: {}", i, c.f_cannot_decrease_syndrome_given_g(i));
    }
    // let mut s = 0.0;
    // for i in 1..100 {
    //     let t = c.f_zones(i, 100);
    //     s += t;
    //     println!("{} {}", s, t);
    // }
}
