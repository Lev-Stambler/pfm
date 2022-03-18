use crate::U;
use memoize::memoize;
use num::bigint::{BigInt, ToBigUint};
use num::integer::{binomial, gcd, multinomial};
use num::{BigUint, FromPrimitive, ToPrimitive};
use probability::prelude::*;

const P_FRACTION_DENOM: u64 = 1_000_000_000_000;

/// Calculate r * a / b, avoiding overflows and fractions.
///
/// Assumes that b divides r * a evenly.
fn multiply_and_divide_f(r: f64, a: f64, b: f64) -> f64 {
    // See http://blog.plover.com/math/choose-2.html for the idea.
    let g = gcd(r as u128, b as u128) as f64;
    r / g.clone() * (a / (b / g))
}

fn binomial_coeff(mut n: f64, k: f64) -> f64 {
    // See http://blog.plover.com/math/choose.html for the idea.
    if k > n {
        return 0.;
    }
    if k > n.clone() - k.clone() {
        return binomial_coeff(n.clone(), n - k);
    }
    let mut r = 1.;
    let mut d = 1.;
    loop {
        if d > k {
            break;
        }
        r = multiply_and_divide_f(r, n.clone(), d.clone());
        n = n - 1.;
        d = d + 1.;
    }
    r
}

fn multinomial_coeef(k: &[f64]) -> f64 {
    let mut r = 1.;
    let mut p = 0.;
    for i in k {
        p = p + i;
        r = r * binomial_coeff(p.clone(), i.clone());
    }
    r
}

// TODO: f128
#[memoize]
/// p is represented as p_fraction / P_FRACTION_DENOM
/// if k = 1, then it is the smallest of n_samples, if k = n_samples, then it is the largest
// Following http://www.math.wm.edu/~leemis/2005informsjoc.pdf
pub(crate) fn binomial_order_stat_pdf(
    N: U,
    p_fraction: u64,
    n_samples: U,
    r_order_stat: U,
    val: U,
) -> f64 {
    let distr = Binomial::new(N as usize, p_fraction as f64 / P_FRACTION_DENOM as f64);

    if val == 0 {
        (0..(n_samples + 1 - r_order_stat))
            .map(|w| {
                binomial_coeff(n_samples as f64, w as f64) as f64
                    * distr.mass(0).powi((n_samples - w) as i32)
                    * (1.0 - distr.distribution(1.0)).powi(w as i32)
            })
            .sum::<f64>()
            .min(1.)
    } else if val == N {
        let mut sum = 0.0;
        for u in 0..r_order_stat {
            sum += binomial_coeff(n_samples as f64, u as f64) as f64
                * distr.distribution(N as f64 - 1.).powi(u as i32)
                * distr.mass(N as usize).powi((n_samples - u) as i32);
        }
        sum.min(1.)
    } else {
        let mut sum = 0.0;
        for u in 0..r_order_stat {
            for w in 0..(n_samples - r_order_stat + 1) {
                // TODO: hmmmmmm... f64 implementation of multinomial and binomial
                sum += multinomial(&vec![
                    BigUint::from_u128(n_samples).unwrap(),
                    BigUint::from_u128(n_samples - u - w).unwrap(),
                    BigUint::from_u128(w).unwrap(),
                ])
                .to_f64().unwrap()
                    * distr.distribution((val - 1) as f64).powi(u as i32)
                    * distr.mass(val as usize).powi((n_samples - u - w) as i32)
                    * (1. - distr.distribution((val + 1) as f64)).powi(w as i32);
            }
        }
        sum
    }

    // println!(
    //     "distwithin: {}, distr below: {}, distr above: {}, b: {}",
    //     distr.mass(val as usize),
    //     distr
    //         .distribution(val as f64)
    //         .powi((k_th_order_stat - 1) as i32),
    //     (1.0 - distr.distribution(val as f64)).powi((n_samples - k_th_order_stat) as i32),
    //     binomial(n_samples as u128, k_th_order_stat as u128 - 1)
    // );
    // (n_samples - k_th_order_stat + 1) as f64 // out of the RVs not choosen to be below the order stat, choose the order stat
    //     * distr.mass(val as usize)
    //     * binomial(n_samples as u128, (k_th_order_stat - 1) as u128) as u128 as f64 // There are n_samples choose k-1 ways to select the variables which must be lower than
    //     // the k_th order stat. Then, the remainder must be above
    //     * distr
    //         .distribution(val as f64)
    //         .powi((k_th_order_stat - 1) as i32)
    //     * (1.0 - distr.distribution(val as f64)).powi((n_samples - k_th_order_stat) as i32)
    // // TODO: this is wrong cause it gives above w/o equality. need equality
}

#[cfg(test)]
mod test {
    use crate::utils::{binomial_order_stat_pdf, P_FRACTION_DENOM};

    #[test]
    fn test_order_stat_pdf() {
        let n = 8;
        let p_fraction = (0.5 * P_FRACTION_DENOM as f64) as u64;
        println!(
            "AAAA {}",
            binomial_order_stat_pdf(n, p_fraction, 100, 50, 4)
        )
    }
}
