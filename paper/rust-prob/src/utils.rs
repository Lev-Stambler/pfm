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

pub(crate) fn binomial_coeff(mut n: f64, k: f64) -> f64 {
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
// hmmmmm... seems like it over predicts
#[memoize]
// Following http://www.math.wm.edu/~leemis/2005informsjoc.pdf
// Approximate with a guassian
// See https://people.sc.fsu.edu/~jburkardt/presentations/truncated_normal.pdf for truncated Guassian estimates
pub(crate) fn binomial_order_stat_pdf(
    N: U,
    p_fraction: u64,
    n_samples: U,
    r_order_stat: U,
    val: U,
) -> f64 {
    let p = p_fraction as f64 / P_FRACTION_DENOM as f64;
    let mu = N as f64 * p;
    let sigma = (N as f64 * p * (1. - p)).sqrt();
    let distr = Gaussian::new(mu, sigma);

    // TODO: on val - 1???
    (r_order_stat as f64
        * binomial_coeff(n_samples as f64, r_order_stat as f64) as f64
        * distr.density(val as f64)
        * distr.distribution(val as f64).powi((r_order_stat - 1) as i32)
        * (1. - distr.distribution(val as f64)).powi((n_samples - r_order_stat) as i32)).min(1.)
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
            binomial_order_stat_pdf(n, p_fraction, 1_000, 500, 4)
        )
    }
}
