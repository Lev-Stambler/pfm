use crate::U;
use memoize::memoize;
use num_integer::binomial;
use probability::prelude::*;

const P_FRACTION_DENOM: u64 = 1_000_000_000_000;

#[memoize]
/// p is represented as p_fraction / P_FRACTION_DENOM
/// if k = 1, then it is the smallest of n_samples, if k = n_samples, then it is the largest
pub(crate) fn binomial_order_stat_pdf(
    n: U,
    p_fraction: u64,
    n_samples: U,
    k_th_order_stat: U,
    val: U,
) -> f64 {
    let distr = Binomial::new(n as usize, p_fraction as f64 / P_FRACTION_DENOM as f64);
    n_samples as f64
        * distr.mass(val as usize)
        * binomial(n_samples, k_th_order_stat - 1) as f64 // There are n_samples choose k-1 ways to select the variables which must be lower than
        // the k_th order stat. Then, the remainder must be above
        * distr
            .distribution(val as f64)
            .powi((k_th_order_stat - 1) as i32)
        * (1.0 - distr.distribution(val as f64).powi((n_samples - k_th_order_stat) as i32))
}

#[cfg(test)]
mod test {
    use crate::utils::{binomial_order_stat_pdf, P_FRACTION_DENOM};

    #[test]
    fn test_order_stat_pdf() {
        let n = 100;
        let p_fraction = (0.5 * P_FRACTION_DENOM as f64) as u64;
        println!("{}", binomial_order_stat_pdf(n, p_fraction, 100, 100, 70))
    }
}