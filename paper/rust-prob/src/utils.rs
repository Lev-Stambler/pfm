use memoize::memoize;
use crate::U;


#[memoize]
/// p is represented as p_fraction / u128::max
pub(crate) fn binomial_order_stat_pdf(n: U, p_fraction: u64, k_th_order_stat: U, val: U) -> f64 {
    todo!()
}
