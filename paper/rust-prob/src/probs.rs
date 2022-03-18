use factorial::Factorial;

use num::{integer::binomial, BigUint, FromPrimitive, ToPrimitive};
use probability::distribution::{self, Discrete};

use crate::{utils::binomial_coeff, U};

const DEBUG: bool = false;

type BinomParams = (U, f64);

pub struct Construction {
    pub n: U,
    pub m: U,
    pub deg_stab: U,
    pub deg_bit: U,
    pub p: f64,
}

impl Construction {
    pub fn default() -> Self {
        Self {
            n: 10_000_000,
            m: 5_000_000,
            deg_stab: 11,
            deg_bit: 6,
            p: 0.01,
        }
    }
}

impl Construction {
    pub fn f_cannot_decrease_syndrome_given_g(&self, g: U) -> f64 {
        let mut p = 0.;
        // In all honesty its not additive because these events are not independent. Its prob like Pr[n errors and lying | not n-1 errors and not lying and not n-2 errors and not lying ....]
        // But an upper bound should be fine...
        for i in 0..self.deg_stab + 1 {
            let t = self.f_cannot_decr_syndrome_and_n_errors_given_g(i, g);
            // println!("for {} errors the prob is {}", i, t);
            p += t;
        }
        p.min(1.)
    }

    fn f_cannot_decr_syndrome_and_n_errors_given_g(&self, n: U, g: U) -> f64 {
        if n == 0 {
            self.f_no_error(g)
        } else {
            return 0.;
            let binom_mult = binomial(self.deg_stab, n);
            let mut s = 0.0;
            let min_lying_per_bit = ((self.deg_bit as f64) / 2.).ceil().to_u128().unwrap();
            let lying_n_lower_bound = min_lying_per_bit * n;
            for l in lying_n_lower_bound..n * self.deg_bit {
                // TODO: this may overflow...
                // See the following
                // Basically you can think of having l balls and n bins.
                // Start with putting n * ceil(degbit/2) balls so each bin has that many
                // then, stars and bars for n-1 bars with (1 + l - lying_n_lower_bound) positions
                // https://math.stackexchange.com/questions/794716/the-number-of-different-ways-to-choose-k-out-of-n-unique-items
                // TODO: verify
                let numb_distribution_of_lying = binomial(
                    BigUint::from_u128(1 + l - lying_n_lower_bound + n - 1 - 1).unwrap(),
                    BigUint::from_u128(n - 1).unwrap(),
                );
                // println!("Max {}, Given l {}, n {}, num_dist {}", self.deg_bit * n, l, n, numb_distribution_of_lying);
                let g_to_check = (g + l).checked_sub(n * self.deg_bit);
                if let Some(g_prime) = g_to_check {
                    let pr_all_lying = self.lying_stabilizer_prob_given_1_error().powi(l as i32);
                    // println!("pr all lying {}", pr_all_lying * numb_distribution_of_lying);
                    s += (numb_distribution_of_lying.to_f64().unwrap() * pr_all_lying)
                        * self.f_no_error(g_prime);
                }
            }
            let ret = s * binom_mult as f64;
            ret
        }
    }

    // This seems wrong... its like multiplicative me thinks. Like p_no_error
    pub fn f_no_error(&self, g: U) -> f64 {
        let mut p_no_error = 0.;
        for zone_i in 1..(g + 1) {
            let p_distinct_zone = self.f_zones(zone_i, g);
            if p_distinct_zone > 0. {
                let zone_size = (g * (self.deg_stab - 1)) as f64 / zone_i as f64;
                if DEBUG {
                    println!(
                        "g {}, zone_i: {}: self.p_error_in_zone(zone_size: {}) {} {}",
                        g,
                        zone_i,
                        zone_size,
                        self.p_error_in_zone(zone_size),
                        p_distinct_zone
                    );
                }
                let t = self.p_error_in_zone(zone_size).powi(zone_i as i32) * p_distinct_zone;
                // TODO: this is wrong... (I think)
                p_no_error += t;
            }
        }
        p_no_error
    }

    pub fn f_zones(&self, n_zones: U, g: U) -> f64 {
        // if n_zones
        let zone_size = g as f64 * (self.deg_stab - 1) as f64 / n_zones as f64;

        let n_pairs = binomial(n_zones, 2);
        let p_at_least_zone_n = self
            .f_overlap_n_zones(n_zones, zone_size.round().to_u128().unwrap(), zone_size, 0)
            .powi(n_pairs as i32);

        let p_at_least_zone_n_plus_1 = if n_zones == g {
            0.
        } else {
            let zone_size_plus_1 = (g as f64 * (self.deg_stab - 1) as f64) / (n_zones + 1) as f64;
            let n_plus_1_pairs = binomial(n_zones + 1, 2);
            self.f_overlap_n_zones(
                n_zones + 1,
                zone_size.round().to_u128().unwrap(),
                zone_size_plus_1,
                0,
            )
            .powi(n_plus_1_pairs as i32)
        };

        // TODO: this is from rounding err?
        if p_at_least_zone_n_plus_1 > p_at_least_zone_n {
            return 0.;
        }
        let p = p_at_least_zone_n - p_at_least_zone_n_plus_1;
        assert!(p >= 0.);
        p
    }

    fn overlap_params(&self, g: U) -> BinomParams {
        (g * (self.deg_stab - 1), 1.0 / self.n as f64)
    }

    // probability neighouring - 1 have odd errors
    fn lying_stabilizer_prob_given_1_error(&self) -> f64 {
        0.5 - 0.5 * (1. - 2. * self.p).powi(self.deg_stab as i32 - 1)
    }

    fn lying_stabilizer_prob(&self) -> f64 {
        0.5 + 0.5 * (1. - 2. * self.p).powi(self.deg_stab as i32)
            - (1. - self.p).powi(self.deg_stab as i32)
    }

    fn F_overlap_stab_neighb(&self, n_stabs: U, degree: Option<U>, x: U) -> f64 {
        todo!()
    }

    // An upper bound as we consider the Pr that there is an odd or even # of errors
    fn p_error_in_zone(&self, zone_size: f64) -> f64 {
        // TODO: rounding
        let distr = distribution::Binomial::new(zone_size.round().to_usize().unwrap(), self.p);
        1. - distr.mass(0)
    }

    fn f_overlap_n_zones(&self, n_zones: U, zone_size: U, zone_size_not_rounded: f64, x: U) -> f64 {
        // For each pair of zones, there are zone_size many events where there can be an overlap
        // This seems to be an upper bound (i.e. the probabilities will be too high for a larger x) as
        // zone size / N should decrease with every successful event
        // i.e. we are sampling w/o replacement...
        // TODO: update to reflect this?

        let distr = distribution::Binomial::new(
            binomial(n_zones, 2) as usize * zone_size as usize,
            zone_size_not_rounded / self.n as f64,
        );
        distr.mass(x as usize)
    }
}
