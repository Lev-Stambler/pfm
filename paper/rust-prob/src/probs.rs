use num::integer::binomial;
use probability::distribution::{self, Discrete};

use crate::{utils::binomial_coeff, U};

//
type BinomParams = (U, f64);

pub struct Construction {
    pub n: U,
    pub m: U,
    pub deg_stab: U,
    pub deg_bit: U,
    pub p: f64,
}

impl Construction {
    fn f_cannot_decrease_syndrome_given_g(&self, g: U) -> f64 {
        todo!()
    }

    fn f_cannot_decr_syndrome_and_n_errors_given_g(&self, n: U, g: U) -> f64 {
        if n == 0 {
            todo!()
        } else {
            todo!()
        }
    }

    fn f_no_error(&self, g: U) -> f64 {
        let mut p_no_error = 0.;
        for zone_i in 1..(g + 1) {
            let p_distinct_zone = self.f_zones(zone_i, g);
            println!("Pr zone: {}", p_distinct_zone);
            let zone_size = (g * self.deg_stab - 1) / zone_i;
            p_no_error += self.p_error_in_zone(zone_size).powi(zone_i as i32) * p_distinct_zone;
        }
        p_no_error
    }

    fn f_zones(&self, n_zones: U, g: U) -> f64 {
        // TODO: rounding
        let zone_size = g * (self.deg_stab - 1) / n_zones;
        let zone_size_plus_1 = g * (self.deg_stab - 1) / (n_zones + 1);

        let p_at_least_zone_n = self.f_overlap_n_zones(n_zones, zone_size, 0);

        let p_at_least_zone_n_plus_1 = self.f_overlap_n_zones(n_zones + 1, zone_size_plus_1, 0);

        let p = p_at_least_zone_n - p_at_least_zone_n_plus_1;
        assert!(p >= 0.);
        p
    }

    fn overlap_params(&self, g: U) -> BinomParams {
        (g * (self.deg_stab - 1), 1.0 / self.n as f64)
    }

    fn lying_stabilizer_prob(&self) -> f64 {
        0.5 + 0.5 * (1. - 2. * self.p).powi(self.deg_stab as i32)
            - (1. - self.p).powi(self.deg_stab as i32)
    }

    fn F_overlap_stab_neighb(&self, n_stabs: U, degree: Option<U>, x: U) -> f64 {
        todo!()
    }

    // An upper bound as we consider the Pr that there is an odd or even # of errors
    fn p_error_in_zone(&self, zone_size: U) -> f64 {
        let distr = distribution::Binomial::new(zone_size as usize, self.p);
        1. - distr.mass(0)
    }

    fn f_overlap_n_zones(&self, n_zones: U, zone_size: U, x: U) -> f64 {
        // For each pair of zones, there are zone_size many events where there can be an overlap
        // This seems to be an upper bound (i.e. the probabilities will be too high for a larger x) as
        // zone size / N should decrease with every successful event
        // i.e. we are sampling w/o replacement...
        // TODO: update to reflect this?

        let distr = distribution::Binomial::new(
            binomial(n_zones, 2) as usize * zone_size as usize,
            zone_size as f64 / self.n as f64,
        );
        distr.mass(x as usize)
    }
}
