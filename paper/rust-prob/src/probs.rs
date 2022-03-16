use crate::U;

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
        todo!()
    }

    fn f_zones(&self, n_zones: U, g: U) -> f64 {
        todo!()
    }

		fn overlap_params(&self, g: U) -> BinomParams {
			todo!()
		}

		fn lying_stabilizer_prob(&self) -> f64 {
		todo!()
		}
}
