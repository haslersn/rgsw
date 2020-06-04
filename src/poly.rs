use crate::residue::*;
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

trait Poly {
    fn eval(self, x: &Residue) -> Residue;
}

pub struct ChremPoly(pub [Residue; DEGREE]);

impl AddAssign<&ChremPoly> for ChremPoly {
    fn add_assign(&mut self, other: &ChremPoly) {
        for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
            *a += b;
        }
    }
}

impl Add for ChremPoly {
    type Output = ChremPoly;

    fn add(mut self, other: ChremPoly) -> ChremPoly {
        self += &other;
        self
    }
}

impl SubAssign<&ChremPoly> for ChremPoly {
    fn sub_assign(&mut self, other: &ChremPoly) {
        for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
            *a -= b;
        }
    }
}

impl Sub for ChremPoly {
    type Output = ChremPoly;

    fn sub(mut self, other: ChremPoly) -> ChremPoly {
        self -= &other;
        self
    }
}

impl MulAssign<&ChremPoly> for ChremPoly {
    fn mul_assign(&mut self, other: &ChremPoly) {
        for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
            *a *= b;
        }
    }
}

impl Mul for ChremPoly {
    type Output = ChremPoly;

    fn mul(mut self, other: ChremPoly) -> ChremPoly {
        self *= &other;
        self
    }
}

impl Poly for ChremPoly {
    fn eval(self, x: &Residue) -> Residue {
        Residue::zero() // TODO: Implement this!
    }
}

pub struct CoeffPoly(pub [Residue; DEGREE]);

impl AddAssign<&CoeffPoly> for CoeffPoly {
    fn add_assign(&mut self, other: &CoeffPoly) {
        for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
            *a += b;
        }
    }
}

impl Add for CoeffPoly {
    type Output = CoeffPoly;

    fn add(mut self, other: CoeffPoly) -> CoeffPoly {
        self += &other;
        self
    }
}

impl SubAssign<&CoeffPoly> for CoeffPoly {
    fn sub_assign(&mut self, other: &CoeffPoly) {
        for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
            *a -= b;
        }
    }
}

impl Sub for CoeffPoly {
    type Output = CoeffPoly;

    fn sub(mut self, other: CoeffPoly) -> CoeffPoly {
        self -= &other;
        self
    }
}

impl MulAssign<&CoeffPoly> for CoeffPoly {
    fn mul_assign(&mut self, other: &CoeffPoly) {
        // TODO: Implement this!
    }
}

impl Mul for CoeffPoly {
    type Output = CoeffPoly;

    fn mul(mut self, other: CoeffPoly) -> CoeffPoly {
        self *= &other;
        self
    }
}

impl Poly for CoeffPoly {
    fn eval(self, x: &Residue) -> Residue {
        let mut result = Residue::zero();
        for coeff in self.0.iter().rev() {
            result *= x;
            result += coeff
        }
        result
    }
}
