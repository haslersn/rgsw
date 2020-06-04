use crate::residue::*;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

trait Poly {
    fn eval(self, x: &Residue) -> Residue;
}

pub struct ChremPoly(pub [Residue; DEGREE]);

impl Neg for ChremPoly {
    type Output = ChremPoly;

    fn neg(mut self) -> ChremPoly {
        for a in self.0.iter_mut() {
            *a = -*a;
        }
        self
    }
}

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

pub struct PowerPoly(pub [Residue; DEGREE]);

impl Neg for PowerPoly {
    type Output = PowerPoly;

    fn neg(mut self) -> PowerPoly {
        for a in self.0.iter_mut() {
            *a = -*a;
        }
        self
    }
}

impl AddAssign<&PowerPoly> for PowerPoly {
    fn add_assign(&mut self, other: &PowerPoly) {
        for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
            *a += b;
        }
    }
}

impl Add for PowerPoly {
    type Output = PowerPoly;

    fn add(mut self, other: PowerPoly) -> PowerPoly {
        self += &other;
        self
    }
}

impl SubAssign<&PowerPoly> for PowerPoly {
    fn sub_assign(&mut self, other: &PowerPoly) {
        for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
            *a -= b;
        }
    }
}

impl Sub for PowerPoly {
    type Output = PowerPoly;

    fn sub(mut self, other: PowerPoly) -> PowerPoly {
        self -= &other;
        self
    }
}

impl MulAssign<&PowerPoly> for PowerPoly {
    fn mul_assign(&mut self, other: &PowerPoly) {
        // TODO: Implement this!
    }
}

impl Mul for PowerPoly {
    type Output = PowerPoly;

    fn mul(mut self, other: PowerPoly) -> PowerPoly {
        self *= &other;
        self
    }
}

impl Poly for PowerPoly {
    fn eval(self, x: &Residue) -> Residue {
        let mut result = Residue::zero();
        for coeff in self.0.iter().rev() {
            result *= x;
            result += coeff
        }
        result
    }
}
