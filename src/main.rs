use std::convert::{TryFrom, TryInto};
use itertools::izip;
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

#[macro_use]
extern crate uint;

construct_uint! {
    pub struct U384(6); // 6 x 64-bit word
}

construct_uint! {
    pub struct U448(7);
}

construct_uint! {
    pub struct U768(12);
}

impl From<U384> for U448 {
    fn from(other: U384) -> U448 {
        let mut data = [0; 7];
        data[..6].clone_from_slice(other.0.as_ref());
        U448(data)
    }
}

impl TryFrom<U448> for U384 {
    type Error = &'static str;

    fn try_from(other: U448) -> Result<U384, &'static str> {
        let data = other.as_ref();
        if data[6..].iter().all(|&x| x == 0) {
            let mut new_data = [0; 6];
            new_data.clone_from_slice(&data[..6]);
            Ok(U384(new_data))
        } else {
            Err("Overflow during conversion from U448 to U384.")
        }
    }
}

impl From<U384> for U768 {
    fn from(other: U384) -> U768 {
        let mut data = [0; 12];
        data[..6].clone_from_slice(other.0.as_ref());
        U768(data)
    }
}

impl TryFrom<U768> for U384 {
    type Error = &'static str;

    fn try_from(other: U768) -> Result<U384, &'static str> {
        let data = other.as_ref();
        if data[6..].iter().all(|&x| x == 0) {
            let mut new_data = [0; 6];
            new_data.clone_from_slice(&data[..6]);
            Ok(U384(new_data))
        } else {
            Err("Overflow during conversion from U768 to U384.")
        }
    }
}

impl From<U448> for U768 {
    fn from(other: U448) -> Self {
        let mut data = [0; 12];
        data[..7].clone_from_slice(other.0.as_ref());
        Self(data)
    }
}

impl TryFrom<U768> for U448 {
    type Error = &'static str;

    fn try_from(other: U768) -> Result<U448, &'static str> {
        let data = other.as_ref();
        if data[7..].iter().all(|&x| x == 0) {
            let mut new_data = [0; 7];
            new_data.clone_from_slice(&data[..7]);
            Ok(U448(new_data))
        } else {
            Err("Overflow during conversion from U768 to U448.")
        }
    }
}

const DEGREE: usize = 24320;
const MODULUS: U384 = U384([
    0xFFFFFFFFFFEF8001,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
]); // 2^(384) − 1081343

#[derive(Clone, Copy, Debug)]
struct Residue(U384);

impl Residue {
    fn zero() -> Residue {
        Residue(U384::zero())
    }
}

impl AddAssign<&Residue> for Residue {
    fn add_assign(&mut self, other: &Residue) {
        let q = Into::<U448>::into(MODULUS);
        let sum = Into::<U448>::into(self.0) + Into::<U448>::into(other.0);
        self.0 = (sum % q).try_into().unwrap();
    }
}

impl Add for Residue {
    type Output = Residue;

    fn add(mut self, other: Residue) -> Residue {
        self += &other;
        self
    }
}

impl SubAssign<&Residue> for Residue {
    fn sub_assign(&mut self, other: &Residue) {
        let q = Into::<U448>::into(MODULUS);
        let diff = q + Into::<U448>::into(self.0) - Into::<U448>::into(other.0);
        self.0 = (diff % q).try_into().unwrap();
    }
}

impl Sub for Residue {
    type Output = Residue;

    fn sub(mut self, other: Residue) -> Residue {
        self -= &other;
        self
    }
}

impl MulAssign<&Residue> for Residue {
    fn mul_assign(&mut self, other: &Residue) {
        let q = Into::<U768>::into(MODULUS);
        let prod = Into::<U768>::into(self.0) * Into::<U768>::into(other.0);
        self.0 = (prod % q).try_into().unwrap();
    }
}

impl Mul for Residue {
    type Output = Residue;

    fn mul(mut self, other: Residue) -> Residue {
        self *= &other;
        self
    }
}

impl MulAssign<u64> for Residue {
    fn mul_assign(&mut self, other: u64) {
        let q = Into::<U448>::into(MODULUS);
        let prod = Into::<U448>::into(self.0) * other;
        self.0 = (prod % q).try_into().unwrap();
    }
}

impl Mul<u64> for Residue {
    type Output = Residue;

    fn mul(mut self, other: u64) -> Residue {
        self *= other;
        self
    }
}

impl Mul<Residue> for u64 {
    type Output = Residue;

    fn mul(self, other: Residue) -> Residue {
        other * self
    }
}

impl MulAssign<i64> for Residue {
    fn mul_assign(&mut self, other: i64) {
        let q = Into::<U448>::into(MODULUS);
        let prod = if other >= 0 {
            Into::<U448>::into(self.0) * other
        } else {
            (q - Into::<U448>::into(self.0)) * ((-other) as u64)
        };
        self.0 = (prod % q).try_into().unwrap();
    }
}

impl Mul<i64> for Residue {
    type Output = Residue;

    fn mul(mut self, other: i64) -> Residue {
        self *= other;
        self
    }
}

impl Mul<Residue> for i64 {
    type Output = Residue;

    fn mul(self, other: Residue) -> Residue {
        other * self
    }
}

struct ChremPoly([Residue; DEGREE]);

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

struct CoeffPoly([Residue; DEGREE]);

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

fn main() {
    let one = Residue(U384::one());
    let big = 2 as u64 * one;
    let big = big * big;
    let big = big * big;
    let big = big * big;
    let big = big * big;
    let big = big * big;
    let big = big * big;
    let big = big * big;
    let big = big * big;
    let big = big * big; // 2^512 mod q
    println!("Hello, world! {:?}", big);
}
