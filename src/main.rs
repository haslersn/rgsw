use std::ops::{Add, Mul, Sub};

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

const DEGREE: usize = 24320;
const MODULUS: U384 = U384([
    0xFFFFFFFFFFEF8001,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
]); // 2^(384) − 1081343

#[derive(Clone, Debug)]
struct Residue(U384);

impl Add for Residue {
    type Output = Residue;

    fn add(self, other: Residue) -> Residue {
        let mut a = [0; 7];
        let mut b = [0; 7];
        let mut q = [0; 7];
        a[..6].clone_from_slice(self.0.as_ref());
        b[..6].clone_from_slice(other.0.as_ref());
        q[..6].clone_from_slice(MODULUS.as_ref());
        let mut r = [0; 6];
        r.clone_from_slice(&((U448(a) + U448(b)) % U448(q)).as_ref()[..6]);
        Residue(U384(r))
    }
}

impl Sub for Residue {
    type Output = Residue;

    fn sub(self, other: Residue) -> Residue {
        let mut a = [0; 7];
        let mut b = [0; 7];
        let mut q = [0; 7];
        a[..6].clone_from_slice(self.0.as_ref());
        b[..6].clone_from_slice(other.0.as_ref());
        q[..6].clone_from_slice(MODULUS.as_ref());
        let mut r = [0; 6];
        // First add the modulus q, so that intermediate values are positive.
        r.clone_from_slice(&((U448(q) + U448(a) - U448(b)) % U448(q)).as_ref()[..6]);
        Residue(U384(r))
    }
}

impl Mul for Residue {
    type Output = Residue;

    fn mul(self, other: Residue) -> Residue {
        let mut a = [0; 12];
        let mut b = [0; 12];
        let mut q = [0; 12];
        a[..6].clone_from_slice(self.0.as_ref());
        b[..6].clone_from_slice(other.0.as_ref());
        q[..6].clone_from_slice(MODULUS.as_ref());
        let mut r = [0; 6];
        r.clone_from_slice(&((U768(a) * U768(b)) % U768(q)).as_ref()[..6]);
        Residue(U384(r))
    }
}

struct RingElem {
    coeffs: [Residue; DEGREE],
}

struct CrtElem {
    values: [Residue; DEGREE],
}

fn main() {
    let one = Residue(U384::one());
    let big = one.clone() + one;
    let big = big.clone() * big;
    let big = big.clone() * big;
    let big = big.clone() * big;
    let big = big.clone() * big;
    let big = big.clone() * big;
    let big = big.clone() * big;
    let big = big.clone() * big;
    let big = big.clone() * big;
    let big = big.clone() * big; // 2^512 mod q
    println!("Hello, world! {:?}", big);
}
