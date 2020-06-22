use std::convert::{TryFrom, TryInto};
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

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

pub const DEGREE: usize = 16384;
pub const INDEX_BASE: usize = 2;
pub const INDEX_POWER: u32 = 15;
pub const MODULUS: U384 = U384([
    0xFFFFFFFFFFEF8001,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
]); // 2^(384) − 1081343
pub const INDEX_TH_ROOT: Residue = Residue(U384([
    0xC65EA8DC3E81B1F8,
    0x187343D949A162D5,
    0xF74AC8D88DE2DC56,
    0x4225431217768CF2,
    0x9C79D95C1EF48D17,
    0x0001C062055D83EA,
])); // 1053046320810374670051386479886365060104392590871858199730685648701498648418052922079079421152377431379245380088
     // which is a primitive 32768th root of unity mod MODULUS.

#[derive(Clone, Copy, Debug)]
pub struct Residue(pub U384);

impl Residue {
    pub fn zero() -> Residue {
        Residue(U384::zero())
    }

    pub fn one() -> Residue {
        Residue(U384::one())
    }

    pub fn from_u64(x: u64) -> Residue {
        Residue(U384([x, 0, 0, 0, 0, 0]))
    }

    pub fn from_i64(x: i64) -> Residue {
        if x >= 0 {
            Residue(U384([x as u64, 0, 0, 0, 0, 0]))
        } else {
            -Residue(U384([-x as u64, 0, 0, 0, 0, 0]))
        }
    }
}

impl Neg for Residue {
    type Output = Residue;

    fn neg(mut self) -> Residue {
        self.0 = MODULUS - self.0;
        self
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
