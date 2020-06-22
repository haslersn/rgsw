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

fn transpose_impl<T>(
    majors: usize,
    minors: usize,
    src: &[T],
    dest: &mut [T],
    my_majors_begin: usize,
    my_majors_end: usize,
    my_minors_begin: usize,
    my_minors_end: usize,
) where
    T: Clone,
{
    let my_majors = my_majors_end - my_majors_begin;
    let my_minors = my_minors_end - my_minors_begin;
    if my_majors * my_minors <= 256 {
        for maj in my_majors_begin..my_majors_end {
            for (d, s) in dest
                .iter_mut()
                .skip(my_minors_begin * majors + maj)
                .step_by(majors)
                .take(my_minors)
                .zip(
                    src.iter()
                        .skip(maj * minors + my_minors_begin)
                        .take(my_minors),
                )
            {
                *d = s.clone();
            }
        }
    } else if my_majors >= my_minors {
        transpose_impl(
            majors,
            minors,
            src,
            dest,
            my_majors_begin,
            my_majors_begin + my_majors / 2,
            my_minors_begin,
            my_minors_end,
        );
        transpose_impl(
            majors,
            minors,
            src,
            dest,
            my_majors_begin + my_majors / 2,
            my_majors_end,
            my_minors_begin,
            my_minors_end,
        );
    } else {
        transpose_impl(
            majors,
            minors,
            src,
            dest,
            my_majors_begin,
            my_majors_end,
            my_minors_begin,
            my_minors_begin + my_minors / 2,
        );
        transpose_impl(
            majors,
            minors,
            src,
            dest,
            my_majors_begin,
            my_majors_end,
            my_minors_begin + my_minors / 2,
            my_minors_end,
        );
    }
}

fn transpose<T>(majors: usize, minors: usize, src: &[T], dest: &mut [T])
where
    T: Clone + PartialEq,
{
    assert!(src.len() == majors * minors);
    assert!(src.len() == dest.len());
    transpose_impl(majors, minors, src, dest, 0, majors, 0, minors);
}

// DFT for prime-power index
fn dft(p: usize, power: u32, data: &mut [Residue], extra_buffer: &mut [Residue]) {
    assert!(power >= 1);
    let m_ = p.pow(power - 1);
    let m = m_ * p;
    assert!(data.len() == m);
    // T_m \cdot (DFT_p \otimes I_{[m']})
    for j1 in 0..m_ {
        for (j0_out, output) in extra_buffer.iter_mut().skip(j1).step_by(m_).enumerate() {
            *output = Residue::zero();
            // DFT_p \otimes I_{[m']}
            for (j0_in, input) in data.iter().skip(j1).step_by(m_).enumerate() {
                // TODO: Precompute the matrix with (i,j)th entry INDEX_TH_ROOT.pow(i*j).
                *output += &(*input * INDEX_TH_ROOT.pow((j0_out * j0_in) as u32));
            }
            // T_m
            if j1 > 0 && j0_out > 0 {
                *output *= &INDEX_TH_ROOT.pow((j0_out * j1) as u32);
            }
        }
    }
    // I_{[p]} \otimes DFT_{m'}
    if power > 1 {
        for start in (0..m).step_by(m_) {
            let end = start + m_;
            // Swap extra_buffer and data.
            dft(
                p,
                power - 1,
                &mut extra_buffer[start..end],
                &mut data[start..end],
            );
        }
    }
    transpose(p, m_, extra_buffer, data);
}

// CRT for prime-power index
fn crt(p: usize, power: u32, data: &mut [Residue], extra_buffer: &mut [Residue]) {
    assert!(power >= 1);
    let m_ = p.pow(power - 1);
    let totient = m_ * (p - 1);
    assert!(data.len() == totient);
    // \hat T_m \cdot (CRT_p \otimes I_{[m']})
    for j1 in 0..m_ {
        for (j0_out, output) in extra_buffer.iter_mut().skip(j1).step_by(m_).enumerate() {
            *output = Residue::zero();
            // CRT_p \otimes I_{[m']}
            for (j0_in, input) in data.iter().skip(j1).step_by(m_).enumerate() {
                // TODO: Precompute the matrix with (i,j)th entry INDEX_TH_ROOT.pow(i*j).
                *output += &(*input * INDEX_TH_ROOT.pow(((j0_out + 1) * j0_in) as u32));
            }
            // \hat T_m
            if j1 > 0 {
                *output *= &INDEX_TH_ROOT.pow(((j0_out + 1) * j1) as u32);
            }
        }
    }
    // I_{\mathbb Z_p^*} \otimes DFT_{m'}
    if power > 1 {
        for start in (0..totient).step_by(m_) {
            let end = start + m_;
            // Swap extra_buffer and data.
            dft(
                p,
                power - 1,
                &mut extra_buffer[start..end],
                &mut data[start..end],
            );
        }
    }
    transpose(p - 1, m_, extra_buffer, data);
}

// Inverse DFT for prime-power index
fn inv_dft(p: usize, power: u32, data: &mut [Residue], extra_buffer: &mut [Residue]) {
    assert!(power >= 1);
    let m_ = p.pow(power - 1);
    let m = m_ * p;
    assert!(data.len() == m);
    transpose(m_, p, data, extra_buffer);
    // I_{[p]} \otimes DFT_{m'}^{-1}
    if power > 1 {
        for start in (0..m).step_by(m_) {
            let end = start + m_;
            // Swap extra_buffer and data.
            inv_dft(
                p,
                power - 1,
                &mut extra_buffer[start..end],
                &mut data[start..end],
            );
        }
    }
    // (DFT_p^{-1} \otimes I_{[m']}) \cdot T_m^{-1}
    for j1 in 0..m_ {
        // T_m^{-1}
        if j1 > 0 {
            for (j0_in, input) in extra_buffer
                .iter_mut()
                .skip(j1)
                .step_by(m_)
                .enumerate()
                .skip(1)
            {
                *input *= &INV_INDEX_TH_ROOT.pow((j0_in * j1) as u32);
            }
        }
        // DFT_p^{-1} \otimes I_{[m']}
        for (j0_out, output) in data.iter_mut().skip(j1).step_by(m_).enumerate() {
            *output = Residue::zero();
            for (j0_in, input) in extra_buffer.iter().skip(j1).step_by(m_).enumerate() {
                // TODO: This is a HACK(!) for p = 2
                if j0_out + j0_in == 1 {
                    *output -= input;
                } else if j0_out + j0_in == 0 {
                    *output += &(*input * INDEX_TH_ROOT);
                } else {
                    *output += input;
                }
            }
            *output *= &INV_INDEX_TH_ROOT_MINUS_ONE;
        }
    }
}

// Inverse CRT for prime-power index
fn inv_crt(p: usize, power: u32, data: &mut [Residue], extra_buffer: &mut [Residue]) {
    assert!(power >= 1);
    let m_ = p.pow(power - 1);
    let totient = m_ * (p - 1);
    assert!(data.len() == totient);
    transpose(m_, p - 1, data, extra_buffer);
    // I_{\mathbb Z_p^*} \otimes DFT_{m'}^{-1}
    if power > 1 {
        for start in (0..totient).step_by(m_) {
            let end = start + m_;
            // Swap extra_buffer and data.
            inv_dft(
                p,
                power - 1,
                &mut extra_buffer[start..end],
                &mut data[start..end],
            );
        }
    }
    // (CRT_p^{-1} \otimes I_{[m']}) \cdot \hat T_m^{-1}
    for j1 in 0..m_ {
        // \hat T_m^{-1}
        if j1 > 0 {
            for (j0_in, input) in extra_buffer.iter_mut().skip(j1).step_by(m_).enumerate() {
                *input *= &INV_INDEX_TH_ROOT.pow(((j0_in + 1) * j1) as u32);
            }
        }
        // CRT_p^{-1} \otimes I_{[m']}
        for (j0_out, output) in data.iter_mut().skip(j1).step_by(m_).enumerate() {
            *output = Residue::zero();
            for (j0_in, input) in extra_buffer.iter().skip(j1).step_by(m_).enumerate() {
                // TODO: This is a HACK(!) for p = 2
                *output += input;
            }
        }
    }
}

impl From<PowerPoly> for ChremPoly {
    fn from(mut other: PowerPoly) -> ChremPoly {
        let mut result = ChremPoly(other.0.clone());
        crt(INDEX_BASE, INDEX_POWER, &mut result.0, &mut other.0);
        result
    }
}

impl From<ChremPoly> for PowerPoly {
    fn from(mut other: ChremPoly) -> PowerPoly {
        let mut result = PowerPoly(other.0.clone());
        inv_crt(INDEX_BASE, INDEX_POWER, &mut result.0, &mut other.0);
        result
    }
}
