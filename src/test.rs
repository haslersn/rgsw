use crate::poly::*;
use crate::residue::*;
use rand::thread_rng;
use rand::Rng;
use std::iter::{once, repeat};

fn random_residues<'a, R: Rng>(rng: &'a mut R) -> impl Iterator<Item = Residue> + 'a {
    repeat(()).map(move |_| rng.gen())
}

fn example_residues<'a, R: Rng>(rng: &'a mut R) -> impl Iterator<Item = Residue> + 'a {
    once(Residue::zero())
        .chain(once(Residue::one()))
        .chain(once(INDEX_TH_ROOT))
        .chain(once(INV_INDEX_TH_ROOT))
        .chain(once(INV_INDEX_TH_ROOT_MINUS_ONE))
        .chain(random_residues(rng))
}

#[test]
fn test_residue_neutral_add() {
    for elem in example_residues(&mut thread_rng()).take(100) {
        assert_eq!(elem + Residue::zero(), elem);
    }
}

#[test]
fn test_residue_neutral_mult() {
    for elem in example_residues(&mut thread_rng()).take(100) {
        assert_eq!(elem * Residue::one(), elem);
    }
}

#[test]
fn test_residue_commutative_add() {
    for a in example_residues(&mut thread_rng()).take(10) {
        for b in example_residues(&mut thread_rng()).take(10) {
            assert_eq!(a + b, b + a);
        }
    }
}

#[test]
fn test_residue_commutative_mult() {
    for a in example_residues(&mut thread_rng()).take(10) {
        for b in example_residues(&mut thread_rng()).take(10) {
            assert_eq!(a * b, b * a);
        }
    }
}

#[test]
fn test_index_th_root_order() {
    assert_eq!(INDEX_TH_ROOT.pow(DEGREE as u32), -Residue::one());
    assert_eq!(INDEX_TH_ROOT.pow(DEGREE as u32 * 2), Residue::one());
}

#[test]
fn test_inv_index_th_root() {
    assert_eq!(INV_INDEX_TH_ROOT * INDEX_TH_ROOT, Residue::one());
}

#[test]
fn test_inv_index_th_root_minus_one() {
    assert_eq!(
        INV_INDEX_TH_ROOT_MINUS_ONE * (INDEX_TH_ROOT - Residue::one()),
        Residue::one()
    );
}

#[test]
fn test_crt_roundtrip() {
    let mut rng = thread_rng();
    for _ in 0..3 {
        let mut p1 = PowerPoly(Box::new([Residue::zero(); DEGREE]));
        for (d, s) in p1.0.iter_mut().zip(random_residues(&mut rng)) {
            *d = s;
        }
        let c = ChremPoly::from(p1.clone());
        let p2 = PowerPoly::from(c);
        assert_eq!(p1, p2);
    }
}

#[test]
fn test_inv_crt_roundtrip() {
    let mut rng = thread_rng();
    for _ in 0..3 {
        let mut c1 = ChremPoly(Box::new([Residue::zero(); DEGREE]));
        for (d, s) in c1.0.iter_mut().zip(random_residues(&mut rng)) {
            *d = s;
        }
        let p = PowerPoly::from(c1.clone());
        let c2 = ChremPoly::from(p);
        assert_eq!(c1, c2);
    }
}
