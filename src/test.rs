use crate::residue::*;
use rand::thread_rng;
use rand::Rng;
use std::iter::{once, repeat};

fn example_residues<'a, R: Rng>(rng: &'a mut R) -> impl Iterator<Item = Residue> + 'a {
    once(Residue::zero())
        .chain(once(Residue::one()))
        .chain(once(INDEX_TH_ROOT))
        .chain(once(INV_INDEX_TH_ROOT))
        .chain(once(INV_INDEX_TH_ROOT_MINUS_ONE))
        .chain(repeat(()).map(move |_| rng.gen()))
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
