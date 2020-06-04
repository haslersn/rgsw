#[macro_use]
extern crate uint;

mod poly;
mod residue;

use crate::residue::*;

fn main() {
    let one = Residue::one();
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
