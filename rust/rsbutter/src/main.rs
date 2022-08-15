use std::f64::consts::PI;
const N_SOS_COEFFS: u8 = 6;

// linear spaced vector for spacing poles.
fn linspace(start: i32, stop: i32, step: i32) {
    m = ((stop - start)/step).round()
    let y = Vec::with_capacity(m)
    for d in 0..m-1 {
        y.push(start + d*step)
    }
    y
}

// simple in-line from euler's identity.
fn euler(x: f64) -> num::complex::Complex {
    num::complex::Complex::new(x.cos(), x.sin())
}

fn main() {
    println!("Hello, world!");
}
