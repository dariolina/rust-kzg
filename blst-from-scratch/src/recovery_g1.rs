extern crate alloc;

use alloc::string::String;
use alloc::vec;
use alloc::vec::Vec;

use kzg::{FFTFr, FFTG1, Fr, G1, G1Mul, Poly, PolyG1, ZeroPoly, PolyRecoverG1};
use crate::consts::G1_IDENTITY;

use crate::types::fft_settings::FsFFTSettings;
use crate::types::fr::FsFr;
use crate::types::g1::FsG1;
use crate::types::poly::{FsPolyG1};


impl PolyRecoverG1<FsG1, FsFr, FsPolyG1, FsFFTSettings> for FsPolyG1 {
    fn recover_poly_from_samples_g1(
    samples: &[Option<FsG1>],
    fs: &FsFFTSettings,
) -> Result<Self, String> {
    let len_samples = samples.len();
    if !len_samples.is_power_of_two() {
        return Err(String::from("len_samples must be a power of two"));
    }

    let mut missing: Vec<usize> = Vec::new();
    for (i, sample) in samples.iter().enumerate() {
        if sample.is_none() {
            missing.push(i);
        }
    }

    // Calculate `Z_r,I`
    let (zero_eval, mut zero_poly) = fs
        .zero_poly_via_multiplication(len_samples, &missing)
        .unwrap();

    for i in 0..len_samples {
        if samples[i].is_none() != zero_eval[i].is_zero() {
            return Err(String::from(
                "recovery error: samples should be none when and only when zero_eval is zero",
            ));
        }
    }

    let mut poly_evaluations_with_zero = FsPolyG1::default();

    // Construct E * Z_r,I: the loop makes the evaluation polynomial
    for i in 0..len_samples {
        if samples[i].is_none() {
            poly_evaluations_with_zero.coeffs.push(G1_IDENTITY);
        } else {
            poly_evaluations_with_zero
                .coeffs
                .push(samples[i].unwrap().mul(&zero_eval[i]));
        }
    }
    // Now inverse FFT so that poly_with_zero is (E * Z_r,I)(x) = (D * Z_r,I)(x)
    let mut poly_with_zero = FsPolyG1 {
        coeffs: fs.fft_g1(&poly_evaluations_with_zero.coeffs, true).unwrap(),
    };

    // x -> k * x
    let len_zero_poly = zero_poly.coeffs.len();

    let _zero_poly_scale = (len_zero_poly - 1).next_power_of_two();

    poly_with_zero.scale();
    zero_poly.scale();

    // Q1 = (D * Z_r,I)(k * x)
    let scaled_poly_with_zero = poly_with_zero.coeffs;

    // Q2 = Z_r,I(k * x)
    let scaled_zero_poly = zero_poly.coeffs;

    #[allow(unused_assignments)]
    let mut eval_scaled_poly_with_zero = vec![];
    #[allow(unused_assignments)]
    let mut eval_scaled_zero_poly = vec![];

    // Polynomial division by convolution: Q3 = Q1 / Q2
    eval_scaled_poly_with_zero = fs.fft_g1(&scaled_poly_with_zero, false).unwrap();
    eval_scaled_zero_poly = fs.fft_fr(&scaled_zero_poly, false).unwrap();

    let mut eval_scaled_reconstructed_poly = FsPolyG1 {
        coeffs: eval_scaled_poly_with_zero.clone(),
    };
    for i in 0..len_samples {
        eval_scaled_reconstructed_poly.coeffs[i] = eval_scaled_poly_with_zero[i]
            .mul(&eval_scaled_zero_poly[i].eucl_inverse());
    }

    // The result of the division is D(k * x):
    let scaled_reconstructed_poly: Vec<FsG1> = fs
        .fft_g1(&eval_scaled_reconstructed_poly.coeffs, true)
        .unwrap();

    // k * x -> x
    let mut reconstructed_poly = FsPolyG1{coeffs:scaled_reconstructed_poly};

    // Finally we have D(x) which evaluates to our original data at the powers of roots of unity
    reconstructed_poly.unscale();

    // The evaluation polynomial for D(x) is the reconstructed data:
    let reconstructed_data = fs.fft_g1(&reconstructed_poly.coeffs, false).unwrap();

    // Check all is well
    for i in 0..len_samples {
        if !(samples[i].is_none() || reconstructed_data[i].equals(&samples[i].unwrap())) {
            return Err(String::from(
                "recovery error: samples should be none or equal reconstructed data",
            ));
        }
    }

    Ok(FsPolyG1{coeffs:reconstructed_data})
}
}
