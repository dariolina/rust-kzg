use blst_rust::types::fft_settings::FsFFTSettings;
use blst_rust::types::fr::FsFr;
use blst_rust::types::g1::FsG1;
use criterion::{criterion_group, criterion_main, Criterion};
use kzg_bench::benches::fft::{bench_fft_fr, bench_fft_g1,bench_fft_fr_in_place};

fn bench_fft_fr_(c: &mut Criterion) {
    bench_fft_fr::<FsFr, FsFFTSettings>(c);
}

fn bench_fft_fr_in_place_(c: &mut Criterion) {
    bench_fft_fr_in_place::<FsFr, FsFFTSettings>(c);
}

fn bench_fft_g1_(c: &mut Criterion) {
    bench_fft_g1::<FsFr, FsG1, FsFFTSettings>(c);
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = bench_fft_fr_, bench_fft_fr_in_place_,bench_fft_g1_
}

criterion_main!(benches);
