 #[path = "./local_tests/local_recovery.rs"]
 pub mod local_recovery;

#[cfg(test)]
mod tests {
    //use kzg_bench::tests::recover::{recover_random, recover_simple};
    // uncomment to use the local tests
    use crate::local_recovery::{recover_random, recover_simple, recover_g1};

    use blst_from_scratch::types::fft_settings::FsFFTSettings;
    use blst_from_scratch::types::fr::FsFr;
    use blst_from_scratch::types::g1::FsG1;
    use blst_from_scratch::types::poly::{FsPoly, FsPolyG1};

    // Shared tests
    #[test]
    fn recover_simple_() {
        recover_simple::<FsFr, FsFFTSettings, FsPoly, FsPoly>()
    }

    #[test]
    fn recover_random_() {
        recover_random::<FsFr, FsFFTSettings, FsPoly, FsPoly>()
    }

    #[test]
    fn recover_g1_() {
        recover_g1::<FsG1, FsFr, FsFFTSettings, FsPolyG1, FsPolyG1>()
    }
}
