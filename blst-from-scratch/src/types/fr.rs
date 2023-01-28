use blst::{
    blst_fr, blst_fr_add, blst_fr_cneg, blst_fr_eucl_inverse, blst_fr_from_scalar,
    blst_fr_from_uint64, blst_fr_inverse, blst_fr_mul, blst_fr_sqr, blst_fr_sub, blst_scalar,
    blst_scalar_from_fr, blst_uint64_from_fr, blst_scalar_fr_check, blst_scalar_from_lendian,
};
use kzg::Fr;

#[derive(Debug)]
pub struct FsFr(pub blst::blst_fr);

impl Fr for FsFr {
    fn default() -> Self {
        Self(blst_fr::default())
    }

    fn null() -> Self {
        Self::from_u64_arr(&[u64::MAX, u64::MAX, u64::MAX, u64::MAX])
    }

    fn zero() -> Self {
        Self::from_u64(0)
    }

    fn one() -> Self {
        Self::from_u64(1)
    }

    fn rand() -> Self {
        let val: [u64; 4] = [
            rand::random(),
            rand::random(),
            rand::random(),
            rand::random(),
        ];
        let mut ret = Self::default();
        unsafe {
            blst_fr_from_uint64(&mut ret.0, val.as_ptr());
        }

        ret
    }

    fn from_u64_arr(u: &[u64; 4]) -> Self {
        let mut ret = Self::default();
        unsafe {
            blst_fr_from_uint64(&mut ret.0, u.as_ptr());
        }

        ret
    }

    fn from_u64(val: u64) -> Self {
        Self::from_u64_arr(&[val, 0, 0, 0])
    }

    fn to_u64_arr(&self) -> [u64; 4] {
        let mut val: [u64; 4] = [0; 4];
        unsafe {
            blst_uint64_from_fr(val.as_mut_ptr(), &self.0);
        }

        val
    }

    fn is_one(&self) -> bool {
        let mut val: [u64; 4] = [0; 4];
        unsafe {
            blst_uint64_from_fr(val.as_mut_ptr(), &self.0);
        }

        val[0] == 1 && val[1] == 0 && val[2] == 0 && val[3] == 0
    }

    fn is_zero(&self) -> bool {
        let mut val: [u64; 4] = [0; 4];
        unsafe {
            blst_uint64_from_fr(val.as_mut_ptr(), &self.0);
        }

        val[0] == 0 && val[1] == 0 && val[2] == 0 && val[3] == 0
    }

    fn is_null(&self) -> bool {
        self.equals(&Self::null())
    }

    fn sqr(&self) -> Self {
        let mut ret = Self::default();
        unsafe {
            blst_fr_sqr(&mut ret.0, &self.0);
        }

        ret
    }

    fn mul(&self, b: &Self) -> Self {
        let mut ret = Self::default();
        unsafe {
            blst_fr_mul(&mut ret.0, &self.0, &b.0);
        }

        ret
    }

    fn add(&self, b: &Self) -> Self {
        let mut ret = Self::default();
        unsafe {
            blst_fr_add(&mut ret.0, &self.0, &b.0);
        }

        ret
    }

    fn sub(&self, b: &Self) -> Self {
        let mut ret = Self::default();
        unsafe {
            blst_fr_sub(&mut ret.0, &self.0, &b.0);
        }

        ret
    }

    fn eucl_inverse(&self) -> Self {
        let mut ret = Self::default();
        unsafe {
            blst_fr_eucl_inverse(&mut ret.0, &self.0);
        }

        ret
    }

    fn negate(&self) -> Self {
        let mut ret = Self::default();
        unsafe {
            blst_fr_cneg(&mut ret.0, &self.0, true);
        }

        ret
    }

    fn inverse(&self) -> Self {
        let mut ret = Self::default();
        unsafe {
            blst_fr_inverse(&mut ret.0, &self.0);
        }

        ret
    }

    fn pow(&self, n: usize) -> Self {
        let mut out = Self::one();

        let mut temp = *self;
        let mut n = n;
        loop {
            if (n & 1) == 1 {
                out = out.mul(&temp);
            }
            n >>= 1;
            if n == 0 {
                break;
            }

            temp = temp.sqr();
        }

        out
    }

    fn div(&self, b: &Self) -> Result<Self, String> {
        let tmp = b.eucl_inverse();
        let out = self.mul(&tmp);

        Ok(out)
    }

    fn equals(&self, b: &Self) -> bool {
        let mut val_a: [u64; 4] = [0; 4];
        let mut val_b: [u64; 4] = [0; 4];

        unsafe {
            blst_uint64_from_fr(val_a.as_mut_ptr(), &self.0);
            blst_uint64_from_fr(val_b.as_mut_ptr(), &b.0);
        }

        val_a[0] == val_b[0] && val_a[1] == val_b[1] && val_a[2] == val_b[2] && val_a[3] == val_b[3]
    }
}

impl FsFr {
    pub fn to_scalar(&self) -> [u8; 32usize] {
        let mut scalar = blst_scalar::default();
        unsafe {
            blst_scalar_from_fr(&mut scalar, &self.0);
        }

        scalar.b
    }

    pub fn from_scalar(scalar: [u8; 32usize]) -> Result<Self, u8> {
        let mut bls_scalar = blst_scalar::default();
        
        let mut fr = blst_fr::default();
        unsafe {
            blst_scalar_from_lendian(& mut bls_scalar, scalar.as_ptr());
            if !blst_scalar_fr_check(&bls_scalar)
            {
                return Err(1);
            }
            blst_fr_from_scalar(&mut fr, &bls_scalar);
        }
        let mut ret = Self::default();
        ret.0 = fr;
        Ok(ret)
    }

    pub fn hash_to_bls_field(scalar: [u8; 32usize]) -> Self {
        let bls_scalar = blst_scalar { b: scalar };
        let mut fr = blst_fr::default();
        unsafe {
            blst_fr_from_scalar(&mut fr, &bls_scalar);
        }
        let mut ret = Self::default();
        ret.0 = fr;
        ret
    }
}

impl Clone for FsFr {
    fn clone(&self) -> Self {
        FsFr(self.0)
    }
}

impl Copy for FsFr {}
