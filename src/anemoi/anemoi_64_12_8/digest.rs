// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! Digest trait implementation for Anemoi

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

use super::DIGEST_SIZE;
use crate::traits::Digest;

use cheetah::Fp;

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
/// A Anemoi Digest for the Anemoi Hash over Fp
pub struct AnemoiDigest([Fp; DIGEST_SIZE]);

impl AnemoiDigest {
    /// Returns a new Digest from a provided array
    pub fn new(value: [Fp; DIGEST_SIZE]) -> Self {
        Self(value)
    }

    /// Returns a reference to the wrapped digest
    pub fn as_elements(&self) -> &[Fp; DIGEST_SIZE] {
        &self.0
    }

    /// Returns the wrapped digest
    pub fn to_elements(&self) -> [Fp; DIGEST_SIZE] {
        self.0
    }

    /// Returns a `Vec<Fp>` from the provided digest slice
    pub fn digests_to_elements(digests: &[Self]) -> Vec<Fp> {
        let mut res = Vec::with_capacity(digests.len() * DIGEST_SIZE);
        for digest in digests {
            for i in 0..DIGEST_SIZE {
                res.push(digest.0[i]);
            }
        }

        res
    }
}

impl Default for AnemoiDigest {
    fn default() -> Self {
        AnemoiDigest([Fp::default(); DIGEST_SIZE])
    }
}

impl Digest for AnemoiDigest {
    fn to_bytes(&self) -> [u8; 32] {
        let mut digest = [0u8; 32];
        digest[0..8].copy_from_slice(&self.0[0].to_bytes());
        digest[8..16].copy_from_slice(&self.0[1].to_bytes());
        digest[16..24].copy_from_slice(&self.0[2].to_bytes());
        digest[24..32].copy_from_slice(&self.0[3].to_bytes());

        digest
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::OsRng;

    #[test]
    fn digest_elements() {
        let mut rng = OsRng;

        for _ in 0..100 {
            let mut array = [Fp::zero(); DIGEST_SIZE];
            for item in array.iter_mut() {
                *item = Fp::random(&mut rng);
            }

            let digest = AnemoiDigest::new(array);
            assert_eq!(digest.to_elements(), array);
            assert_eq!(&digest.to_elements(), digest.as_elements());
            assert_eq!(
                digest.as_elements(),
                &AnemoiDigest::digests_to_elements(&[digest])[..]
            );
        }

        let digest = AnemoiDigest::default();
        assert_eq!(digest.to_elements(), [Fp::zero(); DIGEST_SIZE]);
        assert_eq!(digest.as_elements(), &vec![Fp::zero(); DIGEST_SIZE][..]);
        assert_eq!(digest.to_bytes(), [0u8; 32]);
    }
}
