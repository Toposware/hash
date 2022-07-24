// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! Hasher trait implementation for Anemoi

use core::convert::TryInto;

use super::digest::AnemoiDigest;
use super::AnemoiHasher;
use super::{apply_permutation, DIGEST_SIZE, RATE_WIDTH, STATE_WIDTH};
use crate::error::SerializationError;
use crate::traits::Hasher;

use cheetah::Fp;

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
/// A Anemoi Hash over Fp
pub struct AnemoiHash {
    state: [Fp; STATE_WIDTH],
    idx: usize,
}

impl Default for AnemoiHash {
    fn default() -> Self {
        Self {
            state: [Fp::zero(); STATE_WIDTH],
            idx: 0,
        }
    }
}

impl AnemoiHash {
    /// Serializes the current state to an array of bytes
    pub fn to_bytes(&self) -> [u8; 104] {
        let mut res = [0u8; 104];
        assert_eq!(res.len(), STATE_WIDTH * 8 + 8);

        for (index, elem) in self.state.iter().enumerate() {
            res[index * 8..index * 8 + 8].copy_from_slice(&elem.to_bytes());
        }
        res[96..104].copy_from_slice(&(self.idx as u64).to_le_bytes());

        res
    }

    /// Returns a AnemoiHash from an array of bytes
    pub fn from_bytes(bytes: &[u8; 104]) -> Result<Self, SerializationError> {
        let mut state = [Fp::zero(); STATE_WIDTH];
        let mut array = [0u8; 8];
        for index in 0..STATE_WIDTH {
            array.copy_from_slice(&bytes[index * 8..index * 8 + 8]);
            let value = Fp::from_bytes(&array);
            state[index] = match value.is_some().into() {
                true => value.unwrap(),
                false => return Err(SerializationError::InvalidFieldElement),
            };
        }

        array.copy_from_slice(&bytes[96..104]);
        let idx = u64::from_le_bytes(array) as usize;

        Ok(Self { state, idx })
    }
}

impl Hasher<Fp> for AnemoiHash {
    type Digest = AnemoiDigest;

    fn hash(bytes: &[u8]) -> Self::Digest {
        // compute the number of elements required to represent the string; we will be processing
        // the string in 7-byte chunks, thus the number of elements will be equal to the number
        // of such chunks (including a potential partial chunk at the end).
        let num_elements = if bytes.len() % 7 == 0 {
            bytes.len() / 7
        } else {
            bytes.len() / 7 + 1
        };

        // initialize state to all zeros, except for the last element of the capacity part, which
        // is set to the number of elements to be hashed. this is done so that adding zero elements
        // at the end of the list always results in a different hash.
        let mut state = [Fp::zero(); STATE_WIDTH];
        state[STATE_WIDTH - 1] = Fp::new(num_elements as u64);

        // break the string into 7-byte chunks, convert each chunk into a field element, and
        // absorb the element into the rate portion of the state. we use 7-byte chunks because
        // every 7-byte chunk is guaranteed to map to some field element.
        let mut i = 0;
        let mut num_hashed = 0;
        let mut buf = [0u8; 8];
        for chunk in bytes.chunks(7) {
            if num_hashed + i < num_elements - 1 {
                buf[..7].copy_from_slice(chunk);
            } else {
                // if we are dealing with the last chunk, it may be smaller than 7 bytes long, so
                // we need to handle it slightly differently. we also append a byte with value 1
                // to the end of the string; this pads the string in such a way that adding
                // trailing zeros results in different hash
                let chunk_len = chunk.len();
                buf = [0u8; 8];
                buf[..chunk_len].copy_from_slice(chunk);
                buf[chunk_len] = 1;
            }

            // convert the bytes into a field element and absorb it into the rate portion of the
            // state; if the rate is filled up, apply the Anemoi permutation and start absorbing
            // again from zero index.
            state[i] += Fp::new(u64::from_le_bytes(buf));
            i += 1;
            if i % RATE_WIDTH == 0 {
                apply_permutation(&mut state);
                i = 0;
                num_hashed += RATE_WIDTH;
            }
        }

        // if we absorbed some elements but didn't apply a permutation to them (would happen when
        // the number of elements is not a multiple of RATE_WIDTH), apply the Anemoi permutation.
        // we don't need to apply any extra padding because we injected total number of elements
        // in the input list into the capacity portion of the state during initialization.
        if i > 0 {
            apply_permutation(&mut state);
        }

        // return the first DIGEST_SIZE elements of the state as hash result
        AnemoiDigest::new(state[..DIGEST_SIZE].try_into().unwrap())
    }

    fn hash_field(bytes: &[Fp]) -> Self::Digest {
        // initialize state to all zeros
        let mut state = [Fp::zero(); STATE_WIDTH];

        let mut i = 0;
        for &element in bytes.iter() {
            state[i] += element;
            i += 1;
            if i % RATE_WIDTH == 0 {
                apply_permutation(&mut state);
                i = 0;
            }
        }

        // Apply padding specification
        if i > 0 {
            state[i] += Fp::one();
            i += 1;

            while i % RATE_WIDTH != 0 {
                state[i] = Fp::zero();
                i += 1;
            }

            apply_permutation(&mut state);
        }

        AnemoiDigest::new(state[..DIGEST_SIZE].try_into().unwrap())
    }

    fn merge(values: &[Self::Digest; 2]) -> Self::Digest {
        let mut state = [Fp::zero(); STATE_WIDTH];
        state[..DIGEST_SIZE].copy_from_slice(values[0].as_elements());
        state[DIGEST_SIZE..RATE_WIDTH].copy_from_slice(values[1].as_elements());
        apply_permutation(&mut state);

        AnemoiDigest::new(state[..DIGEST_SIZE].try_into().unwrap())
    }
}

impl AnemoiHasher<Fp> for AnemoiHash {
    /// Initializes a new instance of the permutation.
    fn new() -> Self {
        Self::default()
    }

    /// Absorbs a sequence of bytes.
    fn absorb(&mut self, input: &[u8]) {
        // compute the number of elements required to represent the string; we will be processing
        // the string in 7-byte chunks, thus the number of elements will be equal to the number
        // of such chunks (including a potential partial chunk at the end).
        let num_elements = if input.len() % 7 == 0 {
            input.len() / 7
        } else {
            input.len() / 7 + 1
        };

        // break the string into 7-byte chunks, convert each chunk into a field element, and
        // absorb the element into the rate portion of the state. we use 7-byte chunks because
        // every 7-byte chunk is guaranteed to map to some field element.
        let mut num_hashed = 0;
        let mut buf = [0u8; 8];
        for chunk in input.chunks(7) {
            if num_hashed + self.idx < num_elements - 1 {
                buf[..7].copy_from_slice(chunk);
            } else {
                // if we are dealing with the last chunk, it may be smaller than 7 bytes long, so
                // we need to handle it slightly differently. we also append a byte with value 1
                // to the end of the string; this pads the string in such a way that adding
                // trailing zeros results in different hash

                // Compatibility with the binary hash() is not possible because this would require
                // knowing the total input sequence length at initialization, to write in the capacity
                // registers. Hence, we prevent length-extension attacks on every absorbed chunk
                let chunk_len = chunk.len();
                buf = [0u8; 8];
                buf[..chunk_len].copy_from_slice(chunk);
                buf[chunk_len] = 1;
            }

            // convert the bytes into a field element and absorb it into the rate portion of the
            // state; if the rate is filled up, apply the Anemoi permutation and start absorbing
            // again from zero index.
            self.state[self.idx] += Fp::new(u64::from_le_bytes(buf));
            self.idx += 1;
            if self.idx % RATE_WIDTH == 0 {
                apply_permutation(&mut self.state);
                self.idx = 0;
                num_hashed += RATE_WIDTH;
            }
        }
    }

    /// Absorbs a sequence of field elements.
    fn absorb_field(&mut self, input: &[Fp]) {
        for &element in input {
            self.state[self.idx] += element;
            self.idx += 1;
            if self.idx % RATE_WIDTH == 0 {
                apply_permutation(&mut self.state);
                self.idx = 0;
            }
        }
    }

    /// Returns hash of the data absorbed into the hasher.
    fn finalize(&mut self) -> Self::Digest {
        // Apply padding specification
        if self.idx > 0 {
            self.state[self.idx] += Fp::one();
            self.idx += 1;

            while self.idx % RATE_WIDTH != 0 {
                self.state[self.idx] += Fp::zero();
                self.idx += 1;
            }

            apply_permutation(&mut self.state);
            self.idx = 0;
        }

        AnemoiDigest::new(self.state[..DIGEST_SIZE].try_into().unwrap())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::OsRng;

    #[test]
    fn test_anemoi_hash() {
        // Hardcoded input / output list generated from the
        // Sagemath code at https://github.com/vesselinux/anemoi-hash/

        let input_data = [
            [Fp::zero(); 8],
            [Fp::one(); 8],
            [
                Fp::new(17984842363570004217),
                Fp::new(10852236884266818836),
                Fp::new(10676388419824119985),
                Fp::new(11301224415378875499),
                Fp::new(1358243899414205366),
                Fp::new(4292865320243257527),
                Fp::new(18270623391838632347),
                Fp::new(8444434354056942962),
            ],
            [
                Fp::new(8712806641306796433),
                Fp::new(10268031888921461830),
                Fp::new(10559025976246513094),
                Fp::new(14162087600471753434),
                Fp::new(7676579746637058731),
                Fp::new(17334364963435480710),
                Fp::new(1630109610567533489),
                Fp::new(4166170351851415947),
            ],
            [
                Fp::new(4794268187954129850),
                Fp::new(12059731652215576642),
                Fp::new(8647552991167658020),
                Fp::new(2790028884491538111),
                Fp::new(1570356508244096438),
                Fp::new(11826893810771925487),
                Fp::new(2120176184902868026),
                Fp::new(15690904052937280985),
            ],
            [
                Fp::new(7467591281694523315),
                Fp::new(14835812323019073600),
                Fp::new(2529421407595381399),
                Fp::new(18005513861558773540),
                Fp::new(9148631752947046111),
                Fp::new(15994963396633558467),
                Fp::new(1249975758071850329),
                Fp::new(11768757815743855626),
            ],
            [
                Fp::new(184118663655219488),
                Fp::new(5175577586329124244),
                Fp::new(10092167467118086843),
                Fp::new(14507513766515168398),
                Fp::new(10202650830708651764),
                Fp::new(5068557337482218631),
                Fp::new(6951691851729231799),
                Fp::new(14876954462140975976),
            ],
            [
                Fp::new(17045445415261233849),
                Fp::new(635417663701319433),
                Fp::new(9574544570370214838),
                Fp::new(14261736828690442921),
                Fp::new(10807201351559864903),
                Fp::new(18269877510678660252),
                Fp::new(2480683274878739233),
                Fp::new(8758175479803205565),
            ],
            [
                Fp::new(8886495645822610101),
                Fp::new(1591626565372620151),
                Fp::new(7763431313759214136),
                Fp::new(670517855511071710),
                Fp::new(9791495431883229480),
                Fp::new(8993521367454606087),
                Fp::new(17468155451152355939),
                Fp::new(941657751076229314),
            ],
            [
                Fp::new(13607910476244079099),
                Fp::new(10101750930321417662),
                Fp::new(15483623777947699986),
                Fp::new(17160114433809822048),
                Fp::new(12307130960029289407),
                Fp::new(16923943649005112423),
                Fp::new(773584190033609472),
                Fp::new(14973933375664266091),
            ],
        ];

        let output_data = [
            [
                Fp::new(6155559623514467584),
                Fp::new(11046471219574225279),
                Fp::new(6806677178655161577),
                Fp::new(2233980693767959800),
            ],
            [
                Fp::new(5586115931723519951),
                Fp::new(17178352092800175766),
                Fp::new(11367055701620428725),
                Fp::new(11725629169683932889),
            ],
            [
                Fp::new(13463880730571272115),
                Fp::new(9400118443279441563),
                Fp::new(14470818421305374514),
                Fp::new(17812478870892306776),
            ],
            [
                Fp::new(2028467862697714118),
                Fp::new(7605561620699906335),
                Fp::new(7395304823865194828),
                Fp::new(15683379984392346780),
            ],
            [
                Fp::new(137778207588373523),
                Fp::new(834554531589735485),
                Fp::new(4823896296740388981),
                Fp::new(12840137232119835898),
            ],
            [
                Fp::new(12151096573622576980),
                Fp::new(254356417596022964),
                Fp::new(12431210437747696798),
                Fp::new(5740682134354625308),
            ],
            [
                Fp::new(6662632532307331099),
                Fp::new(13584961009394474970),
                Fp::new(641591515456263539),
                Fp::new(17062130891955436756),
            ],
            [
                Fp::new(15038388590439401742),
                Fp::new(7426490138705799526),
                Fp::new(2804710318407548938),
                Fp::new(6277819945416926803),
            ],
            [
                Fp::new(13216766884435952254),
                Fp::new(17300400485696324883),
                Fp::new(1577389083909716875),
                Fp::new(11225934879501371052),
            ],
            [
                Fp::new(7521806294696492405),
                Fp::new(14120968937053013606),
                Fp::new(6207947916212339848),
                Fp::new(18399401376790620484),
            ],
        ];

        for (input, expected) in input_data.iter().zip(output_data) {
            let mut hasher = AnemoiHash::new();
            hasher.absorb_field(input);

            assert_eq!(expected, hasher.finalize().to_elements());
            assert_eq!(expected, AnemoiHash::hash_field(input).to_elements());
        }
    }

    #[test]
    fn test_sequential_hashing() {
        let mut rng = OsRng;

        for _ in 0..100 {
            let mut data = [Fp::zero(); 160];
            for e in data.iter_mut() {
                *e = Fp::random(&mut rng);
            }

            let mut hasher = AnemoiHash::new();
            for chunk in data.chunks(8) {
                hasher.absorb_field(chunk);
            }

            assert_eq!(hasher.finalize(), AnemoiHash::hash_field(&data));
        }
    }

    #[test]
    fn test_serialization() {
        let mut rng = OsRng;

        for _ in 0..100 {
            let mut data = [Fp::zero(); DIGEST_SIZE];
            for e in data.iter_mut() {
                *e = Fp::random(&mut rng);
            }

            let mut hasher = AnemoiHash::new();
            hasher.absorb_field(&data);

            let bytes = hasher.to_bytes();

            assert_eq!(hasher, AnemoiHash::from_bytes(&bytes).unwrap());
        }

        // Test invalid encoding
        let mut data = [Fp::zero(); DIGEST_SIZE];
        for e in data.iter_mut() {
            *e = Fp::random(&mut rng);
        }

        let mut hasher = AnemoiHash::new();
        hasher.absorb_field(&data);

        let bytes = [255u8; 104];

        assert!(AnemoiHash::from_bytes(&bytes).is_err());
    }
}
