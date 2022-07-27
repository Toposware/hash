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
                Fp::new(1954858876411629165),
                Fp::new(2818494584435130042),
                Fp::new(6783907228711662769),
                Fp::new(6646437732877752266),
                Fp::new(665371416343971011),
                Fp::new(17922757374239957108),
                Fp::new(11623914693490540367),
                Fp::new(10401638767690912515),
            ],
            [
                Fp::new(6150065054911601402),
                Fp::new(5449191287923068830),
                Fp::new(5734734849509679934),
                Fp::new(16510419923240473245),
                Fp::new(11493189615583258437),
                Fp::new(9735729831484839646),
                Fp::new(5094076912951934488),
                Fp::new(7489135876442975633),
            ],
            [
                Fp::new(18107789852855413844),
                Fp::new(12506344626766034983),
                Fp::new(167374782600459863),
                Fp::new(5574074256491935097),
                Fp::new(3537762868466561649),
                Fp::new(5485970734139262805),
                Fp::new(4898425914480226622),
                Fp::new(1393194757633584103),
            ],
            [
                Fp::new(10707575619529978689),
                Fp::new(15887976604584119102),
                Fp::new(8057066235251795766),
                Fp::new(17502385575944172508),
                Fp::new(6207801137788422023),
                Fp::new(10584312182628920065),
                Fp::new(15487777107392015324),
                Fp::new(5382489007368416990),
            ],
            [
                Fp::new(9050759089401277385),
                Fp::new(11412116684310639371),
                Fp::new(8841625585153623580),
                Fp::new(13822245615727198631),
                Fp::new(2513317181586056033),
                Fp::new(12996821361402367226),
                Fp::new(11405878517965701115),
                Fp::new(7983362093716068058),
            ],
        ];

        let output_data = [
            [
                Fp::new(7978594262959769940),
                Fp::new(4891354612598643906),
                Fp::new(1424109007685561175),
                Fp::new(11040610241142283265),
            ],
            [
                Fp::new(10808514591873376767),
                Fp::new(4593493611920905725),
                Fp::new(12049638854682712424),
                Fp::new(17832871376549447047),
            ],
            [
                Fp::new(15786155975842809504),
                Fp::new(10708298196492519241),
                Fp::new(3814101100908864102),
                Fp::new(8153390099050784510),
            ],
            [
                Fp::new(13660325862315895822),
                Fp::new(1558771665817158253),
                Fp::new(6504833200094474751),
                Fp::new(17103652055620766798),
            ],
            [
                Fp::new(12743823521645607064),
                Fp::new(8580579411597161414),
                Fp::new(2346072150915238250),
                Fp::new(17723342489889136508),
            ],
            [
                Fp::new(12436689980941808290),
                Fp::new(589333582780954485),
                Fp::new(16476806977505234761),
                Fp::new(15516724128038112601),
            ],
            [
                Fp::new(7751128593403983592),
                Fp::new(15238706365551341165),
                Fp::new(1566422974340761854),
                Fp::new(6002247307586655584),
            ],
            [
                Fp::new(11475867083295323660),
                Fp::new(15546443306660494809),
                Fp::new(13515839470720320970),
                Fp::new(12587452541853591902),
            ],
            [
                Fp::new(17677008490425618360),
                Fp::new(17147166035375234080),
                Fp::new(414188505166238040),
                Fp::new(9985702984736094118),
            ],
            [
                Fp::new(13197551130977259157),
                Fp::new(9856180230431691217),
                Fp::new(14350247721079505569),
                Fp::new(5144649016578241967),
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
