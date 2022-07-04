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
    pub fn to_bytes(&self) -> [u8; 72] {
        let mut res = [0u8; 72];
        assert_eq!(res.len(), STATE_WIDTH * 8 + 8);

        for (index, elem) in self.state.iter().enumerate() {
            res[index * 8..index * 8 + 8].copy_from_slice(&elem.to_bytes());
        }
        res[64..72].copy_from_slice(&(self.idx as u64).to_le_bytes());

        res
    }

    /// Returns a AnemoiHash from an array of bytes
    pub fn from_bytes(bytes: &[u8; 72]) -> Result<Self, SerializationError> {
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

        array.copy_from_slice(&bytes[64..72]);
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
        state[..RATE_WIDTH].copy_from_slice(values[0].as_elements());
        state[RATE_WIDTH..STATE_WIDTH].copy_from_slice(values[1].as_elements());
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
                Fp::new(11537543289854239743),
                Fp::new(15211617747006487828),
                Fp::new(2664939301251259440),
                Fp::new(7626873324951209813),
                Fp::new(9941180307347331261),
                Fp::new(2840344071115909949),
                Fp::new(5823229232739669681),
                Fp::new(628981740748381054),
            ],
            [
                Fp::new(6503150440414379279),
                Fp::new(3167450922377693239),
                Fp::new(2453832157088976847),
                Fp::new(2571822692977591383),
                Fp::new(17074233473382357693),
                Fp::new(812474399550104630),
                Fp::new(6247720078278335874),
                Fp::new(930644360578429399),
            ],
            [
                Fp::new(9705739363138877133),
                Fp::new(15653212624154705200),
                Fp::new(8669609798514298023),
                Fp::new(14648065771611451658),
                Fp::new(12326404382605042768),
                Fp::new(16982707784029275247),
                Fp::new(5152778279860118363),
                Fp::new(3589464785313198640),
            ],
            [
                Fp::new(533098780191878770),
                Fp::new(10577067137382791622),
                Fp::new(4882085344241734715),
                Fp::new(6008379056067445454),
                Fp::new(13972064846016069751),
                Fp::new(12192562426017919124),
                Fp::new(3761341220749356588),
                Fp::new(8248728561637501882),
            ],
            [
                Fp::new(7157563675420492164),
                Fp::new(11306681171822683565),
                Fp::new(1840812812002368146),
                Fp::new(14033180946347444391),
                Fp::new(7161981864287297146),
                Fp::new(2053301154803094047),
                Fp::new(12161358187387510936),
                Fp::new(13428022460351489085),
            ],
            [
                Fp::new(9052137509714315417),
                Fp::new(11014936732502188986),
                Fp::new(10979612237679423201),
                Fp::new(17137565686958396927),
                Fp::new(3410562732024104069),
                Fp::new(10439202231038923597),
                Fp::new(2336259394386946893),
                Fp::new(18144339585556906347),
            ],
            [
                Fp::new(10415015670289733822),
                Fp::new(16697825496378845836),
                Fp::new(17931480217273612690),
                Fp::new(15548492952630527142),
                Fp::new(10088247497895431502),
                Fp::new(18383763776271801686),
                Fp::new(3807632813275665888),
                Fp::new(17325442000004973583),
            ],
            [
                Fp::new(11440126530297973171),
                Fp::new(14845681488757718696),
                Fp::new(2471489167613578153),
                Fp::new(9871436222121912093),
                Fp::new(15475160933618006632),
                Fp::new(7667919044216829439),
                Fp::new(3550950114022886363),
                Fp::new(12331697617309768066),
            ],
        ];

        let output_data = [
            [
                Fp::new(11474871846796477803),
                Fp::new(4144704552913832387),
                Fp::new(10542094298089029483),
                Fp::new(9110399001298944547),
            ],
            [
                Fp::new(2845656333646540438),
                Fp::new(17047782626247437315),
                Fp::new(124493436826911145),
                Fp::new(16015697058842868295),
            ],
            [
                Fp::new(6852127272490302260),
                Fp::new(10660478178540059490),
                Fp::new(17320165645914424207),
                Fp::new(13496544386923020985),
            ],
            [
                Fp::new(2047498148869288349),
                Fp::new(1535703037067685060),
                Fp::new(2071132230268933947),
                Fp::new(14407727394151338471),
            ],
            [
                Fp::new(7790084708662227226),
                Fp::new(16693722354583284304),
                Fp::new(2393265998621550888),
                Fp::new(5400577202186225718),
            ],
            [
                Fp::new(5560062451032549491),
                Fp::new(9399054561030152082),
                Fp::new(4903195221769924979),
                Fp::new(6023972987201240771),
            ],
            [
                Fp::new(16075988500875386140),
                Fp::new(6373161351136927623),
                Fp::new(3015320823670087660),
                Fp::new(9401741029174706432),
            ],
            [
                Fp::new(8821262909878733903),
                Fp::new(13312822718767200372),
                Fp::new(14357917015068908588),
                Fp::new(3241003008577219579),
            ],
            [
                Fp::new(1628238088292292851),
                Fp::new(1018867324057487080),
                Fp::new(11864223016012536912),
                Fp::new(5890826092818440359),
            ],
            [
                Fp::new(16527427468407157944),
                Fp::new(495606115022430282),
                Fp::new(11945455405446995073),
                Fp::new(10716176196518003164),
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
            let mut data = [Fp::zero(); 120];
            for e in data.iter_mut() {
                *e = Fp::random(&mut rng);
            }

            let mut hasher = AnemoiHash::new();
            for chunk in data.chunks(4) {
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

        let bytes = [255u8; 72];

        assert!(AnemoiHash::from_bytes(&bytes).is_err());
    }
}
