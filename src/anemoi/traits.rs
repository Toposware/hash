// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::traits::Hasher;
use group::ff::Field;

pub trait AnemoiHasher<F: Field>: Hasher<F> {
    /// Initializes a new instance of the permutation.
    fn new() -> Self;

    /// Absorbs a sequence of bytes.
    fn absorb(&mut self, input: &[u8]);

    /// Absorbs a sequence of field elements.
    fn absorb_field(&mut self, input: &[F]);

    /// Returns hash of the data absorbed into the hasher.
    fn finalize(&mut self) -> Self::Digest;
}
