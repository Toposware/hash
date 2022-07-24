// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use cheetah::Fp;

#[allow(unused)]
/// Exponent of the Anemoi S-Box
pub(crate) const ALPHA: u32 = 7;

#[allow(unused)]
/// Inverse exponent
pub(crate) const INV_ALPHA: u64 = 10540996611094048183;

/// Multiplier of the Anemoi S-Box
pub(crate) const BETA: u32 = 7;

/// First added constant of the Anemoi S-Box
pub(crate) const DELTA: Fp = Fp::new(2635249152773512046);

#[allow(unused)]
/// Second added constant of the Anemoi S-Box
pub(crate) const QUAD: u32 = 2;
