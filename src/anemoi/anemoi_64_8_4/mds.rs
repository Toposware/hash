// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use super::NUM_COLUMNS;
use cheetah::Fp;

/// Maximum Diffusion Layer matrix for Anemoi,
pub(crate) const MDS: [Fp; NUM_COLUMNS * NUM_COLUMNS] = [
    Fp::one(),
    Fp::new(49),
    Fp::new(49),
    Fp::new(8),
    Fp::new(8),
    Fp::new(56),
    Fp::new(49),
    Fp::new(15),
    Fp::new(7),
    Fp::new(8),
    Fp::one(),
    Fp::new(7),
    Fp::new(7),
    Fp::new(15),
    Fp::new(8),
    Fp::new(8),
];
