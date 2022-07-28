// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use super::NUM_COLUMNS;

/// Maximum Diffusion Layer matrix for Anemoi,
// Because of the small sizes of the coefficients of the matrix,
// those are stored as u32 to enjoy faster multiplication with
// `Fp` elements.
#[allow(unused)]
pub(crate) const MDS: [u32; NUM_COLUMNS * NUM_COLUMNS] =
    [1, 49, 49, 8, 8, 56, 49, 15, 7, 8, 1, 7, 7, 15, 8, 8];
