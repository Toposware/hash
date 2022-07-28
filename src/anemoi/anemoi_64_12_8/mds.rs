// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use super::NUM_COLUMNS;

/// Maximum Diffusion Layer matrix for Anemoi,
/// This matrix is the first MDS matrix being outputted by this exhaustive
/// SageMath search among all circulant matrices:
///
/// sage: from itertools import combinations_with_replacement
/// sage: for l in combinations_with_replacement(range(0,7), 6):
///           if is_mds(matrix.circulant(list(l))):
///               print(i)
///               break
// Because of the small sizes of the coefficients of the matrix,
// those are stored as u32 to enjoy faster multiplication with
// `Fp` elements.
#[allow(unused)]
pub(crate) const MDS: [u32; NUM_COLUMNS * NUM_COLUMNS] = [
    1, 1, 3, 4, 5, 6, 6, 1, 1, 3, 4, 5, 5, 6, 1, 1, 3, 4, 4, 5, 6, 1, 1, 3, 3, 4, 5, 6, 1, 1, 1, 3,
    4, 5, 6, 1,
];
