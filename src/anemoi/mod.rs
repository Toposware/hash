// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

mod traits;

/// The Anemoi hash function over Cheetah's small
/// primefield with state width 8 and rate 4.
#[cfg(feature = "f64")]
pub mod anemoi_64_8_4;
