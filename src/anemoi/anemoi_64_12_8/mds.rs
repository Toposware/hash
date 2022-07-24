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
    Fp::new(7),
    Fp::new(7683234439643377518),
    Fp::new(4969464316397294165),
    Fp::new(1486652383518389446),
    Fp::new(10374160003859427328),
    Fp::new(2080639838189532388),
    Fp::new(49),
    Fp::new(16889152938674473991),
    Fp::new(5575996515595268031),
    Fp::new(15376031001026020287),
    Fp::new(318796132876043458),
    Fp::new(6491894801771569723),
    Fp::new(2401),
    Fp::new(15911754940807515141),
    Fp::new(5868215007390030733),
    Fp::new(17435674188209071342),
    Fp::new(13491306927794819838),
    Fp::new(13001985994131448344),
    Fp::new(5764801),
    Fp::new(916645121239609502),
    Fp::new(7270777467373168997),
    Fp::new(17628207846428173141),
    Fp::new(11665471366669951938),
    Fp::new(10255705198543237153),
    Fp::new(33232930569601),
    Fp::new(3958698637022038721),
    Fp::new(17681628264411007595),
    Fp::new(9238126393220571964),
    Fp::new(982080019696549145),
    Fp::new(2867565295223319292),
    Fp::new(3732854072722565977),
    Fp::new(11250920084865090207),
    Fp::new(5482877217478222584),
    Fp::new(1482022857694203904),
    Fp::new(16108729179239576902),
    Fp::new(17911738215939738154),
];
