// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use super::traits::AnemoiHasher;
use crate::f64_utils::apply_rescue_inv_sbox;
use cheetah::Fp;

/// Digest for Anemoi
mod digest;
/// Hasher for Anemoi
mod hasher;
/// MDS matrix for Anemoi
mod mds;
/// Round constants for Anemoi
mod round_constants;
/// S-Box for Anemoi
mod sbox;

pub use digest::AnemoiDigest;
pub use hasher::AnemoiHash;

// ANEMOI CONSTANTS
// ================================================================================================

/// Function state is set to 8 field elements or 64 bytes.
/// 4 element of the state is reserved for capacity.
pub const STATE_WIDTH: usize = 8;
/// 4 elements of the state are reserved for rate.
pub const RATE_WIDTH: usize = 4;

/// The state is divided into two even-length rows.
pub const NUM_COLUMNS: usize = 4;

/// Four elements (32-bytes) are returned as digest.
pub const DIGEST_SIZE: usize = 4;

/// The number of rounds is set to 10 to provide 128-bit security level.
pub const NUM_HASH_ROUNDS: usize = 10;

// HELPER FUNCTIONS
// ================================================================================================

#[inline(always)]
/// Applies exponentiation of the current hash
/// state elements with the Anemoi S-Box.
pub(crate) fn apply_sbox(state: &mut [Fp; STATE_WIDTH]) {
    let mut x: [Fp; NUM_COLUMNS] = state[..NUM_COLUMNS].try_into().unwrap();
    let mut y: [Fp; NUM_COLUMNS] = state[NUM_COLUMNS..].try_into().unwrap();

    x.iter_mut()
        .enumerate()
        .for_each(|(i, t)| *t -= y[i].square().mul_by_u32(sbox::BETA));

    let mut x_alpha_inv = x;
    apply_rescue_inv_sbox(&mut x_alpha_inv);

    y.iter_mut()
        .enumerate()
        .for_each(|(i, t)| *t -= x_alpha_inv[i]);

    x.iter_mut()
        .enumerate()
        .for_each(|(i, t)| *t += y[i].square().mul_by_u32(sbox::BETA) + sbox::DELTA);

    state[..NUM_COLUMNS].copy_from_slice(&x);
    state[NUM_COLUMNS..].copy_from_slice(&y);
}

#[inline(always)]
/// Applies matrix-vector multiplication of the current
/// hash state with the Anemoi MDS matrix.
pub(crate) fn apply_mds(state: &mut [Fp; STATE_WIDTH]) {
    let x: [Fp; NUM_COLUMNS] = state[..NUM_COLUMNS].try_into().unwrap();
    let mut y: [Fp; NUM_COLUMNS] = state[NUM_COLUMNS..].try_into().unwrap();
    y[NUM_COLUMNS - 1] = y[0];

    let mut result = [Fp::zero(); STATE_WIDTH];
    for (i, r) in result.iter_mut().enumerate().take(NUM_COLUMNS) {
        for (j, s) in x.into_iter().enumerate().take(NUM_COLUMNS) {
            *r += mds::MDS[i * NUM_COLUMNS + j] * s;
        }
    }
    for (i, r) in result.iter_mut().enumerate().skip(NUM_COLUMNS) {
        for (j, s) in y.into_iter().enumerate() {
            *r += mds::MDS[(i - NUM_COLUMNS) * NUM_COLUMNS + j] * s;
        }
    }

    state.copy_from_slice(&result);
}

// ANEMOI PERMUTATION
// ================================================================================================

/// Applies Anemoi permutation to the provided state.
pub(crate) fn apply_permutation(state: &mut [Fp; STATE_WIDTH]) {
    for i in 0..NUM_HASH_ROUNDS {
        apply_round(state, i);
    }

    apply_mds(state)
}

/// Anemoi round function;
/// implementation based on algorithm 3 of <https://eprint.iacr.org/2020/1143.pdf>
#[inline(always)]
pub(crate) fn apply_round(state: &mut [Fp; STATE_WIDTH], step: usize) {
    // determine which round constants to use
    let c = &round_constants::C[step % NUM_HASH_ROUNDS];
    let d = &round_constants::D[step % NUM_HASH_ROUNDS];

    for i in 0..NUM_COLUMNS {
        state[i] += c[i];
        state[NUM_COLUMNS + i] += d[i];
    }

    apply_mds(state);
    apply_sbox(state);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sbox() {
        let mut input = [
            [Fp::zero(); 8],
            [Fp::one(); 8],
            [
                Fp::new(1712038687461379489),
                Fp::new(5447692412612352510),
                Fp::new(18204172923278645249),
                Fp::new(10591898010320281136),
                Fp::new(283191858142073309),
                Fp::new(3864336007846896293),
                Fp::new(1522005613891560833),
                Fp::new(5418759553085432934),
            ],
            [
                Fp::new(1541062408453381426),
                Fp::new(3503646097076694849),
                Fp::new(2712578047928052391),
                Fp::new(6290432062219367636),
                Fp::new(14441238844480069048),
                Fp::new(5445570216321326792),
                Fp::new(5572102543327390213),
                Fp::new(3272656753080261139),
            ],
            [
                Fp::new(17750372577009177710),
                Fp::new(18045199254715991266),
                Fp::new(431807966204104042),
                Fp::new(4422564390015731990),
                Fp::new(11440835893337523470),
                Fp::new(3925012699882739076),
                Fp::new(9237338518877523097),
                Fp::new(14847589435329253892),
            ],
            [
                Fp::new(13256793983916649931),
                Fp::new(10010473922917342458),
                Fp::new(7626915519600657249),
                Fp::new(1668659765529187623),
                Fp::new(4177593121955056366),
                Fp::new(1932519349258285252),
                Fp::new(9286001150948808617),
                Fp::new(818375147643082505),
            ],
            [
                Fp::new(18111899899168656778),
                Fp::new(11637224498565880810),
                Fp::new(3926654893371677306),
                Fp::new(1084585967661441637),
                Fp::new(8380925107748770016),
                Fp::new(18393066941408618868),
                Fp::new(15681834392335680385),
                Fp::new(2175405604383493142),
            ],
            [
                Fp::new(5632466407657881660),
                Fp::new(11484723426861672517),
                Fp::new(13407788852044947478),
                Fp::new(10114073292181808127),
                Fp::new(10549748759110055225),
                Fp::new(14369187703589314050),
                Fp::new(15515036773990820875),
                Fp::new(1757341577814839666),
            ],
            [
                Fp::new(8731916342278781126),
                Fp::new(4236018243673470033),
                Fp::new(8418751023433932793),
                Fp::new(13237799449732800883),
                Fp::new(14143145891008054791),
                Fp::new(10651931506109164093),
                Fp::new(2587335361737723073),
                Fp::new(14668526727523831857),
            ],
            [
                Fp::new(11061481425177694301),
                Fp::new(7713462632215422249),
                Fp::new(9518551352150821073),
                Fp::new(15946096186210043374),
                Fp::new(13464420365798483472),
                Fp::new(6419978687819605030),
                Fp::new(9767046550469939299),
                Fp::new(13372081843923258540),
            ],
        ];

        // Generated from https://github.com/vesselinux/anemoi-hash/
        let output = [
            [
                Fp::new(2635249152773512046),
                Fp::new(2635249152773512046),
                Fp::new(2635249152773512046),
                Fp::new(2635249152773512046),
                Fp::new(0),
                Fp::new(0),
                Fp::new(0),
                Fp::new(0),
            ],
            [
                Fp::new(17136669903572771321),
                Fp::new(17136669903572771321),
                Fp::new(17136669903572771321),
                Fp::new(17136669903572771321),
                Fp::new(9739452640566982996),
                Fp::new(9739452640566982996),
                Fp::new(9739452640566982996),
                Fp::new(9739452640566982996),
            ],
            [
                Fp::new(11065593969245908309),
                Fp::new(2286127855641544998),
                Fp::new(15689296389681640249),
                Fp::new(14314911564785099190),
                Fp::new(18041492866718121394),
                Fp::new(12942790176000775814),
                Fp::new(11279554698674879841),
                Fp::new(13271526095498080630),
            ],
            [
                Fp::new(6464878638129929081),
                Fp::new(17524289919992718869),
                Fp::new(3859655523574647086),
                Fp::new(10586293103549178185),
                Fp::new(13448306645114613967),
                Fp::new(11801119000867309936),
                Fp::new(15773183161778793863),
                Fp::new(6267068725400869727),
            ],
            [
                Fp::new(1071212341937840955),
                Fp::new(8972563262231276501),
                Fp::new(16370595962907239102),
                Fp::new(4729734543326922629),
                Fp::new(6672399665957534038),
                Fp::new(921357344921007409),
                Fp::new(11387755116491454845),
                Fp::new(14769019070528731400),
            ],
            [
                Fp::new(7983873117537237800),
                Fp::new(11563440969981879555),
                Fp::new(4440141217050288717),
                Fp::new(9497684814048276266),
                Fp::new(8172488069878816018),
                Fp::new(4096679275156481593),
                Fp::new(382314221310022416),
                Fp::new(17633812099912861859),
            ],
            [
                Fp::new(2141561183840690073),
                Fp::new(11351850202517186388),
                Fp::new(3255304693446410421),
                Fp::new(9808407631757139858),
                Fp::new(11200959537443332717),
                Fp::new(15169191600341221135),
                Fp::new(16033662830699201455),
                Fp::new(2828893045206078095),
            ],
            [
                Fp::new(18160027176437807363),
                Fp::new(8044717837117114504),
                Fp::new(14536520699326170092),
                Fp::new(11723066257545918551),
                Fp::new(5183732065004339153),
                Fp::new(6702839739455794500),
                Fp::new(16646447688720198307),
                Fp::new(476231656731423553),
            ],
            [
                Fp::new(6005143981396154000),
                Fp::new(1446080273200796136),
                Fp::new(3398891077492945395),
                Fp::new(5728872265745593742),
                Fp::new(15548717630340029217),
                Fp::new(1348876669730540859),
                Fp::new(10203403694170648433),
                Fp::new(69805472060322904),
            ],
            [
                Fp::new(5185980573057915351),
                Fp::new(3447315659597812394),
                Fp::new(2574211670649801739),
                Fp::new(2959227954784150207),
                Fp::new(10052038333208711028),
                Fp::new(836346455852744406),
                Fp::new(13359109125587763740),
                Fp::new(9740298391612936083),
            ],
        ];

        for i in input.iter_mut() {
            apply_sbox(i);
        }

        assert_eq!(input, output);
    }

    #[test]
    fn test_mds() {
        let mut input = [
            [Fp::zero(); 8],
            [Fp::one(); 8],
            [
                Fp::new(1712038687461379489),
                Fp::new(5447692412612352510),
                Fp::new(18204172923278645249),
                Fp::new(10591898010320281136),
                Fp::new(283191858142073309),
                Fp::new(3864336007846896293),
                Fp::new(1522005613891560833),
                Fp::new(5418759553085432934),
            ],
            [
                Fp::new(1541062408453381426),
                Fp::new(3503646097076694849),
                Fp::new(2712578047928052391),
                Fp::new(6290432062219367636),
                Fp::new(14441238844480069048),
                Fp::new(5445570216321326792),
                Fp::new(5572102543327390213),
                Fp::new(3272656753080261139),
            ],
            [
                Fp::new(17750372577009177710),
                Fp::new(18045199254715991266),
                Fp::new(431807966204104042),
                Fp::new(4422564390015731990),
                Fp::new(11440835893337523470),
                Fp::new(3925012699882739076),
                Fp::new(9237338518877523097),
                Fp::new(14847589435329253892),
            ],
            [
                Fp::new(13256793983916649931),
                Fp::new(10010473922917342458),
                Fp::new(7626915519600657249),
                Fp::new(1668659765529187623),
                Fp::new(4177593121955056366),
                Fp::new(1932519349258285252),
                Fp::new(9286001150948808617),
                Fp::new(818375147643082505),
            ],
            [
                Fp::new(18111899899168656778),
                Fp::new(11637224498565880810),
                Fp::new(3926654893371677306),
                Fp::new(1084585967661441637),
                Fp::new(8380925107748770016),
                Fp::new(18393066941408618868),
                Fp::new(15681834392335680385),
                Fp::new(2175405604383493142),
            ],
            [
                Fp::new(5632466407657881660),
                Fp::new(11484723426861672517),
                Fp::new(13407788852044947478),
                Fp::new(10114073292181808127),
                Fp::new(10549748759110055225),
                Fp::new(14369187703589314050),
                Fp::new(15515036773990820875),
                Fp::new(1757341577814839666),
            ],
            [
                Fp::new(8731916342278781126),
                Fp::new(4236018243673470033),
                Fp::new(8418751023433932793),
                Fp::new(13237799449732800883),
                Fp::new(14143145891008054791),
                Fp::new(10651931506109164093),
                Fp::new(2587335361737723073),
                Fp::new(14668526727523831857),
            ],
            [
                Fp::new(11061481425177694301),
                Fp::new(7713462632215422249),
                Fp::new(9518551352150821073),
                Fp::new(15946096186210043374),
                Fp::new(13464420365798483472),
                Fp::new(6419978687819605030),
                Fp::new(9767046550469939299),
                Fp::new(13372081843923258540),
            ],
        ];

        // Generated from https://github.com/vesselinux/anemoi-hash/
        let output = [
            [
                Fp::new(0),
                Fp::new(0),
                Fp::new(0),
                Fp::new(0),
                Fp::new(0),
                Fp::new(0),
                Fp::new(0),
                Fp::new(0),
            ],
            [
                Fp::new(107),
                Fp::new(128),
                Fp::new(23),
                Fp::new(38),
                Fp::new(107),
                Fp::new(128),
                Fp::new(23),
                Fp::new(38),
            ],
            [
                Fp::new(9456771577905369261),
                Fp::new(4590966864761370959),
                Fp::new(339316553332415136),
                Fp::new(10473575290158421696),
                Fp::new(8225049216658878461),
                Fp::new(2346599146747010196),
                Fp::new(17954635621241173182),
                Fp::new(601986623308693410),
            ],
            [
                Fp::new(5971364692563835175),
                Fp::new(11530372389151605268),
                Fp::new(11775231841592517333),
                Fp::new(6239000710601362686),
                Fp::new(5754328324182719121),
                Fp::new(6242990758763377164),
                Fp::new(11506335194229375048),
                Fp::new(10837563009177541526),
            ],
            [
                Fp::new(17720293982075903351),
                Fp::new(4099342423628860977),
                Fp::new(4856055662473052934),
                Fp::new(9490462113027361833),
                Fp::new(10052969982707184867),
                Fp::new(13232320694465843769),
                Fp::new(16341701930518921075),
                Fp::new(9238531938690773830),
            ],
            [
                Fp::new(5394439459632089150),
                Fp::new(13925236541443090345),
                Fp::new(7721442454914416581),
                Fp::new(3724619832167511906),
                Fp::new(15456776455890992924),
                Fp::new(13683739330411441528),
                Fp::new(9445483374727542473),
                Fp::new(18365743720473918638),
            ],
            [
                Fp::new(14655426929987996076),
                Fp::new(9130470654784838799),
                Fp::new(10038923116734400839),
                Fp::new(9390200411471243362),
                Fp::new(11104548798238667367),
                Fp::new(17381295994192183494),
                Fp::new(3458160390868646738),
                Fp::new(10555721932437915352),
            ],
            [
                Fp::new(15006069552515823961),
                Fp::new(2710726675864932558),
                Fp::new(12597169402255728592),
                Fp::new(12491387952638033474),
                Fp::new(9748236402592025386),
                Fp::new(18221360052868400713),
                Fp::new(1463859989027341610),
                Fp::new(18288995327637913779),
            ],
            [
                Fp::new(15294987321235138145),
                Fp::new(14267684876884659229),
                Fp::new(11627466822756923910),
                Fp::new(2768186725754038649),
                Fp::new(1249138628157422771),
                Fp::new(15562284673230158002),
                Fp::new(9105668843505038076),
                Fp::new(5243218392953797103),
            ],
            [
                Fp::new(5321500469830995058),
                Fp::new(8561109272663518353),
                Fp::new(2024878717784182296),
                Fp::new(9467864308655838677),
                Fp::new(10453560567059391640),
                Fp::new(4087623600585799285),
                Fp::new(9821088271815951974),
                Fp::new(7467476819739155502),
            ],
        ];

        for i in input.iter_mut() {
            apply_mds(i);
        }

        assert_eq!(input, output);
    }
}
