// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use super::traits::AnemoiHasher;
use crate::f64_utils::apply_rescue_inv_sbox;
use cheetah::fp_arith_utils::reduce_u96;
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

/// Function state is set to 12 field elements or 96 bytes.
/// 4 element of the state is reserved for capacity.
pub const STATE_WIDTH: usize = 12;
/// 8 elements of the state are reserved for rate.
pub const RATE_WIDTH: usize = 8;

/// The state is divided into two even-length rows.
pub const NUM_COLUMNS: usize = 6;

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
    let x: [u128; NUM_COLUMNS] = [
        state[0].output_unreduced_internal() as u128,
        state[1].output_unreduced_internal() as u128,
        state[2].output_unreduced_internal() as u128,
        state[3].output_unreduced_internal() as u128,
        state[4].output_unreduced_internal() as u128,
        state[5].output_unreduced_internal() as u128,
    ];
    let y: [u128; NUM_COLUMNS] = [
        state[6].output_unreduced_internal() as u128,
        state[7].output_unreduced_internal() as u128,
        state[8].output_unreduced_internal() as u128,
        state[9].output_unreduced_internal() as u128,
        state[10].output_unreduced_internal() as u128,
        state[6].output_unreduced_internal() as u128,
    ];

    // Fully unroll the matrix-vector product
    state[0] = Fp::from_raw_unchecked(reduce_u96(
        x[0] + x[1] + 3 * x[2] + 4 * x[3] + 5 * x[4] + 6 * x[5],
    ));
    state[1] = Fp::from_raw_unchecked(reduce_u96(
        x[1] + x[2] + 3 * x[3] + 4 * x[4] + 5 * x[5] + 6 * x[0],
    ));
    state[2] = Fp::from_raw_unchecked(reduce_u96(
        x[2] + x[3] + 3 * x[4] + 4 * x[5] + 5 * x[0] + 6 * x[1],
    ));
    state[3] = Fp::from_raw_unchecked(reduce_u96(
        x[3] + x[4] + 3 * x[5] + 4 * x[0] + 5 * x[1] + 6 * x[2],
    ));
    state[4] = Fp::from_raw_unchecked(reduce_u96(
        x[4] + x[5] + 3 * x[0] + 4 * x[1] + 5 * x[2] + 6 * x[3],
    ));
    state[5] = Fp::from_raw_unchecked(reduce_u96(
        x[5] + x[0] + 3 * x[1] + 4 * x[2] + 5 * x[3] + 6 * x[4],
    ));

    state[6] = Fp::from_raw_unchecked(reduce_u96(
        y[0] + y[1] + 3 * y[2] + 4 * y[3] + 5 * y[4] + 6 * y[5],
    ));
    state[7] = Fp::from_raw_unchecked(reduce_u96(
        y[1] + y[2] + 3 * y[3] + 4 * y[4] + 5 * y[5] + 6 * y[0],
    ));
    state[8] = Fp::from_raw_unchecked(reduce_u96(
        y[2] + y[3] + 3 * y[4] + 4 * y[5] + 5 * y[0] + 6 * y[1],
    ));
    state[9] = Fp::from_raw_unchecked(reduce_u96(
        y[3] + y[4] + 3 * y[5] + 4 * y[0] + 5 * y[1] + 6 * y[2],
    ));
    state[10] = Fp::from_raw_unchecked(reduce_u96(
        y[4] + y[5] + 3 * y[0] + 4 * y[1] + 5 * y[2] + 6 * y[3],
    ));
    state[11] = Fp::from_raw_unchecked(reduce_u96(
        y[5] + y[0] + 3 * y[1] + 4 * y[2] + 5 * y[3] + 6 * y[4],
    ));
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

    fn apply_naive_mds(state: &mut [Fp; STATE_WIDTH]) {
        let x: [Fp; NUM_COLUMNS] = state[..NUM_COLUMNS].try_into().unwrap();
        let mut y: [Fp; NUM_COLUMNS] = state[NUM_COLUMNS..].try_into().unwrap();
        y[NUM_COLUMNS - 1] = y[0];

        let mut result = [Fp::zero(); STATE_WIDTH];
        for (i, r) in result.iter_mut().enumerate().take(NUM_COLUMNS) {
            for (j, s) in x.into_iter().enumerate().take(NUM_COLUMNS) {
                *r += s.mul_by_u32(mds::MDS[i * NUM_COLUMNS + j]);
            }
        }
        for (i, r) in result.iter_mut().enumerate().skip(NUM_COLUMNS) {
            for (j, s) in y.into_iter().enumerate() {
                *r += s.mul_by_u32(mds::MDS[(i - NUM_COLUMNS) * NUM_COLUMNS + j]);
            }
        }

        state.copy_from_slice(&result);
    }

    #[test]
    fn test_sbox() {
        let mut input = [
            [Fp::zero(); 12],
            [Fp::one(); 12],
            [
                Fp::new(2024109077186517620),
                Fp::new(12588825917665470320),
                Fp::new(7760860720677708682),
                Fp::new(16177882880946617764),
                Fp::new(18191798079926199855),
                Fp::new(10100602807731852232),
                Fp::new(11835745083059610080),
                Fp::new(892561951383864769),
                Fp::new(13988037933637127990),
                Fp::new(14324215206530819782),
                Fp::new(13639172193883149028),
                Fp::new(1920157329771195333),
            ],
            [
                Fp::new(1928413705048742896),
                Fp::new(4522418955528112507),
                Fp::new(4604024865047881149),
                Fp::new(16185663939685474472),
                Fp::new(14408780645238915608),
                Fp::new(17033445442757683162),
                Fp::new(502097313377030940),
                Fp::new(6339402882007789789),
                Fp::new(12804686218596198210),
                Fp::new(12281970597440949873),
                Fp::new(12477861761310076503),
                Fp::new(3173909208507164602),
            ],
            [
                Fp::new(11913114315485008487),
                Fp::new(18141626625953537205),
                Fp::new(5466594395386693597),
                Fp::new(6469679528976785669),
                Fp::new(11345274226090193764),
                Fp::new(3081164339061692103),
                Fp::new(16079373693054707259),
                Fp::new(14965997961000837161),
                Fp::new(5860717931539236521),
                Fp::new(16577609164319677147),
                Fp::new(12174141410160517753),
                Fp::new(8824030838750282846),
            ],
            [
                Fp::new(2264436143936496420),
                Fp::new(11584707206727149111),
                Fp::new(10824636595394241357),
                Fp::new(15094969360501619238),
                Fp::new(2776107005139276639),
                Fp::new(10973335413158021045),
                Fp::new(13717364785956962097),
                Fp::new(14063318386091763240),
                Fp::new(17749083672895835715),
                Fp::new(6196407931777178051),
                Fp::new(14164789136584134239),
                Fp::new(9065243134059269077),
            ],
            [
                Fp::new(8814214913357488237),
                Fp::new(15724100487956173268),
                Fp::new(7142235699456086627),
                Fp::new(5068653739418772304),
                Fp::new(698706186255199095),
                Fp::new(17840849325155126713),
                Fp::new(7076372013939180650),
                Fp::new(17338632044032178558),
                Fp::new(15709981066789899118),
                Fp::new(15588725907967777404),
                Fp::new(6497338426358533130),
                Fp::new(15158059627781191733),
            ],
            [
                Fp::new(8408319644013278863),
                Fp::new(8810238424382262638),
                Fp::new(17851815567127402984),
                Fp::new(9171415913091609951),
                Fp::new(1771509058229150979),
                Fp::new(4046639637765528370),
                Fp::new(14693926144730479214),
                Fp::new(2203940597110882227),
                Fp::new(12105410459676786746),
                Fp::new(3801101816898933090),
                Fp::new(7208834587551219240),
                Fp::new(15716706316021131551),
            ],
            [
                Fp::new(4998506365808379243),
                Fp::new(15208719414762900578),
                Fp::new(13371019653976871462),
                Fp::new(5726894230476081554),
                Fp::new(9708279985692718518),
                Fp::new(3750738894470073493),
                Fp::new(10003716305067087534),
                Fp::new(971889734796253583),
                Fp::new(17465456643262127616),
                Fp::new(13906978308237797693),
                Fp::new(14210529306864364232),
                Fp::new(11017535886186269578),
            ],
            [
                Fp::new(410290834919184177),
                Fp::new(7161248779663427531),
                Fp::new(4277646107852266495),
                Fp::new(9438475434566402675),
                Fp::new(17559498806804689829),
                Fp::new(3224533887260173732),
                Fp::new(16423108785506381247),
                Fp::new(2532045072365757750),
                Fp::new(2649484840291090519),
                Fp::new(15987649591986014205),
                Fp::new(5182698842403106977),
                Fp::new(2544904436168575168),
            ],
        ];

        // Generated from https://github.com/vesselinux/anemoi-hash/
        let output = [
            [
                Fp::new(2635249152773512046),
                Fp::new(2635249152773512046),
                Fp::new(2635249152773512046),
                Fp::new(2635249152773512046),
                Fp::new(2635249152773512046),
                Fp::new(2635249152773512046),
                Fp::new(0),
                Fp::new(0),
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
                Fp::new(17136669903572771321),
                Fp::new(17136669903572771321),
                Fp::new(9739452640566982996),
                Fp::new(9739452640566982996),
                Fp::new(9739452640566982996),
                Fp::new(9739452640566982996),
                Fp::new(9739452640566982996),
                Fp::new(9739452640566982996),
            ],
            [
                Fp::new(6132711489636036126),
                Fp::new(5918822188063673680),
                Fp::new(10446678021584155829),
                Fp::new(16964551615429726657),
                Fp::new(8496222303710632002),
                Fp::new(3477327031934251233),
                Fp::new(2202638584542667414),
                Fp::new(2364868573641652395),
                Fp::new(7563762104617216708),
                Fp::new(5131023974248620867),
                Fp::new(15577193992643265678),
                Fp::new(2002744822286548109),
            ],
            [
                Fp::new(8516356956060030828),
                Fp::new(17151260213857643203),
                Fp::new(16728524655330841118),
                Fp::new(16092206439755848833),
                Fp::new(17342345795614464189),
                Fp::new(6081008446748763062),
                Fp::new(6656207662753047206),
                Fp::new(16713933356042988271),
                Fp::new(7610483331539526435),
                Fp::new(1211264143768457202),
                Fp::new(13505035650457764834),
                Fp::new(3605241985596149370),
            ],
            [
                Fp::new(16899467629092983873),
                Fp::new(18132463529798524294),
                Fp::new(14539701963499931800),
                Fp::new(4508255734005078794),
                Fp::new(90397383874200573),
                Fp::new(18271567876927277141),
                Fp::new(848877993221171272),
                Fp::new(16277623761987857321),
                Fp::new(2619198227170866453),
                Fp::new(16206721901239863900),
                Fp::new(5173115140969115829),
                Fp::new(13821866989200063063),
            ],
            [
                Fp::new(4114877935645344089),
                Fp::new(1298473822225565460),
                Fp::new(139956930119091677),
                Fp::new(967955194690387694),
                Fp::new(14435407564648906903),
                Fp::new(870477618814990282),
                Fp::new(12267049565242416325),
                Fp::new(2264702414618527143),
                Fp::new(12912204536829586772),
                Fp::new(7674965166183680652),
                Fp::new(17925123739136397906),
                Fp::new(12549281313697703205),
            ],
            [
                Fp::new(3671976149642471943),
                Fp::new(11290391244065142847),
                Fp::new(4356060842796664305),
                Fp::new(16727875796200994647),
                Fp::new(6433951687777342085),
                Fp::new(8862904692471486319),
                Fp::new(8731898714809703440),
                Fp::new(16309570361572560422),
                Fp::new(11580350213066709247),
                Fp::new(9973117683696248088),
                Fp::new(17982611496779959930),
                Fp::new(17468157642886583319),
            ],
            [
                Fp::new(12887623233962299490),
                Fp::new(790477856098858301),
                Fp::new(12101399607856536077),
                Fp::new(5775740306705668201),
                Fp::new(18357832718944242404),
                Fp::new(6419466442382007069),
                Fp::new(5388801016629289639),
                Fp::new(16533913276800006179),
                Fp::new(13539039011294569797),
                Fp::new(9318087012520816981),
                Fp::new(3365216551120229265),
                Fp::new(4414259699169438764),
            ],
            [
                Fp::new(10998501662383288656),
                Fp::new(17065004905754582509),
                Fp::new(14761807523428159219),
                Fp::new(829246465715379829),
                Fp::new(5711697479944667205),
                Fp::new(9351454135009701661),
                Fp::new(4847699396439423268),
                Fp::new(6170721655567503245),
                Fp::new(11851929158458500207),
                Fp::new(7653497002256175893),
                Fp::new(8215553716999357425),
                Fp::new(7748022938793759377),
            ],
            [
                Fp::new(7061252014131344379),
                Fp::new(9569647199443577686),
                Fp::new(10701392684239731145),
                Fp::new(18417440558244288162),
                Fp::new(10214933021117073600),
                Fp::new(14907240941524323696),
                Fp::new(9765561107268105567),
                Fp::new(12042677452152309940),
                Fp::new(3592335798835096123),
                Fp::new(11354882668610971946),
                Fp::new(1220690774487256542),
                Fp::new(12009901643060524681),
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
            [Fp::zero(); 12],
            [Fp::one(); 12],
            [
                Fp::new(1637432090852250421),
                Fp::new(872212631474738575),
                Fp::new(143805946100521899),
                Fp::new(9044503320473486741),
                Fp::new(12032547000700062896),
                Fp::new(8007565852702113475),
                Fp::new(18192157418611634124),
                Fp::new(17337861707317333629),
                Fp::new(13148590453303840855),
                Fp::new(8286719720344966712),
                Fp::new(17303444288291025940),
                Fp::new(7324933286623983425),
            ],
            [
                Fp::new(4891357452182915814),
                Fp::new(2373477572306854047),
                Fp::new(267626666322887225),
                Fp::new(12127854719478737375),
                Fp::new(14573677191043602772),
                Fp::new(248122872627292144),
                Fp::new(8621825208670832062),
                Fp::new(3690305879057571458),
                Fp::new(2808832248456775773),
                Fp::new(9827036983186682347),
                Fp::new(3887636109777739025),
                Fp::new(12285693129667191665),
            ],
            [
                Fp::new(14097654108905235627),
                Fp::new(2240981132181423892),
                Fp::new(8498850823804890413),
                Fp::new(10973818453710135465),
                Fp::new(13398096907044156547),
                Fp::new(89196488989576870),
                Fp::new(5869304038230422818),
                Fp::new(11611115030688719634),
                Fp::new(15493299010329949294),
                Fp::new(2278580744843544238),
                Fp::new(14629944224573630807),
                Fp::new(17070230625408128505),
            ],
            [
                Fp::new(16785894986009148114),
                Fp::new(14872740749783351432),
                Fp::new(14511060816078700503),
                Fp::new(14315351033678091155),
                Fp::new(2943131107662096154),
                Fp::new(16263639011987418134),
                Fp::new(8368092107586295303),
                Fp::new(9294973412841951895),
                Fp::new(14757757745016701792),
                Fp::new(8798405716138394751),
                Fp::new(1239659786614780649),
                Fp::new(5783714904247965800),
            ],
            [
                Fp::new(17407483403065413018),
                Fp::new(1104412903491921002),
                Fp::new(9744666190165019169),
                Fp::new(9655730978269576173),
                Fp::new(3929911743243983987),
                Fp::new(9444961235307621710),
                Fp::new(8414893042018185610),
                Fp::new(741481346187608237),
                Fp::new(1387260118324200819),
                Fp::new(17422534320210012517),
                Fp::new(8583254756347046661),
                Fp::new(17128914923947700271),
            ],
            [
                Fp::new(4434815463382610094),
                Fp::new(7861326706600265258),
                Fp::new(16099206092059814865),
                Fp::new(6555577537247021458),
                Fp::new(17563815378530411139),
                Fp::new(13523156535424665383),
                Fp::new(13731032034104152468),
                Fp::new(330264921768612172),
                Fp::new(8153592736206639878),
                Fp::new(18307741545190040212),
                Fp::new(2362429009179471570),
                Fp::new(3071438693405684843),
            ],
            [
                Fp::new(1671639519009262502),
                Fp::new(5195361531688627050),
                Fp::new(6099073245301310998),
                Fp::new(3358632638730564693),
                Fp::new(6684886421851447096),
                Fp::new(5319837422278702388),
                Fp::new(744417216002137167),
                Fp::new(12410274699695691896),
                Fp::new(10966915251029523160),
                Fp::new(17877256688199808108),
                Fp::new(10618892414989196791),
                Fp::new(12023665021544044844),
            ],
            [
                Fp::new(5375943932951176091),
                Fp::new(16603680963474285123),
                Fp::new(18356478783738838494),
                Fp::new(10266223074231148311),
                Fp::new(15625369611447264579),
                Fp::new(7734048730075843963),
                Fp::new(10368197512847870954),
                Fp::new(3053565235529016530),
                Fp::new(17031229526319444423),
                Fp::new(16030882407381955756),
                Fp::new(3568162697506100987),
                Fp::new(8126446901302821969),
            ],
        ];

        let mut input2 = input;

        // Generated from https://github.com/vesselinux/anemoi-hash/
        let output = [
            [Fp::zero(); 12],
            [Fp::new(20); 12],
            [
                Fp::new(18199997476333406740),
                Fp::new(15461673933932536256),
                Fp::new(16949673644933997499),
                Fp::new(1533142868421963052),
                Fp::new(9640333027543074591),
                Fp::new(1127449743046388051),
                Fp::new(8644930209711942474),
                Fp::new(11079470899500220287),
                Fp::new(9060836800468076499),
                Fp::new(4921468015279211338),
                Fp::new(16632374408791721005),
                Fp::new(9545365062642856940),
            ],
            [
                Fp::new(1809048606453057513),
                Fp::new(17227671820986814163),
                Fp::new(3572936797124383239),
                Fp::new(5144245988363720457),
                Fp::new(2414579941446412995),
                Fp::new(13837803895360947938),
                Fp::new(2088699080857057477),
                Fp::new(17243662325662412807),
                Fp::new(13356575296872086575),
                Fp::new(17138252092615836129),
                Fp::new(15462079185607820126),
                Fp::new(1330434206453821549),
            ],
            [
                Fp::new(5682172441183442005),
                Fp::new(16264885418941386111),
                Fp::new(14830694806318434239),
                Fp::new(14100923381178336821),
                Fp::new(7060648468029141832),
                Fp::new(5695430306328419222),
                Fp::new(15419487806802599340),
                Fp::new(9448325039061801300),
                Fp::new(18131442332369250230),
                Fp::new(6094837688893001174),
                Fp::new(9668903396229953624),
                Fp::new(4803533517945880215),
            ],
            [
                Fp::new(4943039026586359076),
                Fp::new(7881527101732508876),
                Fp::new(17621863677634421324),
                Fp::new(17921886725509229530),
                Fp::new(10801667318053648988),
                Fp::new(3586612493049443403),
                Fp::new(5962860643306932073),
                Fp::new(18328392150280118622),
                Fp::new(10783859711011133019),
                Fp::new(721939026606785217),
                Fp::new(14003704196095803697),
                Fp::new(7508170178829236540),
            ],
            [
                Fp::new(15114192362879671846),
                Fp::new(4291440832811763765),
                Fp::new(15059919220358950677),
                Fp::new(9519825073416538497),
                Fp::new(10651994997426822753),
                Fp::new(11875264974418568950),
                Fp::new(10393227433131534448),
                Fp::new(15272490287998816072),
                Fp::new(14062019746377167305),
                Fp::new(4707287464494682942),
                Fp::new(9106306266734626678),
                Fp::new(15641518179711667190),
            ],
            [
                Fp::new(15966413797960857599),
                Fp::new(5193017618381349001),
                Fp::new(14313452768965530080),
                Fp::new(15415809697381657071),
                Fp::new(11198780120139327699),
                Fp::new(4291883541875308662),
                Fp::new(3037193832214690249),
                Fp::new(2537221872433924262),
                Fp::new(11535546591167465476),
                Fp::new(1339579194450220297),
                Fp::new(6307845021283543261),
                Fp::new(759814934636472555),
            ],
            [
                Fp::new(11708487637380609521),
                Fp::new(10951926328378170179),
                Fp::new(16535105346220764057),
                Fp::new(3027116186633248594),
                Fp::new(14661274309144051054),
                Fp::new(11642615895247674468),
                Fp::new(9104733167813179043),
                Fp::new(16992654634817429455),
                Fp::new(12735043840488217273),
                Fp::new(13985982064542888940),
                Fp::new(3974847629151965639),
                Fp::new(14326028633168652096),
            ],
            [
                Fp::new(2837421079868360440),
                Fp::new(14718773969636060017),
                Fp::new(11575882225241735735),
                Fp::new(5500375155710925228),
                Fp::new(1027289386550679149),
                Fp::new(4729122687025793613),
                Fp::new(5774794757920346905),
                Fp::new(12032824721199469453),
                Fp::new(7827816499707966338),
                Fp::new(6717446266824475078),
                Fp::new(17235726933927595722),
                Fp::new(15117956365361110682),
            ],
        ];

        for i in input.iter_mut() {
            apply_mds(i);
        }
        for i in input2.iter_mut() {
            apply_naive_mds(i);
        }

        assert_eq!(input, output);
        assert_eq!(input2, output);
    }
}
