use exhact::matrix::{RingSpec, MajorDimension, RingMetadata};
use exhact::csm::{CSM, transpose};
use num::rational::Ratio;

fn main() {
///////////////////////////////////////////////////////////////////////////////////////////////////
/// Example of tranpose()
    let ringmetadata = RingMetadata{
        ringspec: RingSpec::Modulus(2),
        identity_additive: 0,
        identity_multiplicative: 1,
    };

    // 1 2 3 4 0
    // 0 0 0 0 0
    // 0 5 0 0 0

    let matrix: CSM<usize, usize> = CSM {
        ringmetadata,
        nummaj: 3,
        majdim: MajorDimension::Row,
        majptr: vec![0, 4, 4, 5],
        minind: vec![0, 1, 2, 3, 1],
        snzval: vec![1, 2, 3, 4, 5]
    };

    let trans: CSM<usize, usize> = transpose(4, &matrix);
    trans.print();
///////////////////////////////////////////////////////////////////////////////////////////////////





}
