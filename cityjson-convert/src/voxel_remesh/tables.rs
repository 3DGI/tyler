//! Triangulation rules ported from meshoptimizer v1.3 (MIT; see LICENSE.md).
#![forbid(unsafe_code)]

use std::sync::LazyLock;

pub(super) struct Tables {
    pub triangles: [[[u16; 5]; 2]; 256],
    pub counts: [u8; 256],
    pub modes: [u8; 256],
}

pub(super) static TABLES: LazyLock<Tables> = LazyLock::new(Tables::build);

// Cube corners, rather than edge intersections, are encoded in three nibbles.
const BASE: &[(usize, u8, [u16; 5], [u16; 5])] = &[
    (0x00, 0, [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]),
    (0x01, 0, [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]),
    (0x03, 0, [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]),
    (0x06, 0, [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]),
    (0x07, 0, [0x120, 0, 0, 0, 0], [0, 0, 0, 0, 0]),
    (0x0f, 2, [0x130, 0x032, 0, 0, 0], [0x321, 0x120, 0, 0, 0]),
    (0x16, 1, [0x124, 0x142, 0, 0, 0], [0, 0, 0, 0, 0]),
    (0x17, 0, [0x124, 0, 0, 0, 0], [0, 0, 0, 0, 0]),
    (0x18, 0, [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]),
    (0x19, 1, [0x043, 0x403, 0, 0, 0], [0, 0, 0, 0, 0]),
    (0x1b, 1, [0x134, 0x430, 0, 0, 0], [0x410, 0x013, 0, 0, 0]),
    (0x1d, 1, [0x243, 0x340, 0, 0, 0], [0x320, 0x024, 0, 0, 0]),
    (0x1e, 1, [0x142, 0x432, 0x341, 0, 0], [0x132, 0, 0, 0, 0]),
    (0x1f, 2, [0x412, 0x213, 0, 0, 0], [0x134, 0x432, 0, 0, 0]),
    (0x3c, 1, [0x345, 0x324, 0x543, 0x423, 0], [0, 0, 0, 0, 0]),
    (
        0x3d,
        1,
        [0x253, 0x524, 0x350, 0, 0],
        [0x240, 0x504, 0x320, 0, 0],
    ),
    (0x3f, 2, [0x325, 0x524, 0, 0, 0], [0x243, 0x345, 0, 0, 0]),
    (0x69, 1, [0x035, 0x560, 0x306, 0x536, 0], [0, 0, 0, 0, 0]),
    (0x6b, 1, [0x063, 0x605, 0x365, 0, 0], [0x305, 0, 0, 0, 0]),
    (0x6f, 1, [0x365, 0x560, 0, 0, 0], [0x530, 0x036, 0, 0, 0]),
    (0x7e, 0, [0x142, 0x536, 0, 0, 0], [0, 0, 0, 0, 0]),
    (0x7f, 0, [0x365, 0, 0, 0, 0], [0, 0, 0, 0, 0]),
    (0xf0, 2, [0x547, 0x746, 0, 0, 0], [0x465, 0x567, 0, 0, 0]),
    (0xff, 0, [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]),
];

impl Tables {
    fn build() -> Self {
        let mut table = Self {
            triangles: [[[0; 5]; 2]; 256],
            counts: [0; 256],
            modes: [0; 256],
        };
        let mut filled = [0; 256];
        for &(code, mode, primary, alternate) in BASE {
            table.triangles[code] = [primary, if mode == 0 { primary } else { alternate }];
            table.modes[code] = mode;
            filled[code] = 1;
        }
        let rotations = [[2, 3, 6, 7, 0, 1, 4, 5], [4, 0, 6, 2, 5, 1, 7, 3]];
        for _ in 0..3 {
            for code in 0..256 {
                if filled[code] != 1 {
                    continue;
                }
                for rotation in rotations {
                    let mut rotated = 0;
                    for (corner, destination) in rotation.iter().enumerate() {
                        rotated |= ((code >> corner) & 1) << destination;
                    }
                    if filled[rotated] != 0 {
                        continue;
                    }
                    for variant in 0..2 {
                        for index in 0..5 {
                            let triangle = table.triangles[code][variant][index];
                            if triangle == 0 {
                                break;
                            }
                            table.triangles[rotated][variant][index] =
                                (rotation[usize::from((triangle >> 8) & 15)] << 8)
                                    | (rotation[usize::from((triangle >> 4) & 15)] << 4)
                                    | rotation[usize::from(triangle & 15)];
                        }
                    }
                    table.modes[rotated] = table.modes[code];
                    filled[rotated] = 1;
                }
                filled[code] = 2;
            }
        }
        for (code, state) in filled.iter().enumerate() {
            assert_ne!(*state, 0);
            let count = table.triangles[code][0]
                .iter()
                .take_while(|&&triangle| triangle != 0)
                .count();
            table.counts[code] = u8::try_from(count).expect("at most four triangles per cube");
            assert_eq!(table.triangles[code][1][count], 0);
        }
        table
    }
}
