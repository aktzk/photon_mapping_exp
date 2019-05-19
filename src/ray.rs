use glm::DVec3;
pub struct Ray {
    pub pos: DVec3,
    pub dir: DVec3,
}

impl Ray {
    pub fn new(pos: DVec3, dir: DVec3) -> Self {
        Self { pos, dir }
    }

    // 物体に近いところに原点があるレイ（拡散など）が自己衝突するのを防ぐときに使う
    pub fn new_with_offset(pos: DVec3, dir: DVec3, offset_normal: &DVec3) -> Self {
        Self {
            pos: offset_ray(&pos, offset_normal),
            dir: dir,
        }
    }
}

fn offset_ray(pos: &glm::DVec3, normal: &glm::DVec3) -> glm::DVec3 {
    use std::mem;
    const ORIGIN: f64 = 1.0 / 32.0;
    const FLOAT_SCALE: f64 = 1.0 / 65536.0;
    const INT_SCALE: f64 = 256.0;
    let posx_asi: i64 = unsafe { mem::transmute(pos.x) };
    let posy_asi: i64 = unsafe { mem::transmute(pos.y) };
    let posz_asi: i64 = unsafe { mem::transmute(pos.z) };

    let of_f = INT_SCALE * normal;
    let of_i = (of_f.x as i64, of_f.y as i64, of_f.z as i64);
    let p_f = (
        (posx_asi + if pos.x < 0.0 { -of_i.0 } else { of_i.0 }),
        (posy_asi + if pos.y < 0.0 { -of_i.1 } else { of_i.1 }),
        (posz_asi + if pos.z < 0.0 { -of_i.2 } else { of_i.2 }),
    );
    let pix: f64 = unsafe { mem::transmute(p_f.0) };
    let piy: f64 = unsafe { mem::transmute(p_f.1) };
    let piz: f64 = unsafe { mem::transmute(p_f.2) };

    glm::DVec3::new(
        if pos.x.abs() < ORIGIN {
            pos.x + FLOAT_SCALE * normal.x
        } else {
            pix
        },
        if pos.y.abs() < ORIGIN {
            pos.y + FLOAT_SCALE * normal.y
        } else {
            piy
        },
        if pos.z.abs() < ORIGIN {
            pos.z + FLOAT_SCALE * normal.z
        } else {
            piz
        },
    )
}
