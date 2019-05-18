use glm::Vec3;
pub struct Ray {
    pub pos: Vec3,
    pub dir: Vec3,
}

impl Ray {
    pub fn new(pos: Vec3, dir: Vec3) -> Self {
        Self { pos, dir }
    }
}
