use crate::ray::Ray;
use glm::Vec3;
#[derive(Eq, PartialEq)]
pub enum Reflection {
    Diffuse,
    Refraction,
}

pub struct Sphere {
    pub radius: f32,
    pub pos: Vec3,
    pub color: Vec3,
    pub emission: Vec3,
    pub reflection: Reflection,
}

impl Sphere {
    pub fn new(
        radius: f32,
        pos: Vec3,
        emission: Vec3,
        color: Vec3,
        reflection: Reflection,
    ) -> Self {
        Self {
            radius,
            pos,
            color,
            emission,
            reflection,
        }
    }

    pub fn intersect(&self, ray: &Ray) -> Option<f32> {
        let ray_to_pos = self.pos - ray.pos;
        let b = ray_to_pos.dot(&ray.dir);
        let det = b * b - ray_to_pos.dot(&ray_to_pos) + self.radius * self.radius;
        if det < 0.0 {
            None
        } else {
            let det_sqrt = det.sqrt();
            let t1 = b - det_sqrt;
            let t2 = b + det_sqrt;
            //自己衝突防ぐ
            if t1 > 0.01 {
                Some(t1)
            } else if t2 > 0.01 {
                Some(t2)
            } else {
                None
            }
        }
    }
}
