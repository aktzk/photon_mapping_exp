use crate::ray::Ray;
use glm::DVec3;
#[derive(Eq, PartialEq)]
pub enum Reflection {
    Diffuse,
    Refraction,
}

pub struct Sphere {
    pub radius: f64,
    pub pos: DVec3,
    pub color: DVec3,
    pub emission: DVec3,
    pub reflection: Reflection,
}

impl Sphere {
    pub fn new(
        radius: f64,
        pos: DVec3,
        emission: DVec3,
        color: DVec3,
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

    pub fn intersect(&self, ray: &Ray) -> Option<f64> {
        let ray_to_pos = self.pos - ray.pos;
        let b = ray_to_pos.dot(&ray.dir);
        let det = b * b - ray_to_pos.dot(&ray_to_pos) + self.radius * self.radius;
        if det < 0.0 {
            None
        } else {
            let det_sqrt = det.sqrt();
            let t1 = b - det_sqrt;
            let t2 = b + det_sqrt;
            if t1 > 0.00001 {
                Some(t1)
            } else if t2 > 0.00001 {
                Some(t2)
            } else {
                None
            }
        }
    }
}
