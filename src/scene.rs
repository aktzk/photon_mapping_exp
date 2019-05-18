use crate::ray::Ray;
use crate::sphere::{Reflection, Sphere};
use glm::Vec3;

// 光源である球のインデックス
pub const LIGHT_ID: usize = 0;

pub struct Scene {
    pub spheres: Vec<Sphere>,
}

impl Scene {
    // cornell box(smallpt)
    pub fn new() -> Self {
        let mut spheres = Vec::new();
        spheres.push(Sphere::new(
            7.0,
            Vec3::new(60.0, 55.0, 71.6),
            Vec3::new(17.0, 12.0, 4.0),
            Vec3::new(0.78, 0.78, 0.78),
            Reflection::Diffuse,
        ));
        spheres.push(Sphere::new(
            13.0,
            Vec3::new(50.0, 15.0, 81.6),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.99, 0.99, 0.99),
            Reflection::Refraction,
        ));
        spheres.push(Sphere::new(
            1e5,
            Vec3::new(1e5 + 1.0, 40.8, 81.6),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.14, 0.45, 0.091),
            Reflection::Diffuse,
        ));
        spheres.push(Sphere::new(
            1e5,
            Vec3::new(1e5 + 99.0, 40.8, 81.6),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.63, 0.065, 0.05),
            Reflection::Diffuse,
        ));
        spheres.push(Sphere::new(
            1e5,
            Vec3::new(50.0, 40.8, 1e5),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.725, 0.71, 0.68),
            Reflection::Diffuse,
        ));
        spheres.push(Sphere::new(
            1e5,
            Vec3::new(50.0, 1e5, 81.6),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.725, 0.71, 0.68),
            Reflection::Diffuse,
        ));
        spheres.push(Sphere::new(
            1e5,
            Vec3::new(50.0, -1e5 + 81.6, 81.6),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.725, 0.71, 0.68),
            Reflection::Diffuse,
        ));
        Self { spheres }
    }

    pub fn intersect(&self, ray: &Ray) -> Option<(f32, usize)> {
        use std::f32::INFINITY;
        let mut t = INFINITY;
        let mut id = -1;
        for i in 0..self.spheres.len() {
            let p = &self.spheres[i];
            if let Some(d) = p.intersect(ray) {
                if d < t {
                    t = d;
                    id = i as isize;
                }
            }
        }

        if id >= 0 {
            Some((t, id as usize))
        } else {
            None
        }
    }
}
