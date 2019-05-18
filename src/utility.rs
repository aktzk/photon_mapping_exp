use crate::ray::Ray;
use rand::rngs::ThreadRng;
use rand::Rng;
use std::f32::consts::PI;
pub fn random_vec_sphere(rand: &mut rand::rngs::ThreadRng) -> glm::Vec3 {
    let r1 = 2.0 * PI * rand.gen::<f32>();
    let r2 = 1.0 - 2.0 * rand.gen::<f32>();
    let sqrt = (1.0 - r2 * r2).sqrt();
    glm::Vec3::new(sqrt * r1.cos(), sqrt * r1.sin(), r2).normalize()
}

pub fn random_vec_hemisphere(normal: &glm::Vec3, rand: &mut rand::rngs::ThreadRng) -> glm::Vec3 {
    let tangent = if normal.x.abs() > 0.1 {
        glm::Vec3::new(0.0, 1.0, 0.0).cross(&normal).normalize()
    } else {
        glm::Vec3::new(1.0, 0.0, 0.0).cross(&normal).normalize()
    };
    let bitangent = normal.cross(&tangent).normalize();

    let r1 = 2.0 * PI * rand.gen::<f32>();
    let r2 = rand.gen::<f32>();
    let r2s = r2.sqrt();

    (tangent * r1.cos() * r2s + bitangent * r1.sin() * r2s + normal * (1.0 - r2).sqrt()).normalize()
}

pub fn primary_ray(
    camera_pos: &glm::Vec3,
    camera_dir: &glm::Vec3,
    (width, height): (usize, usize),
    (x, y): (usize, usize),
    (sx, sy): (usize, usize),
    (cx, cy): (&glm::Vec3, &glm::Vec3),
    rng: &mut ThreadRng,
) -> Ray {
    use rand::Rng;
    let r1 = 2.0 * rng.gen::<f32>();
    let dx = if r1 < 1.0 {
        r1.sqrt() - 1.0
    } else {
        1.0 - (2.0 - r1).sqrt()
    };
    let r2 = 2.0 * rng.gen::<f32>();
    let dy = if r2 < 1.0 {
        r2.sqrt() - 1.0
    } else {
        1.0 - (2.0 - r2).sqrt()
    };
    let dir = cx * (((sx as f32 + 0.5 + dx) / 2.0 + x as f32) / width as f32 - 0.5)
        + cy * (((sy as f32 + 0.5 + dy) / 2.0 + y as f32) / height as f32 - 0.5)
        + camera_dir;
    Ray::new(camera_pos + dir * 130.0, dir.normalize())
}
