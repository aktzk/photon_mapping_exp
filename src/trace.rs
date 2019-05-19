use crate::photon_map::{Photon, PhotonCollectArg, PhotonMap};
use crate::ray::Ray;
use crate::scene::{Scene, LIGHT_ID};
use crate::sphere::Reflection;
use crate::utility::{random_vec_hemisphere, random_vec_sphere};
use rand::Rng;
use std::f64::consts::PI;
use std::sync::Arc;

const REFRACTION_ATM: f64 = 1.000_293; //1atm0度の空気の絶対屈折率
const REFRACTION_DIA: f64 = 2.419; //ダイヤモンドの絶対屈折率
pub fn create_photon_map_caustics(
    photon_num: usize,
    scene: &Scene,
    rand: &mut rand::rngs::ThreadRng,
) -> PhotonMap {
    println!(
        "Constructing PhotonMap(Caustics) Photon Num: {}",
        photon_num
    );
    let mut photons = Vec::new();
    for _ in 0..photon_num {
        let light = &scene.spheres[LIGHT_ID];
        let random_dir = random_vec_sphere(rand);
        let light_pos = light.pos + random_dir * light.radius;
        let light_dir = random_vec_hemisphere(&random_dir, rand);

        let mut current_ray = Ray::new_with_offset(light_pos, light_dir, &random_dir);
        let mut current_flux =
            light.emission * 4.0 * PI * PI * light.radius * light.radius / (photon_num as f64);

        let mut reflect_count = 0;

        for _ in 0..8 {
            if let Some((dist, idx)) = scene.intersect(&current_ray) {
                let sphere = &scene.spheres[idx];
                let hit_pos = current_ray.pos + current_ray.dir * dist;
                let geometry_normal = (hit_pos - sphere.pos).normalize();
                let surface_normal = if geometry_normal.dot(&current_ray.dir) < 0.0 {
                    geometry_normal
                } else {
                    -geometry_normal
                };

                match sphere.reflection {
                    Reflection::Diffuse => {
                        if reflect_count > 0 {
                            photons.push(Arc::new(Photon::new(
                                hit_pos,
                                current_flux,
                                current_ray.dir,
                            )));
                        }
                        break; // 最初に拡散面に当たる場合photonを格納しない
                    }
                    Reflection::Refraction => {
                        reflect_count += 1;
                        match calc_refraction_behavior(
                            &current_ray.dir,
                            &geometry_normal,
                            &surface_normal,
                        ) {
                            Refraction::TotalReflection(dir) => {
                                current_ray = Ray::new_with_offset(hit_pos, dir, &surface_normal);
                                current_flux = glm::matrix_comp_mult(&current_flux, &sphere.color);
                            }
                            Refraction::ReflectionOrTransmission(
                                (reflect_dir, rp),
                                (refract_dir, tp),
                                probability,
                            ) => {
                                if rand.gen::<f64>() < probability {
                                    current_ray = Ray::new_with_offset(
                                        hit_pos,
                                        reflect_dir,
                                        &geometry_normal,
                                    );
                                    current_flux *= rp;
                                } else {
                                    current_ray = Ray::new_with_offset(
                                        hit_pos,
                                        refract_dir,
                                        &(-geometry_normal),
                                    );
                                    current_flux *= tp;
                                }
                            }
                        }
                        continue;
                    }
                }
            } else {
                break;
            }
        }
    }
    println!("Done. photon_num: {}", photons.len());
    PhotonMap::new(&mut photons)
}

// 一回以上拡散反射したフォトンのみ格納
pub fn create_photon_map_indirect(
    photon_num: usize,
    scene: &Scene,
    rand: &mut rand::rngs::ThreadRng,
) -> PhotonMap {
    println!(
        "Constructing PhotonMap(Indirect) Photon Num: {}",
        photon_num
    );
    let mut photons = Vec::new();
    for _ in 0..photon_num {
        let light = &scene.spheres[LIGHT_ID];
        let random_dir = random_vec_sphere(rand);
        let (light_pos, light_dir) = {
            (
                light.pos + random_dir * light.radius,
                random_vec_hemisphere(&random_dir, rand),
            )
        };

        let mut current_ray = Ray::new_with_offset(light_pos, light_dir, &random_dir);
        let mut current_flux =
            light.emission * 4.0 * PI * PI * light.radius * light.radius / (photon_num as f64);
        let mut scatter_count = 0;
        let max_depth = 8;
        for _ in 0..max_depth {
            if let Some((dist, idx)) = scene.intersect(&current_ray) {
                if idx == LIGHT_ID {
                    break;
                }
                let sphere = &scene.spheres[idx];
                let hit_pos = current_ray.pos + current_ray.dir * dist;
                let geometry_normal = (hit_pos - sphere.pos).normalize();
                let surface_normal = if geometry_normal.dot(&current_ray.dir) < 0.0 {
                    geometry_normal
                } else {
                    -geometry_normal
                };

                match sphere.reflection {
                    Reflection::Diffuse => {
                        let probability = (sphere.color.x + sphere.color.y + sphere.color.z) / 3.0;
                        if probability < rand.gen::<f64>() && scatter_count > 0 {
                            photons.push(Arc::new(Photon::new(
                                hit_pos,
                                current_flux,
                                current_ray.dir,
                            )));
                            break;
                        } else {
                            let dir = random_vec_hemisphere(&surface_normal, rand);
                            current_ray = Ray::new_with_offset(hit_pos, dir, &surface_normal);
                            current_flux =
                                glm::matrix_comp_mult(&current_flux, &sphere.color) / probability;
                            scatter_count += 1;
                            continue;
                        }
                    }
                    Reflection::Refraction => {
                        match calc_refraction_behavior(
                            &current_ray.dir,
                            &geometry_normal,
                            &surface_normal,
                        ) {
                            Refraction::TotalReflection(dir) => {
                                current_ray = Ray::new_with_offset(hit_pos, dir, &surface_normal);
                                current_flux = glm::matrix_comp_mult(&current_flux, &sphere.color);
                            }
                            Refraction::ReflectionOrTransmission(
                                (reflect_dir, rp),
                                (refract_dir, tp),
                                probability,
                            ) => {
                                if rand.gen::<f64>() < probability {
                                    current_ray =
                                        Ray::new_with_offset(hit_pos, reflect_dir, &surface_normal);
                                    current_flux *= rp;
                                } else {
                                    current_ray = Ray::new_with_offset(
                                        hit_pos,
                                        refract_dir,
                                        &(-surface_normal),
                                    );
                                    current_flux *= tp;
                                }
                            }
                        }
                        continue;
                    }
                }
            } else {
                break;
            }
        }
    }

    println!("Done. photon_num: {}", photons.len());
    PhotonMap::new(&mut photons)
}
fn estimate_irradiance(
    photon_map: &PhotonMap,
    pos: &glm::DVec3,
    normal: &glm::DVec3,
    collect_radius: f64,
    max_photon_num: usize,
) -> glm::DVec3 {
    let arg = PhotonCollectArg::new(
        collect_radius * collect_radius,
        max_photon_num,
        *pos,
        *normal,
    );
    let nearest_photons = photon_map.collect_photons(&arg);

    let mut accumulated_flux = glm::zero::<glm::DVec3>();
    let max_distance2 = nearest_photons
        .queue
        .peek()
        .and_then(|a| Some(a.distance2.into_inner()));
    if let Some(max_dist2) = max_distance2 {
        for neighbor in &nearest_photons.queue {
            if neighbor.photon.incoming_dir.dot(normal) < 0.0 {
                let weight = 1.0 - (neighbor.distance2.into_inner() / max_dist2).sqrt();
                accumulated_flux += neighbor.photon.power * weight;
            }
        }
        accumulated_flux / (max_dist2 * PI)
    } else {
        glm::zero()
    }
}

pub fn estimate_radiance_caustics(
    scene: &Scene,
    ray: &Ray,
    photon_map: &PhotonMap,
    collect_radius: f64,
    collect_max_photon_num: usize,
) -> glm::DVec3 {
    if let Some((dist, idx)) = scene.intersect(ray) {
        let sphere = &scene.spheres[idx];
        let hit_pos = ray.pos + ray.dir * dist;
        let geometry_normal = (hit_pos - sphere.pos).normalize();
        let surface_normal = if geometry_normal.dot(&ray.dir) < 0.0 {
            geometry_normal
        } else {
            -geometry_normal
        };

        if sphere.reflection == Reflection::Diffuse {
            if idx == 0 {
                return glm::zero();
            }
            let irradiance = estimate_irradiance(
                photon_map,
                &hit_pos,
                &surface_normal,
                collect_radius,
                collect_max_photon_num,
            );
            return irradiance;
        }
    }
    glm::zero()
}
pub fn estimate_radiance_indirect(
    scene: &Scene,
    ray: &Ray,
    photon_map: &PhotonMap,
    depth: usize,
    col_rad: f64,
    col_max_num: usize,
    rand: &mut rand::rngs::ThreadRng,
) -> glm::DVec3 {
    if let Some((dist, idx)) = scene.intersect(ray) {
        if idx == 0 {
            return glm::zero();
        }
        let sphere = &scene.spheres[idx];
        let hit_pos = ray.pos + ray.dir * dist;
        let geometry_normal = (hit_pos - sphere.pos).normalize();
        let surface_normal = if geometry_normal.dot(&ray.dir) < 0.0 {
            geometry_normal
        } else {
            -geometry_normal
        };

        // ロシアンルーレット
        // ロシアンルーレットを実行しない場合1
        let mut russian_roulette_probability =
            glm::max3_scalar(sphere.color.x, sphere.color.y, sphere.color.z);
        if depth > 8 {
            if rand.gen::<f64>() >= russian_roulette_probability {
                return sphere.emission;
            }
        } else {
            russian_roulette_probability = 1.0;
        }
        match sphere.reflection {
            Reflection::Diffuse => {
                if depth == 0 {
                    let mut accum = glm::zero::<glm::DVec3>();
                    for _ in 0..8 {
                        let scatter_dir = random_vec_hemisphere(&surface_normal, rand);
                        let scatter_ray =
                            Ray::new_with_offset(hit_pos, scatter_dir, &surface_normal);
                        accum += estimate_radiance_indirect(
                            scene,
                            &scatter_ray,
                            photon_map,
                            depth + 1,
                            col_rad,
                            col_max_num,
                            rand,
                        );
                    }
                    accum /= 8.0;
                    return accum;
                } else {
                    let irradiance = estimate_irradiance(
                        photon_map,
                        &hit_pos,
                        &surface_normal,
                        col_rad,
                        col_max_num,
                    );
                    return glm::matrix_comp_mult(&sphere.color, &irradiance)
                        * std::f64::consts::FRAC_1_PI
                        / russian_roulette_probability;
                }
            }
            Reflection::Refraction => {
                if depth == 0 {
                    return glm::zero();
                } else {
                    match calc_refraction_behavior(&ray.dir, &geometry_normal, &surface_normal) {
                        Refraction::TotalReflection(dir) => {
                            let ray = Ray::new_with_offset(hit_pos, dir, &surface_normal);
                            return sphere.emission
                                + glm::matrix_comp_mult(
                                    &sphere.color,
                                    &estimate_radiance_indirect(
                                        scene,
                                        &ray,
                                        photon_map,
                                        depth + 1,
                                        col_rad,
                                        col_max_num,
                                        rand,
                                    ),
                                );
                        }
                        Refraction::ReflectionOrTransmission(
                            (reflect_dir, rp),
                            (refract_dir, tp),
                            probability,
                        ) => {
                            if depth > 2 {
                                if rand.gen::<f64>() < probability {
                                    let ray =
                                        Ray::new_with_offset(hit_pos, reflect_dir, &surface_normal);
                                    let rad = glm::matrix_comp_mult(
                                        &sphere.color,
                                        &(estimate_radiance_indirect(
                                            scene,
                                            &ray,
                                            photon_map,
                                            depth + 1,
                                            col_rad,
                                            col_max_num,
                                            rand,
                                        ) * rp),
                                    ) / probability
                                        / russian_roulette_probability;
                                    return sphere.emission + rad;
                                } else {
                                    let ray = Ray::new_with_offset(
                                        hit_pos,
                                        refract_dir,
                                        &(-surface_normal),
                                    );
                                    let rad = glm::matrix_comp_mult(
                                        &sphere.color,
                                        &(estimate_radiance_indirect(
                                            scene,
                                            &ray,
                                            photon_map,
                                            depth + 1,
                                            col_rad,
                                            col_max_num,
                                            rand,
                                        ) * tp),
                                    ) / (1.0 - probability)
                                        / russian_roulette_probability;
                                    return sphere.emission + rad;
                                }
                            } else {
                                let ray_re =
                                    Ray::new_with_offset(hit_pos, reflect_dir, &surface_normal);
                                let ray_tr =
                                    Ray::new_with_offset(hit_pos, refract_dir, &(-surface_normal));
                                let rad = estimate_radiance_indirect(
                                    scene,
                                    &ray_re,
                                    photon_map,
                                    depth + 1,
                                    col_rad,
                                    col_max_num,
                                    rand,
                                ) * rp
                                    + estimate_radiance_indirect(
                                        scene,
                                        &ray_tr,
                                        photon_map,
                                        depth + 1,
                                        col_rad,
                                        col_max_num,
                                        rand,
                                    ) * tp;
                                return sphere.emission
                                    + glm::matrix_comp_mult(&sphere.color, &rad)
                                        / russian_roulette_probability;
                            }
                        }
                    }
                }
            }
        }
    }

    glm::zero()
}

pub fn estimate_radiance_direct(
    scene: &Scene,
    ray: &Ray,
    depth: usize,
    rand: &mut rand::rngs::ThreadRng,
) -> glm::DVec3 {
    if let Some((t, id)) = scene.intersect(ray) {
        let sphere = &scene.spheres[id];
        if id == LIGHT_ID {
            //光源は無視
            return sphere.emission;
        }
        let hit_pos = ray.pos + t * ray.dir;
        let geometry_normal = (hit_pos - sphere.pos).normalize();
        let surface_normal = if geometry_normal.dot(&ray.dir) < 0.0 {
            geometry_normal
        } else {
            -1.0 * geometry_normal
        };

        match sphere.reflection {
            Reflection::Diffuse => {
                let mut accumulated_flux = glm::zero::<glm::DVec3>();
                for _ in 0..128 {
                    let random_dir = random_vec_sphere(rand);

                    let light_pos =
                        scene.spheres[LIGHT_ID].pos + scene.spheres[LIGHT_ID].radius * random_dir;
                    let light_normal = (light_pos - scene.spheres[0].pos).normalize();
                    // シャドウレイ
                    let ray_dir = (light_pos - hit_pos).normalize();
                    let dist2 = glm::distance2(&light_pos, &hit_pos);
                    let dot0 = surface_normal.dot(&ray_dir);
                    let dot1 = light_normal.dot(&(-1.0 * ray_dir));

                    if dot0 >= 0.0 && dot1 >= 0.0 {
                        let g = dot0 * dot1 / dist2;
                        let shadow_ray = Ray::new_with_offset(hit_pos, ray_dir, &surface_normal);
                        if let Some((_, id2)) = scene.intersect(&shadow_ray) {
                            if id2 == LIGHT_ID {
                                // 光源の寄与
                                let col = sphere.emission
                                    + glm::matrix_comp_mult(
                                        &scene.spheres[LIGHT_ID].emission,
                                        &sphere.color,
                                    ) * std::f64::consts::FRAC_1_PI
                                        * g
                                        / (1.0
                                            / (4.0
                                                * PI
                                                * scene.spheres[LIGHT_ID].radius
                                                * scene.spheres[LIGHT_ID].radius));
                                accumulated_flux += col;
                            }
                        }
                    }
                }
                return accumulated_flux / 128.0;
            }
            Reflection::Refraction => {
                if depth > 5 {
                    return glm::zero::<glm::DVec3>();
                }
                match calc_refraction_behavior(&ray.dir, &geometry_normal, &surface_normal) {
                    Refraction::TotalReflection(dir) => {
                        let ray = Ray::new_with_offset(hit_pos, dir, &surface_normal);
                        return sphere.emission
                            + glm::matrix_comp_mult(
                                &sphere.color,
                                &estimate_radiance_direct(scene, &ray, depth + 1, rand),
                            );
                    }
                    Refraction::ReflectionOrTransmission(
                        (reflect_dir, rp),
                        (refract_dir, tp),
                        probability,
                    ) => {
                        if depth > 2 {
                            if rand.gen::<f64>() < probability {
                                let ray =
                                    Ray::new_with_offset(hit_pos, reflect_dir, &surface_normal);
                                let rad = glm::matrix_comp_mult(
                                    &sphere.color,
                                    &(estimate_radiance_direct(scene, &ray, depth + 1, rand) * rp),
                                ) / probability;
                                return sphere.emission + rad;
                            } else {
                                let ray =
                                    Ray::new_with_offset(hit_pos, refract_dir, &(-surface_normal));
                                let rad = glm::matrix_comp_mult(
                                    &sphere.color,
                                    &(estimate_radiance_direct(scene, &ray, depth + 1, rand) * tp),
                                ) / (1.0 - probability);
                                return sphere.emission + rad;
                            }
                        } else {
                            let ray_re =
                                Ray::new_with_offset(hit_pos, reflect_dir, &surface_normal);
                            let ray_tr =
                                Ray::new_with_offset(hit_pos, refract_dir, &(-surface_normal));
                            let rad = estimate_radiance_direct(scene, &ray_re, depth + 1, rand)
                                * rp
                                + estimate_radiance_direct(scene, &ray_tr, depth + 1, rand) * tp;
                            return sphere.emission + glm::matrix_comp_mult(&sphere.color, &rad);
                        }
                    }
                }
            }
        }
    }
    glm::zero()
}

#[allow(dead_code)]
pub fn direct_visalization_photonmap(
    scene: &Scene,
    ray: &Ray,
    photon_map: &PhotonMap,
) -> glm::DVec3 {
    if let Some((dist, idx)) = scene.intersect(ray) {
        let sphere = &scene.spheres[idx];
        let hit_pos = ray.pos + ray.dir * dist;
        let geometry_normal = (hit_pos - sphere.pos).normalize();
        let surface_normal = if geometry_normal.dot(&ray.dir) < 0.0 {
            geometry_normal
        } else {
            -geometry_normal
        };
        let arg = PhotonCollectArg::new(0.1, 1, hit_pos, surface_normal);
        let nearest_photons = photon_map.collect_photons(&arg);
        if !nearest_photons.queue.is_empty() {
            let mut accum = glm::zero::<glm::DVec3>();
            for n in nearest_photons.queue {
                accum += n.photon.power;
            }
            return accum;
        }
    }
    glm::zero()
}

enum Refraction {
    TotalReflection(glm::DVec3),
    ReflectionOrTransmission((glm::DVec3, f64), (glm::DVec3, f64), f64),
}

fn calc_refraction_behavior(
    incoming_dir: &glm::DVec3,
    geometry_normal: &glm::DVec3,
    surface_normal: &glm::DVec3,
) -> Refraction {
    let reflection_dir =
        (incoming_dir - geometry_normal * 2.0 * geometry_normal.dot(&incoming_dir)).normalize();
    let into = geometry_normal.dot(&surface_normal).signum();
    let nr = REFRACTION_ATM;
    let ni = REFRACTION_DIA;

    let nnt = if into > 0.0 { nr / ni } else { ni / nr };

    let ddn = incoming_dir.dot(&surface_normal);
    let cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);
    if cos2t < 0.0 {
        return Refraction::TotalReflection(reflection_dir);
    }
    let refract_dir =
        (incoming_dir * nnt - geometry_normal * into * (ddn * nnt + cos2t.sqrt())).normalize();
    let r0 = {
        let (a, b) = (ni - nr, ni + nr);
        (a * a) / (b * b)
    };
    let c = if into > 0.0 {
        1.0 + ddn
    } else {
        1.0 - refract_dir.dot(&geometry_normal)
    };
    let re = r0 + (1.0 - r0) * c * c * c * c * c;
    let probability = 0.25 + 0.5 * re;
    let tr = 1.0 - re;
    Refraction::ReflectionOrTransmission((reflection_dir, re), (refract_dir, tr), probability)
}
