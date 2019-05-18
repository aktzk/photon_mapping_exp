extern crate nalgebra_glm as glm;
mod photon_map;
mod ray;
mod scene;
mod sphere;
mod trace;
mod utility;
use indicatif::{HumanDuration};
use rayon::prelude::*;
use std::time::Instant;
use trace::*;
use utility::*;

const BYTES_PER_PIXEL: usize = 3;
const WIDTH: usize = 700;
const HEIGHT: usize = 500;
const PHOTON_NUM_CAUSTICS: usize = 500_00;
const PHOTON_NUM_INDIRECT: usize = 100_00;
fn main() {
    let camera_pos = glm::Vec3::new(50.0, 52.0, 325.6);
    let camera_dir = glm::Vec3::new(0.0, -0.05612, -0.5).normalize();
    let elem = WIDTH as f32 * 0.5135 / HEIGHT as f32;
    let cx = glm::Vec3::new(elem, 0.0, 0.0);
    let cy = cx.cross(&camera_dir).normalize() * 0.5135;
    let mut rand = rand::thread_rng();

    let started = Instant::now();

    let scene = scene::Scene::new();
    let photon_map_caustics = Box::new(create_photon_map_caustics(
        PHOTON_NUM_CAUSTICS,
        &scene,
        &mut rand,
    ));

    let photon_map_indirect = Box::new(create_photon_map_indirect(
        PHOTON_NUM_INDIRECT,
        &scene,
        &mut rand,
    ));

    let mut bytes_caustics = vec![0f32; WIDTH as usize * HEIGHT as usize * BYTES_PER_PIXEL];
    let mut bytes_indirect = vec![0f32; WIDTH as usize * HEIGHT as usize * BYTES_PER_PIXEL];
    let mut bytes_direct = vec![0f32; WIDTH as usize * HEIGHT as usize * BYTES_PER_PIXEL];

    let timer = Instant::now();
    println!("Photon Tracing: Caustics");
    bytes_caustics
        .par_chunks_mut(BYTES_PER_PIXEL)
        .into_par_iter()
        .rev()
        .enumerate()
        .for_each(|(idx, chunk)| {
            let mut rng = rand::thread_rng();
            let y = idx / WIDTH as usize;
            let x = idx % WIDTH as usize;
            let mut radiance = glm::zero::<glm::Vec3>();
            for sx in 0..2 {
                for sy in 0..2 {
                    let p_ray = primary_ray(
                        &camera_pos,
                        &camera_dir,
                        (WIDTH, HEIGHT),
                        (x, y),
                        (sx, sy),
                        (&cx, &cy),
                        &mut rng,
                    );

                    radiance +=
                        estimate_radiance_caustics(&scene, &p_ray, &photon_map_caustics, 256.0, 64);
                }
            }
            chunk[0] = radiance.x as f32;
            chunk[1] = radiance.y as f32;
            chunk[2] = radiance.z as f32;
        });
    println!("Done. Took {}", HumanDuration(timer.elapsed()));

    let timer = Instant::now();
    println!("Photon Tracing: Indirect");
    bytes_indirect
        .par_chunks_mut(BYTES_PER_PIXEL)
        .into_par_iter()
        .rev()
        .enumerate()
        .for_each(|(idx, chunk)| {
            let mut rng = rand::thread_rng();
            let y = idx / WIDTH as usize;
            let x = idx % WIDTH as usize;
            let mut radiance = glm::zero::<glm::Vec3>();
            for sx in 0..2 {
                for sy in 0..2 {
                    let p_ray = primary_ray(
                        &camera_pos,
                        &camera_dir,
                        (WIDTH, HEIGHT),
                        (x, y),
                        (sx, sy),
                        (&cx, &cy),
                        &mut rng,
                    );

                    radiance += estimate_radiance_indirect(
                        &scene,
                        &p_ray,
                        &photon_map_indirect,
                        0,
                        128.0,
                        32,
                        &mut rng,
                    );
                }
            }
            chunk[0] = radiance.x as f32;
            chunk[1] = radiance.y as f32;
            chunk[2] = radiance.z as f32;
        });
    println!("Done. Took {}", HumanDuration(timer.elapsed()));

    let timer = Instant::now();
    println!("Ray Tracing: Direct");
    bytes_direct
        .par_chunks_mut(BYTES_PER_PIXEL)
        .into_par_iter()
        .rev()
        .enumerate()
        .for_each(|(idx, chunk)| {
            let mut rng = rand::thread_rng();
            let y = idx / WIDTH as usize;
            let x = idx % WIDTH as usize;
            let mut radiance = glm::zero::<glm::Vec3>();
            for sx in 0..2 {
                for sy in 0..2 {
                    let p_ray = primary_ray(
                        &camera_pos,
                        &camera_dir,
                        (WIDTH, HEIGHT),
                        (x, y),
                        (sx, sy),
                        (&cx, &cy),
                        &mut rng,
                    );

                    radiance += estimate_radiance_direct(&scene, &p_ray, 0, &mut rng);
                }
            }
            chunk[0] = radiance.x as f32;
            chunk[1] = radiance.y as f32;
            chunk[2] = radiance.z as f32;
        });
    println!("Done. Took {}", HumanDuration(timer.elapsed()));

    {
        let mut imagebuf_caustics = image::ImageBuffer::new(WIDTH as u32, HEIGHT as u32);
        let mut idx = 0;
        let mut chunks_caustics = bytes_caustics.chunks(BYTES_PER_PIXEL);
        while let Some([r, g, b]) = chunks_caustics.next() {
            let y = idx / WIDTH as u32;
            let x = idx % WIDTH as u32;

            // powf(0.4545)はガンマ補正です(^(1.0 / 2.2))
            let red = ((r).powf(0.4545).min(1.0).max(0.0) * 255.0) as u8;
            let green = ((g).powf(0.4545).min(1.0).max(0.0) * 255.0) as u8;
            let blue = ((b).powf(0.4545).min(1.0).max(0.0) * 255.0) as u8;
            let pixel = image::Rgb([red, green, blue]);
            imagebuf_caustics.put_pixel(x, y, pixel);
            idx += 1;
        }
        imagebuf_caustics.save("result_caustics.png").unwrap();
    }
    {
        let mut imagebuf_indirect = image::ImageBuffer::new(WIDTH as u32, HEIGHT as u32);
        let mut idx = 0;
        let mut chunks_indirect = bytes_indirect.chunks(BYTES_PER_PIXEL);
        while let Some([r, g, b]) = chunks_indirect.next() {
            let y = idx / WIDTH as u32;
            let x = idx % WIDTH as u32;

            let red = ((r).powf(0.4545).min(1.0).max(0.0) * 255.0) as u8;
            let green = ((g).powf(0.4545).min(1.0).max(0.0) * 255.0) as u8;
            let blue = ((b).powf(0.4545).min(1.0).max(0.0) * 255.0) as u8;
            let pixel = image::Rgb([red, green, blue]);
            imagebuf_indirect.put_pixel(x, y, pixel);
            idx += 1;
        }
        imagebuf_indirect.save("result_indirect.png").unwrap();
    }
    {
        let mut imagebuf_direct = image::ImageBuffer::new(WIDTH as u32, HEIGHT as u32);
        let mut idx = 0;
        let mut chunks_direct = bytes_direct.chunks(BYTES_PER_PIXEL);
        while let Some([r, g, b]) = chunks_direct.next() {
            let y = idx / WIDTH as u32;
            let x = idx % WIDTH as u32;

            let red = ((r).powf(0.4545).min(1.0).max(0.0) * 255.0) as u8;
            let green = ((g).powf(0.4545).min(1.0).max(0.0) * 255.0) as u8;
            let blue = ((b).powf(0.4545).min(1.0).max(0.0) * 255.0) as u8;
            let pixel = image::Rgb([red, green, blue]);
            imagebuf_direct.put_pixel(x, y, pixel);
            idx += 1;
        }
        imagebuf_direct.save("result_direct.png").unwrap();
    }
    {
        let mut imagebuf_sum = image::ImageBuffer::new(WIDTH as u32, HEIGHT as u32);
        let mut idx = 0;
        let chunks_direct = bytes_direct.chunks(BYTES_PER_PIXEL);
        let chunks_indirect = bytes_indirect.chunks(BYTES_PER_PIXEL);
        let chunks_caustics = bytes_caustics.chunks(BYTES_PER_PIXEL);
        let mut zipped = chunks_direct.zip(chunks_indirect).zip(chunks_caustics);
        while let Some(((&[dr, dg, db], &[ir, ig, ib]), &[cr, cg, cb])) = zipped.next() {
            let y = idx / WIDTH as u32;
            let x = idx % WIDTH as u32;

            let red = ((dr + ir + cr).powf(0.4545).min(1.0).max(0.0) * 255.0) as u8;
            let green = ((dg + ig + cg).powf(0.4545).min(1.0).max(0.0) * 255.0) as u8;
            let blue = ((db + ib + cb).powf(0.4545).min(1.0).max(0.0) * 255.0) as u8;
            let pixel = image::Rgb([red, green, blue]);
            imagebuf_sum.put_pixel(x, y, pixel);
            idx += 1;
        }
        imagebuf_sum.save("result.png").unwrap();
    }

    println!("Took {}", HumanDuration(started.elapsed()));
}
