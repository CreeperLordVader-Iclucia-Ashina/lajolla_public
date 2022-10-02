#pragma once
#include "ray.h"
#include "volume.h"
#include "intersection.h"
// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng)
{
    // Homework 2: implememt this!
    Real X = static_cast<Real>(x) / scene.camera.width;
    Real Y = static_cast<Real>(y) / scene.camera.height;
    const Ray ray = sample_primary(scene.camera, Vector2(X, Y));
    std::optional<PathVertex> opt_path_vertex = intersect(scene, ray);
    Spectrum L(0.0, 0.0, 0.0);
    if(opt_path_vertex)
    {
        PathVertex path_vertex = opt_path_vertex.value();
        Real t = distance(ray.org, path_vertex.position);
        Spectrum Le = emission(path_vertex, -ray.dir, scene);
        if(path_vertex.exterior_medium_id >= 0)
        {
          const Medium in_medium = scene.media[path_vertex.exterior_medium_id];
          L += Le * exp(-get_sigma_a(in_medium, path_vertex.position) * t);
        }
    }
    return L;
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum L_scatter(const Scene& scene, const Vector3& pos, const Vector3& dir_out, pcg32_state& rng)
{
    Real u = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, u);
    Vector2 uv(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
    u = next_pcg32_real<Real>(rng);
    PointAndNormal light_src = sample_point_on_light(scene.lights[light_id], pos, uv, u, scene);
    Ray ray;
    ray.org = pos;
    ray.dir = normalize(light_src.position - pos);
    ray.tnear = Real(0);
    ray.tfar = infinity<Real>();
    std::optional<PathVertex> opt_path_vertex = intersect(scene, ray);
    Spectrum L(0.0, 0.0, 0.0);
    if(opt_path_vertex)
    {
        PathVertex path_vertex = opt_path_vertex.value();
        Real t = distance(pos, path_vertex.position);
        Spectrum Le = emission(path_vertex, pos - light_src.position, scene);//emission() has already considered visibility
        if(path_vertex.exterior_medium_id >= 0)
        {
            Medium medium = scene.media[path_vertex.exterior_medium_id];
            PhaseFunction phase = get_phase_function(medium);
            Real sigma_t = get_sigma_a(medium, path_vertex.position).x
                          + get_sigma_s(medium, path_vertex.position).x;
            Spectrum rho = eval(phase, ray.dir, dir_out);
            Spectrum transmit = rho * Le * exp(-sigma_t * t)
                              * abs(dot(ray.dir, light_src.normal)) / (t * t);
            Real pdf = pdf_point_on_light(scene.lights[light_id], light_src, pos, scene)
                            * light_pmf(scene, light_id);
            L += transmit / pdf;
        }
    }
    return L;
}
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    Real X = static_cast<Real>(x) / scene.camera.width;
    Real Y = static_cast<Real>(y) / scene.camera.height;
    const Ray ray = sample_primary(scene.camera, Vector2(X, Y));
    std::optional<PathVertex> opt_path_vertex = intersect(scene, ray);
    Spectrum L(0.0, 0.0, 0.0);
    if(opt_path_vertex)
    {
        PathVertex path_vertex = opt_path_vertex.value();
        Real t_hit = distance(ray.org, path_vertex.position);
        Real u = next_pcg32_real<Real>(rng);
        if(path_vertex.exterior_medium_id >= 0)
        {
            const Medium ext_medium = scene.media[path_vertex.exterior_medium_id];
            Real sigma_s = get_sigma_s(ext_medium, path_vertex.position).x;
            Real sigma_t = get_sigma_a(ext_medium, path_vertex.position).x + sigma_s;
            Real t = -log(1.0 - u) / sigma_t;
            if(t < t_hit)
            {
                Real pdf = sigma_t * exp(-sigma_t * t);
                Spectrum L_sc = L_scatter(scene, ray.org + t * ray.dir, -ray.dir, rng);
                L += sigma_s * exp(-sigma_t * t) * L_sc / pdf;
            }
            else
            {
                Spectrum Le = emission(path_vertex, -ray.dir, scene);
                L += Le;
            }
        }
    }
    else
    {
        Real u = next_pcg32_real<Real>(rng);
        if(scene.camera.medium_id >= 0)
        {
            const Medium& medium = scene.media[scene.camera.medium_id];
            Real sigma_s = get_sigma_s(medium, ray.org).x;
            Real sigma_t = get_sigma_a(medium, ray.org).x + sigma_s;
            Real t = -log(1.0 - u) / sigma_t;
            Real pdf = sigma_t * exp(-sigma_t * t);
            Spectrum L_sc = L_scatter(scene, ray.org + t * ray.dir, -ray.dir, rng);
            L += sigma_s * exp(-sigma_t * t) * L_sc / pdf;
        }
    }
    return L;
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}
