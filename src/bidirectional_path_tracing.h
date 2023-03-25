//
// Created by creeper on 23-3-24.
//
#pragma once
#include "Subpath.h"

Vector3 sample_uniformly_on_hemisphere(const Vector3& normal, const Vector2 &uv)
{
    Real phi = M_2_PI * uv[0];
    Real theta = acos(uv[1]);
    Vector3 ret(cos(theta) * cos(phi), cos(theta) * sin(phi), sin(theta));
    Frame frame(normal);
    return to_world(frame, ret);
}

void generateLightSubpath(const Scene& scene, pcg32_state &rng, Real PRR, Subpath& subpath)
{
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    PointAndNormal light_PN = sample_point_on_light_surface(scene.lights[light_id],
                                                            {next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)},
                                                            next_pcg32_real<Real>(rng), scene);
    if(is_envmap(scene.lights[light_id]))
    {
        std::cerr << "Sorry, bidirectional path tracing for environment mapping is not supported yet." << std::endl;
        std::cerr << "Rendering aborted." << std::endl;
        exit(0);
    }
    Ray ray;
    ray.org = light_PN.position;
    Vector2 dir_uv = {next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    ray.dir = sample_uniformly_on_hemisphere(light_PN.normal, dir_uv);
    ray.tnear = get_intersection_epsilon(scene);
    ray.tfar = infinity<Real>();
    PathVertex pv_light;
    pv_light.position = light_PN.position;
    pv_light.geometry_normal = light_PN.normal;
    subpath.append(pv_light);
    for(int bounce = 0; scene.options.max_depth < 0 || bounce < scene.options.max_depth; bounce++)
    {
        if(bounce > scene.options.rr_depth && next_pcg32_real<Real>(rng) > PRR) return;
        std::optional<PathVertex> opt_pv = intersect(scene, ray);
        if(opt_pv)
        {
            PathVertex pathVertex = opt_pv.value();
            if(is_light(scene.shapes[pathVertex.shape_id]) continue;
            subpath.append(pathVertex);
            ray.org = pathVertex.position;
        }
        else return;
    }
}

void generateEyeSubpath(const Scene& scene, int x, int y, pcg32_state &rng, Real PRR, Subpath& subpath)
{
    Real X = (static_cast<Real>(x) + next_pcg32_real<Real>(rng)) / scene.camera.width;
    Real Y = (static_cast<Real>(y) + next_pcg32_real<Real>(rng))/ scene.camera.height;
    Ray ray = sample_primary(scene.camera, Vector2(X, Y));
    PathVertex pv_camera;
    pv_camera.exterior_medium_id = scene.camera.medium_id;
    pv_camera.position = {0, 0, 0};
    subpath.append(pv_camera);
    Real PRR = 0.95;
    for(int bounce = 0; scene.options.max_depth < 0 || bounce < scene.options.max_depth; bounce++)
    {
        if(bounce > scene.options.rr_depth && next_pcg32_real<Real>(rng) > PRR) return;
        std::optional<PathVertex> opt_pv = intersect(scene, ray);
        if(opt_pv)
        {
            PathVertex pathVertex = opt_pv.value();
            subpath.append(pathVertex);
            if(is_light(scene.shapes[pathVertex.shape_id])) return;
            ray.org = pathVertex.position;
            const Material &mat = scene.materials[pathVertex.material_id];
            auto bRec = sample_bsdf(mat, -ray.dir, pathVertex, scene.texture_pool,
                        {next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)},
                        next_pcg32_real<Real>(rng), TransportDirection::TO_LIGHT);
        }
        else return;
    }
}

Spectrum bidirectional_path_tracing(const Scene& scene,
                                    int x, int y,
                                    pcg32_state &rng) {
    Subpath eyePath, lightPath;
    generateEyeSubpath(scene, x, y, rng, eyePath);
    generateLightSubpath(scene, rng, lightPath);
    uint nE = eyePath.length(), nL = lightPath.length();
    for(uint s = 0; s < nE; s++)
    {

    }
    for(uint t = 0; t < nL; t++)
    {

    }
    for(uint s = 0; s < ; s++)
    {
        for(uint t = 0; t < ; t++)
        {

        }
    }
}