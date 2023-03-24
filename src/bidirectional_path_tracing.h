//
// Created by creeper on 23-3-24.
//
#pragma once
#include "Subpath.h"

void generateLightSubpath(const Scene& scene, pcg32_state &rng, Subpath& subpath)
{
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    PointAndNormal light_PN = sample_point_on_light(scene.lights[light_id], scene.shapes[scene.lights[light_id].shape_id]);
}

void generateEyeSubpath(const Scene& scene, int x, int y, pcg32_state &rng, Subpath& subpath)
{
    Real X = (static_cast<Real>(x) + next_pcg32_real<Real>(rng)) / scene.camera.width;
    Real Y = (static_cast<Real>(y) + next_pcg32_real<Real>(rng))/ scene.camera.height;
    Ray ray = sample_primary(scene.camera, Vector2(X, Y));
    Real PRR = 0.95;
    for(int bounce = 0; scene.options.max_depth >= 0 && bounce < scene.options.max_depth; bounce++)
    {
        if(bounce > scene.options.rr_depth && next_pcg32_real<Real>(rng) > PRR) return;
        std::optional<PathVertex> opt_pv = intersect(scene, ray);
        if(opt_pv)
        {
            PathVertex pathVertex = opt_pv.value();
            subpath.append(pathVertex);
            if(is_light(scene.shapes[pathVertex.shape_id])) return;
        }
        else return;
    }
}

Spectrum bidirectional_path_tracing(const Scene& scene,
                                    int x, int y,
                                    pcg32_state &rng) {

}