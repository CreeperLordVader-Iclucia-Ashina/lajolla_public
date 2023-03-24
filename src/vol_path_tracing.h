#pragma once

#include <assert.h>
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
    // Homework 2: implement this!
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

// The second-simplest volumetric renderer:
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
void update_medium(const PathVertex& pathVertex, const Ray& ray, int& cur_medium_id)
{
    if(pathVertex.exterior_medium_id != pathVertex.interior_medium_id)
    {
        if (dot(ray.dir, pathVertex.geometry_normal) >= 0)
            cur_medium_id = pathVertex.exterior_medium_id;
        else cur_medium_id = pathVertex.interior_medium_id;
    }
}
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng)
{
    // Homework 2: implement this!
    Spectrum cur_path_throughput(1.0, 1.0, 1.0);
    Real X = static_cast<Real>(x) / scene.camera.width;
    Real Y = static_cast<Real>(y) / scene.camera.height;
    Ray ray = sample_primary(scene.camera, Vector2(X, Y));
    int cur_medium_id = scene.camera.medium_id, depth = 0;
    Spectrum L(0.0, 0.0, 0.0);
    while(true)
    {
        if(scene.options.max_depth >= 0 && depth >= scene.options.max_depth)break;
        if(scene.options.rr_depth >= 0 && depth >= scene.options.rr_depth)//first do a Russian Roulette
        {
            Real P_RR = min(cur_path_throughput.x, 0.95);
            Real rnd_rr = next_pcg32_real<Real>(rng);
            if(P_RR < rnd_rr) return L / P_RR;
        }
        std::optional<PathVertex> opt_path_vertex = intersect(scene, ray);
        if(cur_medium_id >= 0)// if we are now in a medium
        {
            Real sigma_s = get_sigma_s(scene.media[cur_medium_id], ray.org).x;
            Real sigma_t = get_sigma_a(scene.media[cur_medium_id], ray.org).x + sigma_s;
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1.0 - u) / sigma_t;
            Real pdf = sigma_t * exp(-sigma_t * t);
            if (opt_path_vertex) //if intersect
            {
                PathVertex path_vertex = opt_path_vertex.value();
                Real t_hit = distance(ray.org, path_vertex.position);
                if (t >= t_hit) // we calc the radiance emitted and continue;
                {
                    if(path_vertex.material_id == -1)//simply go through the interface
                    {
                        depth++;
                        ray.org += t_hit * ray.dir;
                        //cur_path_throughput *= exp(-sigma_t * t_hit) / pdf;
                        update_medium(path_vertex, ray, cur_medium_id);
                        //if(cur_medium_id >= 0)
                        //{
                        //    sigma_t = get_majorant(scene.media[cur_medium_id], ray).x;
                        //    cur_path_throughput *= exp(-sigma_t * (t - t_hit));
                        //}
                        continue;
                    }
                    else
                    {
                        L += cur_path_throughput * emission(path_vertex, -ray.dir, scene);
                        break;
                    }
                }
            }//else we sample a direction and calc the scattered radiance
            cur_path_throughput *= exp(-sigma_t * t) / pdf;
            ray.org += t * ray.dir;//now we go to this new place to calc the radiance scattered there
            PhaseFunction phaseFunction = get_phase_function(scene.media[cur_medium_id]);
            Vector2 rnd_uv(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
            std::optional<Vector3> opt_dir = sample_phase_function(phaseFunction, ray.dir, rnd_uv);
            if (opt_dir)
            {
                Vector3 dir = opt_dir.value();
                Real phase_pdf = pdf_sample_phase(phaseFunction, dir, -ray.dir);
                cur_path_throughput *= sigma_s * eval(phaseFunction, dir, -ray.dir) / phase_pdf;
                ray.dir = dir;
            }
        }
        else
        {
            if(!opt_path_vertex) break;
            else
            {
                PathVertex path_vertex = opt_path_vertex.value();
                L += cur_path_throughput * emission(path_vertex, -ray.dir, scene);
            }
        }
        depth++;
    }
    return L;
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Real Geometry(const PointAndNormal& pointAndNormal, const Vector3& pos)
{
    Vector3 dif = pos - pointAndNormal.position;
    Real dist = length(dif);
    dif /= dist;
    return max(0.0, dot(dif, pointAndNormal.normal)) / (dist * dist);
}
Spectrum next_event_estimation(const Ray& ray, int depth, int cur_medium_id, const Scene& scene, pcg32_state& rng)
{
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    const PointAndNormal& light_PN = sample_point_on_light(scene.lights[light_id], ray.org,
                                                    Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)), next_pcg32_real<Real>(rng), scene);
    Ray nee_ray;
    nee_ray.org = ray.org;
    nee_ray.dir = normalize(light_PN.position - ray.org);
    nee_ray.tnear = get_intersection_epsilon(scene);
    nee_ray.tfar = infinity<Real>();
    const PhaseFunction& phaseFunction = get_phase_function(scene.media[cur_medium_id]);
    Real nee_throughput = eval(phaseFunction, -ray.dir, nee_ray.dir).x;
    int nee_dep = 0;
    Real pdf_nee = pdf_point_on_light(scene.lights[light_id], light_PN, ray.org, scene) * light_pmf(scene, light_id);
    Real p_march = 1.0;//reduce as we pass through different media
    Real G = Geometry(light_PN, ray.org);
    while(true)//we do a next event estimation
    {
        if(scene.options.max_depth >= 0 && depth + nee_dep >= scene.options.max_depth)
            return make_zero_spectrum();
        std::optional<PathVertex> opt_pathVertex = intersect(scene, nee_ray);
        if(opt_pathVertex)
        {
            PathVertex pathVertex = opt_pathVertex.value();
            Real t = distance(pathVertex.position, nee_ray.org);
            if(cur_medium_id >= 0)
            {
                Real sigma_a = get_sigma_a(scene.media[cur_medium_id], nee_ray.org).x;
                Real sigma_s = get_sigma_s(scene.media[cur_medium_id], nee_ray.org).x;
                Real sigma_t = sigma_a + sigma_s;
                nee_throughput *= exp(-sigma_t * t);
                p_march *= exp(-sigma_t * t);
            }
            nee_ray.org = pathVertex.position;
            if(pathVertex.material_id >= 0 && get_area_light_id(scene.shapes[pathVertex.shape_id]) != light_id)
                return make_zero_spectrum();
            if(distance(nee_ray.org, light_PN.position) < get_intersection_epsilon(scene)) // we reach the light source, and we calc the contrib
            {
                Spectrum L_nee = nee_throughput * emission(pathVertex, -nee_ray.dir, scene) * G / pdf_nee;
                Real pdf_phase = pdf_sample_phase(phaseFunction, -nee_ray.dir, ray.dir)
                                    * G * p_march;
                Real w = pdf_nee * pdf_nee / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
                return L_nee * w;
            }
            update_medium(pathVertex, nee_ray, cur_medium_id);
            nee_dep++;
        }
        else return make_zero_spectrum();
    }
}
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implement this!
    Spectrum cur_path_throughput(1.0, 1.0, 1.0);
    Real X = (static_cast<Real>(x) + next_pcg32_real<Real>(rng)) / scene.camera.width;
    Real Y = (static_cast<Real>(y) + next_pcg32_real<Real>(rng))/ scene.camera.height;
    Ray ray = sample_primary(scene.camera, Vector2(X, Y));
    int cur_medium_id = scene.camera.medium_id, depth = 0;
    Spectrum L(0.0, 0.0, 0.0);
    Vector3 pos_nee_cache;
    Real p_march = 1.0;
    Vector3 dir_in_cache;
    int medium_id_cache = -1;
    bool never_scatter = true;
    Real p_path = 1.0;
    while(true)
    {
        if(scene.options.max_depth >= 0 && depth >= scene.options.max_depth)break;
        if(scene.options.rr_depth >= 0 && depth >= scene.options.rr_depth)//first do a Russian Roulette
        {
            Real rnd_rr = next_pcg32_real<Real>(rng);
            Real P_RR = std::min(0.95, cur_path_throughput.x / p_path);
            if(rnd_rr > P_RR) break;
            else cur_path_throughput /= P_RR;
        }
        std::optional<PathVertex> opt_path_vertex = intersect(scene, ray);
        if(cur_medium_id >= 0)// if we are now in a medium
        {
            Real t_hit;
            PathVertex path_vertex;
            if(!opt_path_vertex) t_hit = infinity<Real>();
            else
            {
                path_vertex = opt_path_vertex.value();
                t_hit = distance(ray.org, path_vertex.position);
            }
            Real sigma_s = get_sigma_s(scene.media[cur_medium_id], ray.org).x;
            Real sigma_t = get_sigma_a(scene.media[cur_medium_id], ray.org).x + sigma_s;
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1.0 - u) / sigma_t;
            if (t >= t_hit)
            {
                p_march *= exp(-t_hit * sigma_t);
                p_path *= exp(-t_hit * sigma_t);
                if(path_vertex.material_id == -1)//simply go through the interface
                {
                    depth++;
                    ray.org = path_vertex.position;
                    update_medium(path_vertex, ray, cur_medium_id);
                    continue;
                }
                else // perform MIS
                {
                    PointAndNormal light_PN = {path_vertex.position, path_vertex.geometry_normal};
                    int light_id = get_area_light_id(scene.shapes[path_vertex.shape_id]);
                    if(never_scatter && light_id != -1)
                        L += cur_path_throughput * emission(path_vertex, -ray.dir, scene);
                    else if(medium_id_cache != -1 && light_id != -1)
                    {
                        Real G = Geometry(light_PN, pos_nee_cache);
                        PhaseFunction phaseFunction = get_phase_function(scene.media[medium_id_cache]);
                        Real pdf_dir = pdf_sample_phase(phaseFunction, dir_in_cache, ray.dir)
                                       * p_march * G;
                        Real pdf_nee = pdf_point_on_light(scene.lights[light_id], light_PN,
                                                          pos_nee_cache, scene) * light_pmf(scene, light_id);
                        Real w = pdf_dir * pdf_dir / (pdf_dir * pdf_dir + pdf_nee * pdf_nee);
                        L += cur_path_throughput * emission(path_vertex, -ray.dir, scene) * w;
                    }
                    else break;
                }
            }
            else//we perform nee and phase function sampling
            {
                p_path *= sigma_t;
                cur_path_throughput *= sigma_s / sigma_t;
                ray.org += t * ray.dir;//now we go to this new place to calc the radiance scattered there
                Spectrum L_nee = next_event_estimation(ray, depth, cur_medium_id, scene, rng);
                L += cur_path_throughput * L_nee;
                medium_id_cache = cur_medium_id;
                dir_in_cache = -ray.dir;
                pos_nee_cache = ray.org;
                never_scatter = false;
                p_march = 1.0;
                PhaseFunction phaseFunction = get_phase_function(scene.media[cur_medium_id]);
                Vector2 rnd_uv(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
                std::optional<Vector3> opt_dir = sample_phase_function(phaseFunction, ray.dir, rnd_uv);
                if (opt_dir)
                {
                    Vector3 dir = opt_dir.value();
                    Real phase_pdf = pdf_sample_phase(phaseFunction, dir, -ray.dir);
                    p_path *= phase_pdf;
                    cur_path_throughput *= eval(phaseFunction, dir, -ray.dir) / phase_pdf;
                    ray.dir = dir;
                }
                else break;
            }
        }
        else
        {
            if(!opt_path_vertex) break;
            PathVertex path_vertex = opt_path_vertex.value();
            PointAndNormal light_PN = {path_vertex.position, path_vertex.geometry_normal};
            Real t_hit = distance(ray.org, path_vertex.position);
            int light_id = get_area_light_id(scene.shapes[path_vertex.shape_id]);
            if(light_id >= 0)
            {
                if(never_scatter)
                    L += cur_path_throughput * emission(path_vertex, -ray.dir, scene);
                else if(medium_id_cache != -1)
                {
                    Real G = Geometry(light_PN, pos_nee_cache);
                    PhaseFunction phaseFunction = get_phase_function(scene.media[medium_id_cache]);
                    Real pdf_dir = pdf_sample_phase(phaseFunction, dir_in_cache, ray.dir)
                                   * p_march * G;
                    Real pdf_nee = pdf_point_on_light(scene.lights[light_id], light_PN,
                                                      pos_nee_cache, scene) * light_pmf(scene, light_id);
                    Real w = pdf_dir * pdf_dir / (pdf_dir * pdf_dir + pdf_nee * pdf_nee);
                    L += cur_path_throughput * emission(path_vertex, -ray.dir, scene) * w;
                }
                else break;
            }
            ray.org += t_hit * ray.dir;
            update_medium(path_vertex, ray, cur_medium_id);
        }
        depth++;
    }
    return L;
}

Spectrum next_event_estimation_bsdf(const Ray& ray, int depth, const PathVertex& path_vertex, int cur_medium_id,
                                    const Scene& scene, pcg32_state& rng)
{
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    PointAndNormal light_PN = sample_point_on_light(scene.lights[light_id], ray.org,
                                                    Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)),
                                                    next_pcg32_real<Real>(rng), scene);
    Ray nee_ray;
    nee_ray.org = ray.org;

    nee_ray.dir = normalize(light_PN.position - ray.org);
    nee_ray.tnear = get_intersection_epsilon(scene);
    nee_ray.tfar = infinity<Real>();
    const Material& mat = scene.materials[path_vertex.material_id];
    Spectrum nee_throughput = eval(mat, -ray.dir, nee_ray.dir, path_vertex, scene.texture_pool);
    int nee_dep = 0;
    Real pdf_nee = pdf_point_on_light(scene.lights[light_id], light_PN, ray.org, scene) * light_pmf(scene, light_id);
    Real p_march = 1.0;//reduce as we pass through different media
    Real G = Geometry(light_PN, ray.org);
    update_medium(path_vertex, nee_ray, cur_medium_id);
    while(true)//we do a next event estimation
    {
        if(scene.options.max_depth >= 0 && depth + nee_dep >= scene.options.max_depth)
            return make_zero_spectrum();
        std::optional<PathVertex> opt_pathVertex = intersect(scene, nee_ray);
        if(opt_pathVertex)
        {
            PathVertex pathVertex = opt_pathVertex.value();
            Real t = distance(pathVertex.position, nee_ray.org);
            if(cur_medium_id >= 0)
            {
                Spectrum sigma_a = get_sigma_a(scene.media[cur_medium_id], nee_ray.org);
                Spectrum sigma_s = get_sigma_s(scene.media[cur_medium_id], nee_ray.org);
                Spectrum sigma_t = sigma_a + sigma_s;
                nee_throughput *= exp(-sigma_t.x * t);
                p_march *= exp(-sigma_t.x * t);
            }
            if(pathVertex.material_id >= 0 && get_area_light_id(scene.shapes[pathVertex.shape_id]) != light_id)
                return make_zero_spectrum();
            nee_ray.org = pathVertex.position;
            if(distance(nee_ray.org, light_PN.position) < get_intersection_epsilon(scene)) // we reach the light source, and we calc the contrib
            {
                Spectrum L_nee = nee_throughput * emission(pathVertex, -nee_ray.dir, scene) * G / pdf_nee;
                Real pdf_bsdf = pdf_sample_bsdf(mat, -ray.dir, nee_ray.dir, path_vertex, scene.texture_pool)
                                    * G * p_march;
                Real w = pdf_nee * pdf_nee / (pdf_nee * pdf_nee + pdf_bsdf * pdf_bsdf);
                return L_nee * w;
            }
            update_medium(pathVertex, nee_ray, cur_medium_id);
            nee_dep++;
        }
        else return make_zero_spectrum();
    }
}
// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implement this!
    Spectrum cur_path_throughput(1.0, 1.0, 1.0);
    Real X = (static_cast<Real>(x) + next_pcg32_real<Real>(rng)) / scene.camera.width;
    Real Y = (static_cast<Real>(y) + next_pcg32_real<Real>(rng))/ scene.camera.height;
    Ray ray = sample_primary(scene.camera, Vector2(X, Y));
    ray.tnear = get_intersection_epsilon(scene);
    int cur_medium_id = scene.camera.medium_id, depth = 0;
    Spectrum L(0.0, 0.0, 0.0);
    Vector3 pos_nee_cache;
    Real p_march = 1.0;
    Vector3 dir_in_cache;
    int medium_id_cache = -1;
    PathVertex pv_nee_cache;
    pv_nee_cache.material_id = -1;
    bool never_scatter_or_reflect = true;
    Real p_path = 1.0;
    while(true)
    {
        if(scene.options.max_depth >= 0 && depth >= scene.options.max_depth)break;
        if(scene.options.rr_depth >= 0 && depth >= scene.options.rr_depth)//first do a Russian Roulette
        {
            Real rnd_rr = next_pcg32_real<Real>(rng);
            Real P_RR;
            if(0.95 * p_path <= cur_path_throughput.x) P_RR = 0.95;
            else P_RR = cur_path_throughput.x / p_path;
            if(rnd_rr > P_RR) break;
            else cur_path_throughput /= P_RR;
        }
        std::optional<PathVertex> opt_path_vertex = intersect(scene, ray);
        if(cur_medium_id >= 0)// if we are now in a medium
        {
            Real t_hit;
            PathVertex path_vertex;
            if(!opt_path_vertex) t_hit = infinity<Real>();
            else
            {
                path_vertex = opt_path_vertex.value();
                t_hit = distance(ray.org, path_vertex.position);
            }
            Spectrum sigma_s = get_sigma_s(scene.media[cur_medium_id], ray.org);
            Spectrum sigma_t = get_sigma_a(scene.media[cur_medium_id], ray.org) + sigma_s;
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1.0 - u) / sigma_t.x;
            if (t >= t_hit)
            {
                p_march *= exp(-t_hit * sigma_t.x);
                p_path *= exp(-t_hit * sigma_t.x);
                if(path_vertex.material_id == -1)//it is a index-matching surface, we go through the interface
                {
                    depth++;
                    ray.org = path_vertex.position;
                    update_medium(path_vertex, ray, cur_medium_id);
                    continue;
                }
                else // if it is a light then perform MIS, else perform surface lighting
                {
                    PointAndNormal surface_PN = {path_vertex.position, path_vertex.geometry_normal};
                    int light_id = get_area_light_id(scene.shapes[path_vertex.shape_id]);
                    if(light_id != -1)
                    {
                        if(never_scatter_or_reflect) L += cur_path_throughput * emission(path_vertex, -ray.dir, scene);
                        else
                        {
                            Real pdf_dir, pdf_nee;
                            if(medium_id_cache != -1)// perform MIS for phase function sampling in scattering
                            {
                                Real G = Geometry(surface_PN, pos_nee_cache);
                                PhaseFunction phaseFunction = get_phase_function(scene.media[medium_id_cache]);
                                pdf_dir = pdf_sample_phase(phaseFunction, dir_in_cache, ray.dir)
                                               * p_march * G;
                                pdf_nee = pdf_point_on_light(scene.lights[light_id], surface_PN,
                                                                  pos_nee_cache, scene) * light_pmf(scene, light_id);
                            }
                            else if(pv_nee_cache.material_id != -1)// perform MIS for bsdf sampling in reflecting
                            {
                                Real G = Geometry(surface_PN, pv_nee_cache.position);
                                const Material& mat = scene.materials[pv_nee_cache.material_id];
                                pdf_dir = pdf_sample_bsdf(mat, dir_in_cache, ray.dir, pv_nee_cache, scene.texture_pool)
                                          * p_march * G;
                                pdf_nee = pdf_point_on_light(scene.lights[light_id], surface_PN, pv_nee_cache.position, scene)
                                          * light_pmf(scene, light_id);
                            }
                            else break;
                            Real w = pdf_dir * pdf_dir / (pdf_dir * pdf_dir + pdf_nee * pdf_nee);
                            L += cur_path_throughput * emission(path_vertex, -ray.dir, scene) * w;
                            break;
                        }
                    }
                    else //not light or index-macthing surface, so it must be a solid surface, perform nee
                    {
                        never_scatter_or_reflect = false;
                        ray.org = path_vertex.position;
                        Spectrum L_nee = next_event_estimation_bsdf(ray, depth, path_vertex, cur_medium_id, scene, rng);
                        L += cur_path_throughput * L_nee;
                        medium_id_cache = -1;
                        pv_nee_cache = path_vertex;
                        dir_in_cache = -ray.dir;
                        pos_nee_cache = ray.org;
                        p_march = 1.0;
                        const Material& mat = scene.materials[path_vertex.material_id];
                        Vector2 rnd_uv(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
                        std::optional<BSDFSampleRecord> opt_bRec = sample_bsdf(mat, -ray.dir, path_vertex, scene.texture_pool,
                                                                               rnd_uv, next_pcg32_real<Real>(rng));
                        if (opt_bRec)
                        {
                            BSDFSampleRecord bRec = opt_bRec.value();
                            Real bsdf_pdf = pdf_sample_bsdf(mat, -ray.dir, bRec.dir_out, path_vertex, scene.texture_pool);
                            if(bsdf_pdf == 0.0) break;
                            cur_path_throughput *= eval(mat, -ray.dir, bRec.dir_out, path_vertex, scene.texture_pool) / bsdf_pdf;
                            p_path *= bsdf_pdf;
                            ray.dir = bRec.dir_out;
                            update_medium(path_vertex, ray, cur_medium_id);
                        }
                        else break;
                    }
                }
            }
            else// we perform nee and phase function sampling
            {
                p_path *= sigma_t.x;
                cur_path_throughput *= sigma_s.x / sigma_t.x;
                ray.org += t * ray.dir;//now we go to this new place to calc the radiance scattered there
                Spectrum L_nee = next_event_estimation(ray, depth, cur_medium_id, scene, rng);
                L += cur_path_throughput * L_nee;
                medium_id_cache = cur_medium_id;
                pv_nee_cache.material_id = -1;
                dir_in_cache = -ray.dir;
                pos_nee_cache = ray.org;
                never_scatter_or_reflect = false;
                p_march = 1.0;
                PhaseFunction phaseFunction = get_phase_function(scene.media[cur_medium_id]);
                Vector2 rnd_uv(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
                std::optional<Vector3> opt_dir = sample_phase_function(phaseFunction, -ray.dir, rnd_uv);
                if (opt_dir)
                {
                    Vector3 dir = opt_dir.value();
                    Real phase_pdf = pdf_sample_phase(phaseFunction, -ray.dir, dir);
                    p_path *= phase_pdf;
                    cur_path_throughput *= eval(phaseFunction, -ray.dir, dir) / phase_pdf;
                    ray.dir = dir;
                }
                else break;
            }
        }
        else
        {
            if(!opt_path_vertex) break;
            PathVertex path_vertex = opt_path_vertex.value();
            PointAndNormal surface_PN = {path_vertex.position, path_vertex.geometry_normal};
            int light_id = get_area_light_id(scene.shapes[path_vertex.shape_id]);
            if(light_id >= 0)
            {
                if(never_scatter_or_reflect) L += cur_path_throughput * emission(path_vertex, -ray.dir, scene);
                else
                {
                    Real pdf_dir, pdf_nee;
                    if(medium_id_cache != -1)// perform MIS for phase function sampling in scattering
                    {
                        Real G = Geometry(surface_PN, pos_nee_cache);
                        PhaseFunction phaseFunction = get_phase_function(scene.media[medium_id_cache]);
                        pdf_dir = pdf_sample_phase(phaseFunction, dir_in_cache, ray.dir)
                                    * p_march * G;
                        pdf_nee = pdf_point_on_light(scene.lights[light_id], surface_PN, pos_nee_cache, scene)
                                                        * light_pmf(scene, light_id);
                    }
                    else if(pv_nee_cache.material_id != -1)// perform MIS for bsdf sampling in reflecting
                    {
                        Real G = Geometry(surface_PN, pv_nee_cache.position);
                        const Material& mat = scene.materials[pv_nee_cache.material_id];
                        pdf_dir = pdf_sample_bsdf(mat, dir_in_cache, ray.dir, pv_nee_cache, scene.texture_pool)
                                    * p_march * G;
                        pdf_nee = pdf_point_on_light(scene.lights[light_id], surface_PN, pv_nee_cache.position, scene)
                                    * light_pmf(scene, light_id);
                        //std::cout << dir_in_cache << " " << ray.dir << " " << pv_nee_cache.position << " " << pv_nee_cache.geometry_normal << " " << surface_PN.position << " " << surface_PN.normal << std::endl;
                        //printf("%lf %lf %lf %lf\n", pdf_sample_bsdf(mat, dir_in_cache, ray.dir, pv_nee_cache, scene.texture_pool),
                        //       p_march, G, pdf_nee);
                    }
                    else break;
                    Real w = pdf_dir * pdf_dir / (pdf_dir * pdf_dir + pdf_nee * pdf_nee);
                    L += cur_path_throughput * emission(path_vertex, -ray.dir, scene) * w;
                }
                break;
            }
            else if(path_vertex.material_id >= 0)
            {
                never_scatter_or_reflect = false;
                ray.org = path_vertex.position;
                Spectrum L_nee = next_event_estimation_bsdf(ray, depth, path_vertex, cur_medium_id, scene, rng);
                L += cur_path_throughput * L_nee;
                medium_id_cache = -1;
                pv_nee_cache = path_vertex;
                dir_in_cache = -ray.dir;
                pos_nee_cache = ray.org;
                p_march = 1.0;
                const Material& mat = scene.materials[path_vertex.material_id];
                Vector2 rnd_uv(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
                std::optional<BSDFSampleRecord> opt_bRec = sample_bsdf(mat, -ray.dir, path_vertex, scene.texture_pool,
                                                                       rnd_uv, next_pcg32_real<Real>(rng));
                if (opt_bRec)
                {
                    BSDFSampleRecord bRec = opt_bRec.value();
                    Real bsdf_pdf = pdf_sample_bsdf(mat, -ray.dir, bRec.dir_out, path_vertex, scene.texture_pool);
                    if(bsdf_pdf == 0.0) break;
                    cur_path_throughput *= eval(mat, -ray.dir, bRec.dir_out, path_vertex, scene.texture_pool) / bsdf_pdf;
                    p_path *= bsdf_pdf;
                    ray.dir = bRec.dir_out;
                }
                else break;
            }
            else ray.org = path_vertex.position;
            update_medium(path_vertex, ray, cur_medium_id);
        }
        depth++;
    }
    return L;
}
Real Max(const Spectrum& s)
{
    return std::max(s.x, std::max(s.y, s.z));
}
Real avg(const Spectrum& s)
{
    return (s[0] + s[1] + s[2]) / 3.0;
}
Spectrum nee(const Ray& ray, int depth, int cur_medium_id, const Scene& scene, pcg32_state& rng)
{
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    const PointAndNormal& light_PN = sample_point_on_light(scene.lights[light_id], ray.org,
                                                           Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)), next_pcg32_real<Real>(rng), scene);
    Ray nee_ray;
    nee_ray.org = ray.org;
    nee_ray.dir = normalize(light_PN.position - ray.org);
    nee_ray.tnear = get_intersection_epsilon(scene);
    nee_ray.tfar = infinity<Real>();
    const PhaseFunction& phaseFunction = get_phase_function(scene.media[cur_medium_id]);
    Spectrum nee_throughput = eval(phaseFunction, -ray.dir, nee_ray.dir);
    int nee_dep = 0;
    Real pdf_light = pdf_point_on_light(scene.lights[light_id], light_PN, ray.org, scene) * light_pmf(scene, light_id);
    Real G = Geometry(light_PN, ray.org);
    Real pdf_delta = 1, pdf_ratio = 1.0;
    while(true)//we do a next event estimation
    {
        if(scene.options.max_depth >= 0 && depth + nee_dep >= scene.options.max_depth)
            return make_zero_spectrum();
        std::optional<PathVertex> opt_pathVertex = intersect(scene, nee_ray);
        if(opt_pathVertex)
        {
            PathVertex pathVertex = opt_pathVertex.value();
            if(cur_medium_id >= 0)
            {
                Real t_hit = distance(pathVertex.position, nee_ray.org);
                int channel = floor(next_pcg32_real<Real>(rng) * 3);
                Spectrum sigma_m = get_majorant(scene.media[cur_medium_id], nee_ray);
                int null_collision_cnt = 0;
                Spectrum T = make_const_spectrum(1), trans_pdf_delta = make_const_spectrum(1),
                            trans_pdf_ratio = make_const_spectrum(1);
                while(true)
                {
                    Real u = next_pcg32_real<Real>(rng);
                    if(sigma_m[channel] <= 0.0) break;
                    Real t = -log(1.0 - u) / sigma_m[channel];
                    if (t >= t_hit)
                    {
                        if(pathVertex.material_id == -1)//it is an index-matching surface, we go through the interface
                        {
                            nee_dep++;
                            T *= exp(-sigma_m * t_hit);
                            trans_pdf_delta *= exp(-sigma_m * t_hit);
                            trans_pdf_ratio *= exp(-sigma_m * t_hit);
                        }
                        break;
                    }
                    else// in ratio tracking, we just perform null scattering event
                    {
                        Spectrum sigma_s = get_sigma_s(scene.media[cur_medium_id], nee_ray.org);
                        Spectrum sigma_t = get_sigma_a(scene.media[cur_medium_id], nee_ray.org) + sigma_s;
                        nee_ray.org += t * nee_ray.dir;
                        Spectrum sigma_n = sigma_m - sigma_t;
                        T *= exp(-sigma_m * t) * sigma_n / Max(sigma_m);
                        trans_pdf_delta *= exp(-sigma_m * t) * sigma_n / Max(sigma_m);
                        trans_pdf_ratio *= exp(-sigma_m * t) * sigma_m / Max(sigma_m);
                        null_collision_cnt++;
                        if(null_collision_cnt >= scene.options.max_null_collisions) break;
                        if(Max(T) <= 0.0) return make_zero_spectrum();
                    }
                    t_hit -= t;
                }
                nee_throughput *= (T / avg(trans_pdf_ratio));
                pdf_ratio *= avg(trans_pdf_ratio);
                pdf_delta *= avg(trans_pdf_delta);
            }
            nee_ray.org = pathVertex.position;
            if(pathVertex.material_id >= 0 && get_area_light_id(scene.shapes[pathVertex.shape_id]) != light_id)
                return make_zero_spectrum();
            if(distance(nee_ray.org, light_PN.position) < get_intersection_epsilon(scene)) // we reach the light source, and we calc the contrib
            {
                Real pdf_phase = pdf_sample_phase(phaseFunction, -nee_ray.dir, ray.dir)
                                 * G * pdf_delta;
                Real pdf_nee = pdf_light * pdf_ratio;
                Spectrum L_nee = nee_throughput * emission(pathVertex, -nee_ray.dir, scene) * G / pdf_light;
                Real w = pdf_nee * pdf_nee / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
                return L_nee * w;
            }
            update_medium(pathVertex, nee_ray, cur_medium_id);
            nee_dep++;
        }
        else return make_zero_spectrum();
    }
}
Spectrum nee_bsdf(const Ray& ray, int depth, const PathVertex& path_vertex, int cur_medium_id,
                                    const Scene& scene, pcg32_state& rng)
{
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    const PointAndNormal& light_PN = sample_point_on_light(scene.lights[light_id], ray.org,
                                                           Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)), next_pcg32_real<Real>(rng), scene);
    Ray nee_ray;
    nee_ray.org = ray.org;
    nee_ray.dir = normalize(light_PN.position - ray.org);
    nee_ray.tnear = get_intersection_epsilon(scene);
    nee_ray.tfar = infinity<Real>();
    const Material& mat = scene.materials[path_vertex.material_id];
    Spectrum nee_throughput = eval(mat, -ray.dir, nee_ray.dir, path_vertex, scene.texture_pool);
    int nee_dep = 0;
    Real pdf_light = pdf_point_on_light(scene.lights[light_id], light_PN, ray.org, scene)
                        * light_pmf(scene, light_id);
    Real G = Geometry(light_PN, ray.org);
    Real pdf_delta = 1, pdf_ratio = 1;
    update_medium(path_vertex, nee_ray, cur_medium_id);
    while(true)//we do a next event estimation
    {
        if(scene.options.max_depth >= 0 && depth + nee_dep >= scene.options.max_depth)
            return make_zero_spectrum();
        std::optional<PathVertex> opt_pathVertex = intersect(scene, nee_ray);
        if(opt_pathVertex)
        {
            PathVertex pathVertex = opt_pathVertex.value();
            if(cur_medium_id >= 0)
            {
                Real t_hit = distance(pathVertex.position, nee_ray.org);
                int channel = floor(next_pcg32_real<Real>(rng) * 3);
                Spectrum sigma_m = get_majorant(scene.media[cur_medium_id], nee_ray);
                int null_collision_cnt = 0;
                Spectrum T = make_const_spectrum(1), trans_pdf_delta = make_const_spectrum(1),
                                trans_pdf_ratio = make_const_spectrum(1);
                while(true)
                {
                    Real u = next_pcg32_real<Real>(rng);
                    if(sigma_m[channel] <= 0.0) break;
                    Real t = -log(1.0 - u) / sigma_m[channel];
                    if (t >= t_hit)
                    {
                        if(pathVertex.material_id == -1)//it is a index-matching surface, we go through the interface
                        {
                            depth++;
                            T *= exp(-sigma_m * t_hit);
                            trans_pdf_delta *= exp(-sigma_m * t_hit);
                            trans_pdf_ratio *= exp(-sigma_m * t_hit);
                            nee_ray.org = pathVertex.position;
                        }
                        break;
                    }
                    else// in ratio tracking, we just perform null scattering event
                    {
                        Spectrum sigma_s = get_sigma_s(scene.media[cur_medium_id], nee_ray.org);
                        Spectrum sigma_t = get_sigma_a(scene.media[cur_medium_id], nee_ray.org) + sigma_s;
                        nee_ray.org += t * nee_ray.dir;
                        Spectrum sigma_n = sigma_m - sigma_t;
                        T *= exp(-sigma_m * t) * sigma_n / Max(sigma_m);
                        trans_pdf_delta *= exp(-sigma_m * t) * sigma_n / Max(sigma_m);
                        trans_pdf_ratio *= exp(-sigma_m * t) * sigma_m / Max(sigma_m);
                        null_collision_cnt++;
                        if(null_collision_cnt >= scene.options.max_null_collisions) break;
                        if(Max(T) <= 0.0) return make_zero_spectrum();
                    }
                    t_hit -= t;
                }
                nee_throughput *= (T / avg(trans_pdf_ratio));
                pdf_ratio *= avg(trans_pdf_ratio);
                pdf_delta *= avg(trans_pdf_delta);
            }
            nee_ray.org = pathVertex.position;
            if(pathVertex.material_id >= 0 && get_area_light_id(scene.shapes[pathVertex.shape_id]) != light_id)
                return make_zero_spectrum();
            if(distance(nee_ray.org, light_PN.position) < get_intersection_epsilon(scene)) // we reach the light source, and we calc the contrib
            {
                Real pdf_bsdf = pdf_sample_bsdf(mat, -ray.dir, nee_ray.dir, path_vertex,
                                                     scene.texture_pool) * G * pdf_delta;
                Real pdf_nee = pdf_light * pdf_ratio;
                Spectrum L_nee = nee_throughput * emission(pathVertex, -nee_ray.dir, scene) * G / pdf_light;
                Real w = pdf_nee * pdf_nee / (pdf_nee * pdf_nee + pdf_bsdf * pdf_bsdf);
                return L_nee * w;
            }
            update_medium(pathVertex, nee_ray, cur_medium_id);
            nee_dep++;
        }
        else return make_zero_spectrum();
    }
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implement this!
    Spectrum cur_path_throughput(1.0, 1.0, 1.0);
    Real X = (static_cast<Real>(x) + next_pcg32_real<Real>(rng)) / scene.camera.width;
    Real Y = (static_cast<Real>(y) + next_pcg32_real<Real>(rng))/ scene.camera.height;
    Ray ray = sample_primary(scene.camera, Vector2(X, Y));
    ray.tnear = get_intersection_epsilon(scene);
    int cur_medium_id = scene.camera.medium_id, depth = 0;
    Spectrum L(0.0, 0.0, 0.0);
    Real pdf_ratio = 1, pdf_delta = 1;// p_march for nee
    PathVertex pv_nee_cache;
    Real pdf_dir_cache = 0;
    bool never_scatter_or_reflect = true;
    while(true)
    {
        if(scene.options.max_depth >= 0 && depth >= scene.options.max_depth)break;
        if(scene.options.rr_depth >= 0 && depth >= scene.options.rr_depth)//first do a Russian Roulette
        {
            Real rnd_rr = next_pcg32_real<Real>(rng);
            Real P_RR = 0.95;
            if(rnd_rr > P_RR) break;
            else cur_path_throughput /= P_RR;
        }
        std::optional<PathVertex> opt_path_vertex = intersect(scene, ray);
        bool reach_light = false;
        bool scatter = false;
        bool surface_lighting = false;
        Spectrum sigma_s, sigma_t;
        Real t_hit;
        PathVertex path_vertex;
        PointAndNormal surface_PN;
        int light_id, channel;
        if(!opt_path_vertex) t_hit = infinity<Real>();
        else
        {
            path_vertex = opt_path_vertex.value();
            surface_PN = {path_vertex.position,  path_vertex.geometry_normal};
            t_hit = distance(ray.org, path_vertex.position);
        }
        if(cur_medium_id >= 0)// if we are now in a medium
        {
            channel = floor(next_pcg32_real<Real>(rng) * 3);
            Spectrum sigma_m = get_majorant(scene.media[cur_medium_id], ray);
            int null_collision_cnt = 0;
            Spectrum T = make_const_spectrum(1), trans_pdf_delta = make_const_spectrum(1),
                        trans_pdf_ratio = make_const_spectrum(1);
            while(true)// perform free-flight sampling
            {
                Real u = next_pcg32_real<Real>(rng);
                if(sigma_m[channel] <= 0.0) break;
                Real t = -log(1.0 - u) / sigma_m[channel];
                if (t >= t_hit)
                {
                    trans_pdf_delta *= exp(-t_hit * sigma_m);
                    trans_pdf_ratio *= exp(-t_hit * sigma_m);
                    T *= exp(-t_hit * sigma_m);
                    ray.org = path_vertex.position;
                    if(path_vertex.material_id == -1)//it is a index-matching surface, we go through the interface
                    {
                        depth++;
                        update_medium(path_vertex, ray, cur_medium_id);
                    }
                    else
                    {
                        light_id = get_area_light_id(scene.shapes[path_vertex.shape_id]);
                        if(light_id != -1) reach_light = true;
                        else surface_lighting = true;
                    }
                    break;
                }
                else// we sample to decide whether it is a null collision or scattering event
                {
                    ray.org += t * ray.dir;
                    sigma_s = get_sigma_s(scene.media[cur_medium_id], ray.org);
                    sigma_t = get_sigma_a(scene.media[cur_medium_id], ray.org) + sigma_s;
                    if(next_pcg32_real<Real>(rng) < sigma_t[channel] / sigma_m[channel])
                    {
                        scatter = true;
                        T *= exp(-sigma_m * t) / Max(sigma_m);
                        trans_pdf_delta *= exp(-sigma_m * t) * sigma_t / Max(sigma_m);
                        break;
                    }
                    else
                    {
                        Spectrum sigma_n = sigma_m - sigma_t;
                        T *= exp(-sigma_m * t) * sigma_n / Max(sigma_m);
                        trans_pdf_delta *= exp(-sigma_m * t) * sigma_n / Max(sigma_m);
                        trans_pdf_ratio *= exp(-sigma_m * t) * sigma_m / Max(sigma_m);
                        if(Max(T) <= 0.0) break;
                        null_collision_cnt++;
                        if(null_collision_cnt >= scene.options.max_null_collisions) break;
                    }
                }
                t_hit -= t;
            }
            cur_path_throughput *= (T / avg(trans_pdf_delta));
            pdf_delta *= avg(trans_pdf_delta);
            pdf_ratio *= avg(trans_pdf_ratio);
        }
        else
        {
            if(!opt_path_vertex) break;
            light_id = get_area_light_id(scene.shapes[path_vertex.shape_id]);
            ray.org = path_vertex.position;
            if(light_id >= 0) reach_light = true;
            else if(path_vertex.material_id >= 0) surface_lighting = true;
            else update_medium(path_vertex, ray, cur_medium_id);
        }
        if(scatter)
        {
            cur_path_throughput *= sigma_s;
            Spectrum L_nee = nee(ray, depth, cur_medium_id, scene, rng);
            L += cur_path_throughput * L_nee;
            pv_nee_cache = path_vertex;
            never_scatter_or_reflect = false;
            PhaseFunction phaseFunction = get_phase_function(scene.media[cur_medium_id]);
            Vector2 rnd_uv(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
            std::optional<Vector3> opt_dir = sample_phase_function(phaseFunction, -ray.dir, rnd_uv);
            if (opt_dir)
            {
                Vector3 dir = opt_dir.value();
                Real phase_pdf = pdf_sample_phase(phaseFunction, -ray.dir, dir);
                pdf_dir_cache = phase_pdf;
                pdf_ratio = pdf_delta  = 1.0;
                cur_path_throughput *= eval(phaseFunction, -ray.dir, dir) / phase_pdf;
                ray.dir = dir;
            }
            else break;
        }
        if(reach_light)
        {
            if(never_scatter_or_reflect) L += cur_path_throughput * emission(path_vertex, -ray.dir, scene);
            else
            {
                Real G = Geometry(surface_PN, pv_nee_cache.position);
                Real pdf_dir = pdf_delta * G * pdf_dir_cache;
                Real pdf_nee = pdf_point_on_light(scene.lights[light_id], surface_PN, pv_nee_cache.position,
                                           scene) * light_pmf(scene, light_id) * pdf_ratio;
                Real w = pdf_dir * pdf_dir / (pdf_dir * pdf_dir + pdf_nee * pdf_nee);
                L += cur_path_throughput * emission(path_vertex, -ray.dir, scene) * w;
            }
            break;
        }
        if(surface_lighting)
        {
            never_scatter_or_reflect = false;
            Spectrum L_nee = nee_bsdf(ray, depth, path_vertex, cur_medium_id, scene, rng);
            L += cur_path_throughput * L_nee;
            pv_nee_cache = path_vertex;
            const Material& mat = scene.materials[path_vertex.material_id];
            Vector2 rnd_uv(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
            std::optional<BSDFSampleRecord> opt_bRec = sample_bsdf(mat, -ray.dir, path_vertex, scene.texture_pool,
                                                                   rnd_uv, next_pcg32_real<Real>(rng));
            if (opt_bRec)
            {
                BSDFSampleRecord bRec = opt_bRec.value();
                Real bsdf_pdf = pdf_sample_bsdf(mat, -ray.dir, bRec.dir_out, path_vertex, scene.texture_pool);
                if(bsdf_pdf == 0.0) break;
                cur_path_throughput *= eval(mat, -ray.dir, bRec.dir_out, path_vertex, scene.texture_pool) / bsdf_pdf;
                pdf_dir_cache = bsdf_pdf;
                pdf_delta = pdf_ratio = 1.0;
                ray.dir = bRec.dir_out;
                update_medium(path_vertex, ray, cur_medium_id);
            }
            else break;
        }
        depth++;
    }
    return L;
}
