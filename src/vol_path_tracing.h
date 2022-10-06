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
            cur_medium_id = pathVertex.interior_medium_id;
        else cur_medium_id = pathVertex.exterior_medium_id;
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
                        ray.org += t * ray.dir;
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
    Vector3 dif = pointAndNormal.position - pos;
    return abs(dot(dif, pointAndNormal.normal)) / dot(dif, dif);
}
Spectrum next_event_estimation(const Vector3& pos, const Ray& ray, int depth,
                           int cur_medium_id, const PhaseFunction& phaseFunction, const Scene& scene,
                           pcg32_state& rng)
{
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    PointAndNormal light_PN = sample_point_on_light(scene.lights[light_id],
                                                    ray.org, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)),
                                                    next_pcg32_real<Real>(rng), scene);
    const Vector3& light_N = light_PN.normal;
    Ray nee_ray;
    nee_ray.org = ray.org;
    Real nee_t_hit = distance(ray.org, light_PN.position);
    nee_ray.dir = normalize(light_PN.position - ray.org);
    nee_ray.tnear = 0;
    nee_ray.tfar = infinity<Real>();
    int nee_medium_id = cur_medium_id;
    Real nee_throughput = 1.0;
    int nee_dep = 0;
    Real pdf_nee = pdf_point_on_light(scene.lights[light_id], light_PN, ray.org, scene);
    Real pdf_march = 1.0;//reduce as we pass through different media
    Spectrum L_nee(0.0, 0.0, 0.0);
    Real G = Geometry(light_PN, ray.org);
    while(true)//we do a next event estimation
    {
        if(scene.options.max_depth >= 0 && depth + nee_dep >= scene.options.max_depth)
            return {0.0, 0.0, 0.0};
        PathVertex pathVertex = intersect(scene, nee_ray).value();
        Real t = distance(pathVertex.position, nee_ray.org);
        if(abs(t - nee_t_hit) < 1e-4) // we reach the light source, and we calc the contrib
        {
            L_nee = nee_throughput * emission(pathVertex, -nee_ray.dir, scene) / pdf_nee;
            Real pdf_phase = cur_medium_id >= 0 ?
                             pdf_sample_phase(phaseFunction, -nee_ray.dir, -ray.dir) * G * pdf_march : 0;
            Real w = pdf_nee * pdf_nee / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
            return L_nee * w * eval(phaseFunction, -nee_ray.dir, -ray.dir);
        }
        if(pathVertex.material_id != -1)// we are blocked, so we break
            return {0.0, 0.0, 0.0};
        if(nee_medium_id >= 0)
        {
            Real sigma_t = get_majorant(scene.media[nee_medium_id], nee_ray).x;
            nee_throughput *= exp(-sigma_t * t);
            pdf_march *= exp(-sigma_t * t);
        }
        nee_ray.org += t * nee_ray.dir;
        update_medium(pathVertex, nee_ray, nee_medium_id);
        nee_dep++;
    }
}
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implement this!
    Spectrum cur_path_throughput(1.0, 1.0, 1.0);
    Real X = static_cast<Real>(x) / scene.camera.width;
    Real Y = static_cast<Real>(y) / scene.camera.height;
    Ray ray = sample_primary(scene.camera, Vector2(X, Y));
    int cur_medium_id = scene.camera.medium_id, depth = 0;
    Vector3 p_nee;
    Spectrum L(0.0, 0.0, 0.0);
    Real pdf_dir = 1.0;
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
            PhaseFunction phaseFunction = get_phase_function(scene.media[cur_medium_id]);
            Spectrum L_nee = next_event_estimation(ray.org, ray, depth, cur_medium_id, phaseFunction, scene, rng);
            Real sigma_s = get_sigma_s(scene.media[cur_medium_id], ray.org).x;
            Real sigma_t = get_sigma_a(scene.media[cur_medium_id], ray.org).x + sigma_s;
            L += sigma_s * L_nee;//calc the contribution of L_nee
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1.0 - u) / sigma_t;
            Real pdf = sigma_t * exp(-sigma_t * t);
            Real multi_trans_prob = 1.0;
            if (opt_path_vertex) //if intersect
            {
                PathVertex path_vertex = opt_path_vertex.value();
                Real t_hit = distance(ray.org, path_vertex.position);
                if (t >= t_hit) //we calc the radiance emitted and continue;
                {
                    if(path_vertex.material_id == -1)//simply go through the interface
                    {
                        depth++;
                        ray.org += t_hit * ray.dir;
                        multi_trans_prob *= exp(-sigma_t * t_hit);
                        update_medium(path_vertex, ray, cur_medium_id);
                        continue;
                    }
                    else
                    {
                        if(is_light(scene.shapes[path_vertex.shape_id]))
                        {
                            int light_id = get_area_light_id(scene.shapes[path_vertex.shape_id]);
                            Spectrum Le = emission(path_vertex, -ray.dir, scene);
                            PointAndNormal light_PN;
                            light_PN.position = path_vertex.position;
                            light_PN.normal = path_vertex.geometry_normal;
                            Real pdf_nee = pdf_point_on_light(scene.lights[light_id], light_PN, p_nee, scene) * light_pmf(scene, light_id);
                            Real pdf_dir_ = pdf_dir * multi_trans_prob * Geometry(light_PN, p_nee);
                            Real w = (pdf_dir * pdf_dir) / (pdf_dir * pdf_dir + pdf_nee * pdf_nee);
                            L += cur_path_throughput * Le * w;
                        }
                    }
                }
            }//else we sample a direction and calc the scattered radiance
            cur_path_throughput *= exp(-sigma_t * t) / pdf;
            ray.org += t * ray.dir;//now we go to this new place to calc the radiance scattered there
            Vector2 rnd_uv(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
            std::optional<Vector3> opt_dir = sample_phase_function(phaseFunction, ray.dir, rnd_uv);
            if (opt_dir)
            {
                Vector3 dir = opt_dir.value();
                Real phase_pdf = pdf_sample_phase(phaseFunction, dir, -ray.dir);
                pdf_dir = phase_pdf;
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

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implement this!
    return make_zero_spectrum();
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implement this!
    return make_zero_spectrum();
}
