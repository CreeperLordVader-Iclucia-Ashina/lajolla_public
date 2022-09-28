#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometry_normal, dir_in) < 0 ||
            dot(vertex.geometry_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    const Spectrum baseColor = eval(get_texture_op(), vertex.uv, vertex.uv_screen_size, texture_pool);
    const Real Roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    const Real Anisotropic = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    const Real aspect = sqrt(1.0 - 0.9 * Anisotropic);
    const Real alphax = max(0.0001, Roughness * Roughness / aspect);
    const Real alphay = max(0.0001, Roughness * Roughness * aspect);
    const Vector3 H = normalize(dir_in + dir_out);
    const Vector3 Hl = to_local(frame, H);
    const Spectrum F = baseColor + (Spectrum(1.0) - baseColor) * pow(1 - abs(dot(H, dir_out)), 5);
    const Real tmp = Hl[0] * Hl[0] / (alphax * alphax) + Hl[1] * Hl[1] / (alphay * alphay) + Hl[2] * Hl[2];
    const Real D = M_1_PI / (tmp * tmp);
    const Vector3 wi_l = to_local(frame, dir_in);
    const Vector3 wo_l = to_local(frame, dir_out);
    const Real G_in = 2.0 / (1.0 + sqrt(1.0 + (pow(wi_l.x * alphax, 2) + pow(wi_l.y * alphay, 2)) / (wi_l.z * wi_l.z)));
    const Real G_out = 2.0 / (1.0 + sqrt(1.0 + (pow(wo_l.x * alphax, 2) + pow(wo_l.y * alphay, 2)) / (wo_l.z * wo_l.z)));
    const Real G = G_in * G_out;
    return 0.25 * F * D * G / abs(dot(vertex.geometry_normal, dir_in));
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometry_normal, dir_in) < 0 ||
            dot(vertex.geometry_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    const Real Roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    const Real Anisotropic = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    const Real aspect = sqrt(1.0 - 0.9 * Anisotropic);
    const Real alphax = max(0.0001, Roughness * Roughness / aspect);
    const Real alphay = max(0.0001, Roughness * Roughness * aspect);
    const Vector3 H = normalize(dir_in + dir_out);
    const Vector3 Hl = to_local(frame, H);
    const Real tmp = Hl[0] * Hl[0] / (alphax * alphax) + Hl[1] * Hl[1] / (alphay * alphay) + Hl[2] * Hl[2];
    const Real D = M_1_PI / (tmp * tmp);
    const Vector3 wi_l = to_local(frame, dir_in);
    const Vector3 wo_l = to_local(frame, dir_out);
    const Real G_in = 2.0 / (1.0 + sqrt(1.0 + (pow(wi_l.x * alphax, 2) + pow(wi_l.y * alphay, 2)) / (wi_l.z * wi_l.z)));
    return 0.25 * D * G_in * abs(dot(wi_l, Hl)) / std::max(0.0001, abs(dot(wo_l, Hl)))
            / std::max(0.0001, abs(dot(wi_l, vertex.geometry_normal)));
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometry_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    const Real Roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    const Real Anisotropic = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    const Real aspect = sqrt(1.0 - 0.9 * Anisotropic);
    const Real alphax = max(0.0001, Roughness * Roughness / aspect);
    const Real alphay = max(0.0001, Roughness * Roughness * aspect);
    Vector3 wi = to_local(frame, dir_in); // map the incoming light to local frame
    wi = normalize(Vector3(wi.x * alphax, wi.y * alphay, wi.z));
    const Vector3 T1 = normalize(abs(wi.z) > 0.9999 ?
                        cross(Vector3(0, 1, 0), wi) : cross(Vector3(0, 0, 1), wi));
    const Vector3 T2 = cross(wi, T1);
    const Real R = std::sqrt(rnd_param_uv[0]);
    const Real phi = c_TWOPI * rnd_param_uv[1];
    Real t1 = R * std::cos(phi);
    Real t2 = R * std::sin(phi);
    t2 = 0.5 * (wi.z - 1.0) * std::sqrt(1.0 - t1 * t1) + 0.5 * (wi.z + 1.0) * t2;
    Vector3 d_out = t1 * T1 + t2 * T2 + std::sqrt(std::max(0.0, 1.0 - t1 * t1 - t2 * t2)) * wi;
    d_out = normalize(Vector3(wi.x / alphax, wi.y / alphay, wi.z));
    return BSDFSampleRecord{to_world(frame, d_out), Real(0), Roughness};
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
