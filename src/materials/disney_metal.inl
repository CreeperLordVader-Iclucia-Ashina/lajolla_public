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
    const Spectrum baseColor = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    const Real Roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    const Real Anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    const Real aspect = sqrt(1.0 - 0.9 * Anisotropic);
    const Real alphax = max(0.0001, Roughness * Roughness / aspect);
    const Real alphay = max(0.0001, Roughness * Roughness * aspect);
    const Vector3 H = normalize(dir_in + dir_out);
    const Vector3 Hl = to_local(frame, H);
    const Spectrum F = baseColor + (Spectrum(1.0, 1.0, 1.0) - baseColor) * pow(1 - abs(dot(H, dir_out)), 5);
    const Real tmp = Hl.x * Hl.x / (alphax * alphax) + Hl.y * Hl.y / (alphay * alphay) + Hl.z * Hl.z;
    const Real D = M_1_PI / (alphax * alphay * tmp * tmp);
    const Vector3 wi_l = to_local(frame, dir_in);
    const Vector3 wo_l = to_local(frame, dir_out);
    const Real G_in = 2.0 / (1.0 + sqrt(1.0 + (pow(wi_l.x * alphax, 2) + pow(wi_l.y * alphay, 2)) / (wi_l.z * wi_l.z)));
    const Real G_out = 2.0 / (1.0 + sqrt(1.0 + (pow(wo_l.x * alphax, 2) + pow(wo_l.y * alphay, 2)) / (wo_l.z * wo_l.z)));
    const Real G = G_in * G_out;
    const Spectrum ret = 0.25 * F * D * G / abs(wi_l.z);
    return ret;
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
    const Real Anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    const Real aspect = sqrt(1.0 - 0.9 * Anisotropic);
    const Real alphax = max(0.0001, Roughness * Roughness / aspect);
    const Real alphay = max(0.0001, Roughness * Roughness * aspect);
    const Vector3 Hl = to_local(frame, normalize(dir_in + dir_out));
    const Real tmp = Hl.x * Hl.x / (alphax * alphax) + Hl.y * Hl.y / (alphay * alphay) + Hl.z * Hl.z;
    const Real D = M_1_PI / (tmp * tmp * alphax * alphay);
    const Vector3 wi_l = to_local(frame, dir_in);
    const Real G_in = 2.0 / (1.0 + sqrt(1.0 + (pow(wi_l.x * alphax, 2) + pow(wi_l.y * alphay, 2)) / (wi_l.z * wi_l.z)));
    const Real ret = 0.25 * D * G_in / max(0.0001, abs(wi_l.z));
    return ret;
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
    const Real Anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    const Real aspect = sqrt(1.0 - 0.9 * Anisotropic);
    const Real alphax = max(0.0001, Roughness * Roughness / aspect);
    const Real alphay = max(0.0001, Roughness * Roughness * aspect);
    Vector3 wi_l = to_local(frame, dir_in); // map the incoming light to local frame
    const Vector3 wi = normalize(Vector3(wi_l.x * alphax, wi_l.y * alphay, wi_l.z));
    const Vector3 T1 = normalize(abs(wi.z) > 0.9999 ?
                        Vector3(1, 0, 0) : cross(Vector3(0, 0, 1), wi));
    const Vector3 T2 = cross(wi, T1);
    const Real R = sqrt(rnd_param_uv[0]);
    const Real phi = c_TWOPI * rnd_param_uv[1];
    Real t1 = R * cos(phi);
    Real t2 = R * sin(phi);
    t2 = 0.5 * (1.0 - wi.z) * sqrt(1.0 - t1 * t1) + 0.5 * (wi.z + 1.0) * t2;
    Vector3 Hl = t1 * T1 + t2 * T2 + sqrt(max(0.0, 1.0 - t1 * t1 - t2 * t2)) * wi;// in local frame but scaled
    Hl = normalize(Vector3(Hl.x * alphax, Hl.y * alphay, max(0.0, Hl.z))); // trans to local frame
    Vector3 wo_l = 2.0 * dot(Hl, wi_l) * Hl - wi_l;
    return BSDFSampleRecord{to_world(frame, wo_l), Real(0), Roughness};
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
