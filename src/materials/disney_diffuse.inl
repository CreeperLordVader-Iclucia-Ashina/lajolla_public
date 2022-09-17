#include "vector.h"
Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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
    const TVector3 H = normalize(dir_in + dir_out);
    Real tmp = dot(H, dir_out);
    const Real Roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    const Real FD90 = 0.5 + Real(2) * Roughness * tmp * tmp;
    const Real FD_wi = (1.0 + (FD90 - 1.0) * (1.0 - pow(abs(dot(vertex.geometry_normal, dir_in)), 5)));
    const Real FD_wo = (1.0 + (FD90 - 1.0) * (1.0 - pow(abs(dot(vertex.geometry_normal, dir_out)), 5)));
    const Spectrum BaseColor = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    const Spectrum f_baseDiffuse = BaseColor * M_1_PI * FD_wi * FD_wo * dot(vertex.geometry_normal, dir_out);
    const Real FSS90 = Roughness * tmp * tmp;
    const Real NdotWin = abs(dot(vertex.geometry_normal, dir_in));
    const Real NdotWout = abs(dot(vertex.geometry_normal, dir_out));
    const Real FSS_wi = (1.0 + (FSS90 - 1.0) * (1 - pow(NdotWin, 5)));
    const Real FSS_wo = (1.0 + (FSS90 - 1.0) * (1 - pow(NdotWout, 5)));
    const Spectrum f_subsurface = 1.25 * BaseColor * M_1_PI * (FSS_wi * FSS_wo * (1.0 / (NdotWin + NdotWout) - 0.5) + 0.5) * NdotWout;
    const Real Subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    return (1.0 - Subsurface) * f_baseDiffuse + Subsurface * f_subsurface;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    return BSDFSampleRecord{
            to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
            Real(0) /* eta */,  eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool)/* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
