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
    const Real FD90 = 0.5 + Real(2) * eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool)
                        * tmp * tmp;
    const Real FD_wi = (1.0 + (FD90 - 1.0) * (1.0 - pow(abs(dot(vertex.geometry_normal, dir_in)), 5)));
    const Real FD_wo = (1.0 + (FD90 - 1.0) * (1.0 - pow(abs(dot(vertex.geometry_normal, dir_out)), 5)));
    const Spectrum f_baseDiffuse = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool) * M_1_PI
                                    * FD_wi * FD_wo * dot(vertex.geometry_normal, dir_out);
    const Real FSS90 = bsdf.roughness * tmp * tmp;
    const NdotWin = ;
    const NdotWout = dot(vertex.n)
    const Real FSS_wi = (1.0 + (FSS90 - 1.0) * (1 - ))
    return make_zero_spectrum();
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
    return Real(0);
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
    return {};
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
