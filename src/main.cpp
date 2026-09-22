#include <GLFW/glfw3.h>
#include "gl_api.hpp"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "imgui_internal.h"
#include "../scomplex.hpp"
#include "expr.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace cx = std::complex::experimental::zhangzl;

enum class Mode { Function, Fractal };
enum class FnView { Mapping, Domain, Polya };

struct Pin {
    cx::complex<double> z;
};

struct App {
    Mode mode = Mode::Function;
    FnView fn_view = FnView::Domain;
    double system = -1.0;
    double center_x = 0.0, center_y = 0.0;
    double scale = 4.5;
    int max_iter = 256;
    int family = 0; // 0 M, 1 T, 2 B
    bool L = false;
    bool julia = false;
    double julia_re = -0.8, julia_im = 0.156;
    bool show_grid = true, show_axis = true;
    int map_step = 16;
    int polya_omega = 96;
    float polya_length = 2.5f;
    char expr[256] = "z^2";
    char calc_in[256] = "exp(i*pi)+1";
    char calc_z[64] = "1+i";
    std::string calc_out;
    std::vector<std::string> calc_hist;
    std::vector<Pin> pins;
    double probe_x = 0, probe_y = 0;
    bool has_probe = false;
    std::string status;
    bool expr_ok = true;
    ComplexExpr parser;

    int fbo_w = 0, fbo_h = 0;
    GLuint fbo = 0, color = 0, vao_quad = 0, vao_pts = 0, vbo_quad = 0, vbo_pts = 0;
    GLuint prog_fractal = 0, prog_domain = 0, prog_points = 0;
    bool fractal_double = true;
    bool view_dirty = true;
    bool control_open = false;
};

static const char* kVSQuad = R"(
#version 400
layout(location=0) in vec2 a_pos;
void main(){ gl_Position = vec4(a_pos,0.0,1.0); }
)";

static const char* kFSFractal = R"(
#version 400
out vec4 frag;
uniform dvec2 u_center;
uniform double u_scale;
uniform vec2 u_resolution;
uniform int u_family;
uniform int u_L;
uniform int u_max_iter;
uniform int u_julia;
uniform dvec2 u_seed;

dvec2 cmul(dvec2 a, dvec2 b){
    return dvec2(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
}

void main(){
    double m = double(min(u_resolution.x, u_resolution.y));
    dvec2 pix = dvec2(
        u_center.x + u_scale * (double(gl_FragCoord.x) - double(u_resolution.x)*0.5) / m,
        u_center.y + u_scale * (double(gl_FragCoord.y) - double(u_resolution.y)*0.5) / m);
    dvec2 c = (u_julia != 0) ? u_seed : pix;
    dvec2 z = (u_julia != 0) ? pix : dvec2(0.0);
    int i;
    double mag = 0.0;
    for(i=0;i<u_max_iter;++i){
        if(u_family==1) z = dvec2(z.x, -z.y);
        if(u_family==2) z = dvec2(abs(z.x), abs(z.y));
        dvec2 f = cmul(z,z);
        z = (u_L != 0) ? (z - f + c) : (f + c);
        mag = z.x*z.x + z.y*z.y;
        if(mag > 256.0) break;
    }
    if(i>=u_max_iter){ frag = vec4(0.02,0.02,0.05,1); return; }
    float nu = float(i) + 1.0 - log(log(float(mag))*0.5)/log(2.0);
    float t = nu / float(max(u_max_iter,1));
    vec3 col = vec3(0.5+0.5*cos(0.15*nu+0.0), 0.5+0.5*cos(0.15*nu+2.1), 0.5+0.5*cos(0.15*nu+4.2));
    frag = vec4(col*(0.35+0.65*t), 1.0);
}
)";

static const char* kFSFractalF = R"(
#version 400
out vec4 frag;
uniform vec2 u_center;
uniform float u_scale;
uniform vec2 u_resolution;
uniform int u_family;
uniform int u_L;
uniform int u_max_iter;
uniform int u_julia;
uniform vec2 u_seed;
vec2 cmul(vec2 a, vec2 b){ return vec2(a.x*b.x-a.y*b.y, a.x*b.y+a.y*b.x); }
void main(){
    float m = min(u_resolution.x, u_resolution.y);
    vec2 pix = vec2(
        u_center.x + u_scale * (gl_FragCoord.x - u_resolution.x*0.5) / m,
        u_center.y + u_scale * (gl_FragCoord.y - u_resolution.y*0.5) / m);
    vec2 c = (u_julia != 0) ? u_seed : pix;
    vec2 z = (u_julia != 0) ? pix : vec2(0.0);
    int i;
    float mag = 0.0;
    for(i=0;i<u_max_iter;++i){
        if(u_family==1) z = vec2(z.x, -z.y);
        if(u_family==2) z = vec2(abs(z.x), abs(z.y));
        vec2 f = cmul(z,z);
        z = (u_L != 0) ? (z - f + c) : (f + c);
        mag = dot(z,z);
        if(mag > 256.0) break;
    }
    if(i>=u_max_iter){ frag = vec4(0.02,0.02,0.05,1); return; }
    float nu = float(i) + 1.0 - log(log(mag)*0.5)/log(2.0);
    float t = nu / float(max(u_max_iter,1));
    vec3 col = vec3(0.5+0.5*cos(0.15*nu+0.0), 0.5+0.5*cos(0.15*nu+2.1), 0.5+0.5*cos(0.15*nu+4.2));
    frag = vec4(col*(0.35+0.65*t), 1.0);
}
)";

static const char* kFSDomain = R"(
#version 400
out vec4 frag;
uniform vec2 u_center;
uniform float u_scale;
uniform vec2 u_resolution;
uniform int u_fn;
vec2 cmul(vec2 a, vec2 b){ return vec2(a.x*b.x-a.y*b.y, a.x*b.y+a.y*b.x); }
vec2 cdiv(vec2 a, vec2 b){
    float d = dot(b,b)+1e-12;
    return vec2(a.x*b.x+a.y*b.y, a.y*b.x-a.x*b.y)/d;
}
vec2 cexp(vec2 z){ float e=exp(z.x); return e*vec2(cos(z.y),sin(z.y)); }
vec2 clog(vec2 z){ return vec2(log(length(z)+1e-12), atan(z.y,z.x)); }
vec2 csin(vec2 z){
    vec2 iz=vec2(-z.y,z.x);
    vec2 e=cexp(iz), ei=cexp(-iz);
    vec2 d=e-ei;
    return vec2(d.y, -d.x)*0.5;
}
vec2 ccos(vec2 z){
    vec2 iz=vec2(-z.y,z.x);
    return 0.5*(cexp(iz)+cexp(-iz));
}
vec3 hsv(float h, float s, float v){
    vec3 k = vec3(1.0, 2.0/3.0, 1.0/3.0);
    vec3 p = abs(fract(h+k)*6.0-3.0);
    return v * mix(vec3(1.0), clamp(p-1.0,0.0,1.0), s);
}
vec2 apply_f(vec2 z){
    if(u_fn==1) return cmul(z,z);
    if(u_fn==2) return cmul(z, cmul(z,z));
    if(u_fn==3) return cdiv(vec2(1,0), z);
    if(u_fn==4) return cexp(z);
    if(u_fn==5) return csin(z);
    if(u_fn==6) return ccos(z);
    if(u_fn==7) return clog(z);
    if(u_fn==8) {
        float r=sqrt(length(z)); float t=0.5*atan(z.y,z.x);
        return r*vec2(cos(t),sin(t));
    }
    return z;
}
void main(){
    float m = min(u_resolution.x, u_resolution.y);
    vec2 z = vec2(
        u_center.x + u_scale * (gl_FragCoord.x - u_resolution.x*0.5) / m,
        u_center.y + u_scale * (gl_FragCoord.y - u_resolution.y*0.5) / m);
    vec2 w = apply_f(z);
    float r = length(w);
    float arg = atan(w.y, w.x);
    float h = arg / 6.28318530718;
    float v = 0.25 + 0.75 * fract(log(r+1e-8)/log(2.0));
    frag = vec4(hsv(h, 0.85, v), 1.0);
}
)";

static const char* kVSPoints = R"(
#version 400
layout(location=0) in vec2 a_z;
layout(location=1) in vec4 a_col;
layout(location=2) in float a_size;
uniform vec2 u_center;
uniform float u_scale;
uniform vec2 u_resolution;
out vec4 v_col;
void main(){
    float m = min(u_resolution.x, u_resolution.y);
    vec2 p = vec2(
        u_resolution.x*0.5 + (a_z.x - u_center.x) / u_scale * m,
        u_resolution.y*0.5 + (a_z.y - u_center.y) / u_scale * m);
    vec2 ndc = p / u_resolution * 2.0 - 1.0;
    gl_Position = vec4(ndc, 0.0, 1.0);
    gl_PointSize = a_size;
    v_col = a_col;
}
)";

static const char* kFSPoints = R"(
#version 400
in vec4 v_col;
out vec4 frag;
void main(){ frag = v_col; }
)";

static App g;

static CompiledExpr g_fx;
static std::string g_fx_src;

static bool refresh_compiled() {
    if (g_fx.root >= 0 && g_fx_src == g.expr) return true;
    try {
        g_fx = compile_expr(g.expr);
        g_fx_src = g.expr;
        return true;
    } catch (const std::exception& e) {
        g.status = e.what();
        g_fx = CompiledExpr{};
        g_fx_src.clear();
        return false;
    }
}

static cx::complex<double> eval_f_raw(cx::complex<double> z, cx::complex<double> c = cx::complex<double>()) {
    return eval_compiled(g_fx, z, c);
}

static cx::complex<double> eval_f(cx::complex<double> z) {
    cx::complex<double>::set_system(g.system);
    if (!refresh_compiled()) throw ExprError(g.status.empty() ? "bad expression" : g.status);
    return eval_f_raw(z);
}

static float polya_fade(float t) {
    const float e = 0.30f;
    t = std::clamp(t, 0.f, 1.f);
    if (t < e) {
        float u = t / e;
        return u * u * (3.f - 2.f * u);
    }
    if (t > 1.f - e) {
        float u = (1.f - t) / e;
        return u * u * (3.f - 2.f * u);
    }
    return 1.f;
}

static bool polya_field(cx::complex<double> z, cx::complex<double>& F) {
    try {
        auto fz = eval_f_raw(z);
        double n2 = cx::norm_norm(fz);
        if (!(n2 > 1e-30) || n2 > 1e300) return false;
        F = cx::conj(fz) * (1.0 / std::sqrt(n2));
        return true;
    } catch (...) {
        return false;
    }
}

static bool polya_rk4(cx::complex<double>& z, double h) {
    cx::complex<double> k1, k2, k3, k4;
    if (!polya_field(z, k1)) return false;
    if (!polya_field(z + k1 * (0.5 * h), k2)) return false;
    if (!polya_field(z + k2 * (0.5 * h), k3)) return false;
    if (!polya_field(z + k3 * h, k4)) return false;
    auto dz = (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0);
    double jump2 = cx::norm_norm(dz);
    if (!(jump2 > 1e-24)) return false;
    if (jump2 > g.scale * g.scale * 0.25) return false;
    z = z + dz;
    return true;
}

static int gpu_fn_id() {
    std::string e = g.expr;
    if (e == "z" || e == "x") return 0;
    if (e == "z^2" || e == "x^2") return 1;
    if (e == "z^3" || e == "x^3") return 2;
    if (e == "1/z" || e == "z^(-1)") return 3;
    if (e == "exp(z)") return 4;
    if (e == "sin(z)") return 5;
    if (e == "cos(z)") return 6;
    if (e == "log(z)") return 7;
    if (e == "sqrt(z)") return 8;
    return -1;
}

static bool fractal_use_gpu() {
    if (g.system != -1.0) return false;
    std::string e = g.expr;
    if (e == "z^2" || e == "x^2") return true;
    if (!g.L && (e == "z^2+c" || e == "x^2+c")) return true;
    return false;
}

static cx::complex<double> fractal_apply(cx::complex<double> a, cx::complex<double> c) {
    if (g.family == 1) a = cx::conj(a);
    else if (g.family == 2) a = cx::abs(a);
    auto fa = eval_f_raw(a, c);
    if (g_fx.uses_c) return fa;
    return g.L ? (a - fa + c) : (fa + c);
}

static cx::complex<double> from_pixel(double px, double py, double w, double h) {
    double m = std::min(w, h);
    return cx::complex<double>(
        g.center_x + g.scale * (px - w * 0.5) / m,
        g.center_y + g.scale * (py - h * 0.5) / m);
}

static void hsv_rgb(double h, double s, double v, float& r, float& gch, float& b) {
    double c = v * s;
    double x = c * (1.0 - std::fabs(std::fmod(h * 6.0, 2.0) - 1.0));
    double m = v - c;
    double rp = 0, gp = 0, bp = 0;
    int i = (int)std::floor(h * 6.0) % 6;
    if (i < 0) i += 6;
    if (i == 0) { rp = c; gp = x; }
    else if (i == 1) { rp = x; gp = c; }
    else if (i == 2) { gp = c; bp = x; }
    else if (i == 3) { gp = x; bp = c; }
    else if (i == 4) { rp = x; bp = c; }
    else { rp = c; bp = x; }
    r = (float)(rp + m); gch = (float)(gp + m); b = (float)(bp + m);
}

static void ensure_fbo(int w, int h) {
    if (w < 1) w = 1;
    if (h < 1) h = 1;
    if (g.fbo && g.fbo_w == w && g.fbo_h == h) return;
    if (g.color) glDeleteTextures(1, &g.color);
    if (g.fbo) glDeleteFramebuffers(1, &g.fbo);
    glGenFramebuffers(1, &g.fbo);
    glGenTextures(1, &g.color);
    glBindTexture(GL_TEXTURE_2D, g.color);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
    g.view_dirty = true;
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindFramebuffer(GL_FRAMEBUFFER, g.fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, g.color, 0);
    g.fbo_w = w;
    g.fbo_h = h;
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

static void draw_quad(GLuint prog) {
    glUseProgram(prog);
    glBindVertexArray(g.vao_quad);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
}

struct Pt {
    float x, y;
    float r, gc, b, a;
    float sz;
};

static void upload_pts(const std::vector<Pt>& pts, GLenum mode, bool blend) {
    glUseProgram(g.prog_points);
    glUniform2f(glGetUniformLocation(g.prog_points, "u_center"), (float)g.center_x, (float)g.center_y);
    glUniform1f(glGetUniformLocation(g.prog_points, "u_scale"), (float)g.scale);
    glUniform2f(glGetUniformLocation(g.prog_points, "u_resolution"), (float)g.fbo_w, (float)g.fbo_h);
    glBindVertexArray(g.vao_pts);
    glBindBuffer(GL_ARRAY_BUFFER, g.vbo_pts);
    glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)(pts.size() * sizeof(Pt)), pts.data(), GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(Pt), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Pt), (void*)(2 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(Pt), (void*)(6 * sizeof(float)));
    glEnable(GL_PROGRAM_POINT_SIZE);
    if (blend) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }
    glDrawArrays(mode, 0, (GLint)pts.size());
    if (blend) glDisable(GL_BLEND);
}

static void collect_grid(int w, int h, int step, std::vector<int>& xs, std::vector<int>& ys) {
    xs.clear();
    ys.clear();
    for (int y = step / 2; y < h; y += step)
        for (int x = step / 2; x < w; x += step) {
            xs.push_back(x);
            ys.push_back(y);
        }
}

static void render_mapping() {
    if (!refresh_compiled()) return;
    int w = g.fbo_w, h = g.fbo_h;
    int step = std::max(4, g.map_step);
    std::vector<int> xs, ys;
    collect_grid(w, h, step, xs, ys);
    int n = (int)xs.size();
    std::vector<Pt> pts;
    pts.resize((size_t)n * 2);
    std::vector<char> ok((size_t)n, 0);
    cx::complex<double>::set_system(g.system);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8)
#endif
    for (int i = 0; i < n; ++i) {
        try {
            auto z = from_pixel(xs[i], ys[i], w, h);
            auto fz = eval_f_raw(z);
            float r = (float)(xs[i] * 0.5 / w + ys[i] * 0.5 / h);
            float gc = (float)xs[i] / w;
            float b = (float)ys[i] / h;
            pts[(size_t)i * 2] = { (float)cx::real(z), (float)cx::imag(z), r, gc, b, 1.f, 3.f };
            pts[(size_t)i * 2 + 1] = { (float)cx::real(fz), (float)cx::imag(fz), r, gc, b, 1.f, 6.f };
            ok[i] = 1;
        } catch (...) {
            ok[i] = 0;
        }
    }
    std::vector<Pt> packed;
    packed.reserve(pts.size());
    for (int i = 0; i < n; ++i) if (ok[i]) {
        packed.push_back(pts[(size_t)i * 2]);
        packed.push_back(pts[(size_t)i * 2 + 1]);
    }
    if (!packed.empty()) upload_pts(packed, GL_POINTS, false);
}

static void render_polya() {
    if (!refresh_compiled()) return;
    int w = g.fbo_w, h = g.fbo_h;
    int step = std::max(8, g.map_step);
    const int OmegaI = std::max(8, g.polya_omega);
    const double dt = (double)g.polya_length / (double)OmegaI;
    std::vector<int> xs, ys;
    collect_grid(w, h, step, xs, ys);
    int n = (int)xs.size();
    std::vector<std::vector<Pt>> buckets((size_t)n);
    cx::complex<double>::set_system(g.system);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 2)
#endif
    for (int i = 0; i < n; ++i) {
        try {
            auto z = from_pixel(xs[i], ys[i], w, h);
            float r = (float)(xs[i] * 0.25 / w + ys[i] * 0.25 / h);
            float gc = (float)xs[i] / w * 0.55f;
            float b = (float)ys[i] / h * 0.55f;
            std::vector<cx::complex<double>> back, fwd;
            back.reserve((size_t)OmegaI + 1);
            fwd.reserve((size_t)OmegaI + 1);
            back.push_back(z);
            auto zz = z;
            for (int k = 0; k < OmegaI; ++k) {
                if (!polya_rk4(zz, -dt)) break;
                back.push_back(zz);
            }
            fwd.push_back(z);
            zz = z;
            for (int k = 0; k < OmegaI; ++k) {
                if (!polya_rk4(zz, dt)) break;
                fwd.push_back(zz);
            }
            std::vector<cx::complex<double>> chain;
            chain.reserve(back.size() + fwd.size());
            for (int k = (int)back.size() - 1; k >= 0; --k) chain.push_back(back[(size_t)k]);
            for (size_t k = 1; k < fwd.size(); ++k) chain.push_back(fwd[k]);
            int m = (int)chain.size();
            if (m < 2) continue;
            auto& out = buckets[(size_t)i];
            out.reserve((size_t)(m - 1) * 2);
            float denom = (float)(m - 1);
            for (int k = 0; k < m - 1; ++k) {
                float w0 = polya_fade((float)k / denom);
                float w1 = polya_fade((float)(k + 1) / denom);
                out.push_back({ (float)cx::real(chain[(size_t)k]), (float)cx::imag(chain[(size_t)k]),
                    r * w0, gc * w0, b * w0, w0, 1.f });
                out.push_back({ (float)cx::real(chain[(size_t)k + 1]), (float)cx::imag(chain[(size_t)k + 1]),
                    r * w1, gc * w1, b * w1, w1, 1.f });
            }
        } catch (...) {}
    }
    std::vector<Pt> pts;
    size_t total = 0;
    for (auto& bkt : buckets) total += bkt.size();
    pts.reserve(total);
    for (auto& bkt : buckets) {
        pts.insert(pts.end(), bkt.begin(), bkt.end());
    }
    if (!pts.empty()) upload_pts(pts, GL_LINES, true);
}

static void render_domain_cpu() {
    if (!refresh_compiled()) return;
    int w = g.fbo_w, h = g.fbo_h;
    int step = (w * h > 400000) ? 2 : 1;
    std::vector<unsigned char> img((size_t)w * h * 3, 8);
    cx::complex<double>::set_system(g.system);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int y = 0; y < h; y += step) {
        for (int x = 0; x < w; x += step) {
            try {
                auto z = from_pixel(x + 0.5, y + 0.5, w, h);
                auto fz = eval_f_raw(z);
                auto lg = cx::log(fz);
                double arg = cx::imag(lg);
                double rr = std::sqrt(std::max(0.0, cx::norm_norm(fz)));
                float R, G, B;
                double hue = arg / (2.0 * 3.14159265358979323846);
                hue -= std::floor(hue);
                double val = 0.25 + 0.75 * std::fmod(std::log(rr + 1e-12) / std::log(2.0) + 1000.0, 1.0);
                hsv_rgb(hue, 0.85, val, R, G, B);
                for (int dy = 0; dy < step && y + dy < h; ++dy)
                    for (int dx = 0; dx < step && x + dx < w; ++dx) {
                        size_t o = ((size_t)(y + dy) * w + (x + dx)) * 3;
                        img[o] = (unsigned char)std::clamp(R * 255.0f, 0.f, 255.f);
                        img[o + 1] = (unsigned char)std::clamp(G * 255.0f, 0.f, 255.f);
                        img[o + 2] = (unsigned char)std::clamp(B * 255.0f, 0.f, 255.f);
                    }
            } catch (...) {}
        }
    }
    glBindTexture(GL_TEXTURE_2D, g.color);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, img.data());
}

static void render_fractal_cpu() {
    if (!refresh_compiled()) return;
    int w = g.fbo_w, h = g.fbo_h;
    int maxit = std::max(1, g.max_iter);
    std::vector<unsigned char> img((size_t)w * h * 3, 5);
    cx::complex<double>::set_system(g.system);
    const cx::complex<double> seed(g.julia_re, g.julia_im);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4)
#endif
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            auto pix = from_pixel(x + 0.5, y + 0.5, w, h);
            cx::complex<double> c = g.julia ? seed : pix;
            cx::complex<double> z = g.julia ? pix : cx::complex<double>(0, 0);
            int it = 0;
            double mag = 0.0;
            try {
                for (; it < maxit; ++it) {
                    z = fractal_apply(z, c);
                    mag = cx::norm_norm(z);
                    if (mag > 256.0) break;
                }
            } catch (...) {
                it = maxit;
            }
            size_t o = ((size_t)y * (size_t)w + (size_t)x) * 3;
            if (it >= maxit) {
                img[o] = 5;
                img[o + 1] = 5;
                img[o + 2] = 13;
                continue;
            }
            float nu = (float)it + 1.0f - std::log(std::log((float)mag) * 0.5f) / std::log(2.0f);
            float t = nu / (float)maxit;
            float R = 0.5f + 0.5f * std::cos(0.15f * nu);
            float Gch = 0.5f + 0.5f * std::cos(0.15f * nu + 2.1f);
            float B = 0.5f + 0.5f * std::cos(0.15f * nu + 4.2f);
            float s = 0.35f + 0.65f * t;
            img[o] = (unsigned char)std::clamp(R * s * 255.0f, 0.f, 255.f);
            img[o + 1] = (unsigned char)std::clamp(Gch * s * 255.0f, 0.f, 255.f);
            img[o + 2] = (unsigned char)std::clamp(B * s * 255.0f, 0.f, 255.f);
        }
    }
    glBindTexture(GL_TEXTURE_2D, g.color);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, img.data());
}

static void set_fractal_uniforms() {
    glUseProgram(g.prog_fractal);
    if (g.fractal_double && glUniform2d && glUniform1d) {
        glUniform2d(glGetUniformLocation(g.prog_fractal, "u_center"), g.center_x, g.center_y);
        glUniform1d(glGetUniformLocation(g.prog_fractal, "u_scale"), g.scale);
        glUniform2d(glGetUniformLocation(g.prog_fractal, "u_seed"), g.julia_re, g.julia_im);
    } else {
        glUniform2f(glGetUniformLocation(g.prog_fractal, "u_center"), (float)g.center_x, (float)g.center_y);
        glUniform1f(glGetUniformLocation(g.prog_fractal, "u_scale"), (float)g.scale);
        glUniform2f(glGetUniformLocation(g.prog_fractal, "u_seed"), (float)g.julia_re, (float)g.julia_im);
    }
    glUniform2f(glGetUniformLocation(g.prog_fractal, "u_resolution"), (float)g.fbo_w, (float)g.fbo_h);
    glUniform1i(glGetUniformLocation(g.prog_fractal, "u_family"), g.family);
    glUniform1i(glGetUniformLocation(g.prog_fractal, "u_L"), g.L ? 1 : 0);
    glUniform1i(glGetUniformLocation(g.prog_fractal, "u_max_iter"), g.max_iter);
    glUniform1i(glGetUniformLocation(g.prog_fractal, "u_julia"), g.julia ? 1 : 0);
}

static void render_view() {
    bool gpu = (g.mode == Mode::Fractal && fractal_use_gpu()) ||
               (g.mode == Mode::Function && g.fn_view == FnView::Domain && g.system == -1.0 && gpu_fn_id() >= 0);
    if (!gpu && !g.view_dirty) return;
    glBindFramebuffer(GL_FRAMEBUFFER, g.fbo);
    glViewport(0, 0, g.fbo_w, g.fbo_h);
    glDisable(GL_DEPTH_TEST);
    glClearColor(0.07f, 0.07f, 0.10f, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    if (g.mode == Mode::Fractal) {
        if (fractal_use_gpu()) {
            set_fractal_uniforms();
            draw_quad(g.prog_fractal);
        } else {
            render_fractal_cpu();
        }
    } else if (g.fn_view == FnView::Domain) {
        int id = (g.system == -1.0) ? gpu_fn_id() : -1;
        if (id >= 0) {
            glUseProgram(g.prog_domain);
            glUniform2f(glGetUniformLocation(g.prog_domain, "u_center"), (float)g.center_x, (float)g.center_y);
            glUniform1f(glGetUniformLocation(g.prog_domain, "u_scale"), (float)g.scale);
            glUniform2f(glGetUniformLocation(g.prog_domain, "u_resolution"), (float)g.fbo_w, (float)g.fbo_h);
            glUniform1i(glGetUniformLocation(g.prog_domain, "u_fn"), id);
            draw_quad(g.prog_domain);
        } else {
            render_domain_cpu();
        }
    } else if (g.fn_view == FnView::Mapping) {
        render_mapping();
    } else {
        render_polya();
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    g.view_dirty = false;
}

static void draw_overlay(ImVec2 origin, ImVec2 size) {
    ImDrawList* dl = ImGui::GetWindowDrawList();
    double m = std::min(size.x, size.y);
    auto to_screen = [&](double re, double im) {
        return ImVec2(
            (float)(origin.x + size.x * 0.5 + (re - g.center_x) / g.scale * m),
            (float)(origin.y + size.y * 0.5 - (im - g.center_y) / g.scale * m));
    };
    if (g.show_grid) {
        double step = std::pow(10.0, std::floor(std::log10(g.scale / 4.0)));
        if (step <= 0) step = g.scale / 4.0;
        dl->PushClipRect(origin, ImVec2(origin.x + size.x, origin.y + size.y), true);
        for (int k = -40; k <= 40; ++k) {
            double xv = std::floor(g.center_x / step) * step + k * step;
            double yv = std::floor(g.center_y / step) * step + k * step;
            auto a = to_screen(xv, g.center_y - g.scale);
            auto b = to_screen(xv, g.center_y + g.scale);
            dl->AddLine(ImVec2(a.x, origin.y), ImVec2(b.x, origin.y + size.y), IM_COL32(80, 80, 90, 90));
            a = to_screen(g.center_x - g.scale, yv);
            b = to_screen(g.center_x + g.scale, yv);
            dl->AddLine(ImVec2(origin.x, a.y), ImVec2(origin.x + size.x, b.y), IM_COL32(80, 80, 90, 90));
        }
        dl->PopClipRect();
    }
    if (g.show_axis) {
        auto o = to_screen(0, 0);
        dl->AddLine(ImVec2(o.x, origin.y), ImVec2(o.x, origin.y + size.y), IM_COL32(220, 220, 220, 160));
        dl->AddLine(ImVec2(origin.x, o.y), ImVec2(origin.x + size.x, o.y), IM_COL32(220, 220, 220, 160));
    }
    for (auto& p : g.pins) {
        auto a = to_screen(cx::real(p.z), cx::imag(p.z));
        dl->AddCircleFilled(a, 4.f, IM_COL32(255, 220, 80, 255));
        try {
            auto fz = eval_f(p.z);
            auto b = to_screen(cx::real(fz), cx::imag(fz));
            dl->AddLine(a, b, IM_COL32(255, 180, 60, 180), 1.5f);
            dl->AddCircleFilled(b, 5.f, IM_COL32(80, 200, 255, 255));
        } catch (...) {}
    }
    if (g.has_probe) {
        auto a = to_screen(g.probe_x, g.probe_y);
        dl->AddCircle(a, 6.f, IM_COL32(255, 255, 255, 220), 12, 1.5f);
    }
}

static void run_calc() {
    try {
        cx::complex<double>::set_system(g.system);
        auto z = g.parser.eval(g.calc_z, cx::complex<double>(0, 0));
        auto w = g.parser.eval(g.calc_in, z);
        char buf[256];
        std::snprintf(buf, sizeof(buf), "%.12g %+.12g i", cx::real(w), cx::imag(w));
        g.calc_out = buf;
        g.calc_hist.insert(g.calc_hist.begin(), std::string(g.calc_in) + " | z=" + g.calc_z + " -> " + g.calc_out);
        if (g.calc_hist.size() > 24) g.calc_hist.pop_back();
        g.status.clear();
    } catch (const std::exception& e) {
        g.calc_out = e.what();
        g.status = e.what();
    }
}

static void ui_preset_buttons() {
    const char* presets[] = {
        "z", "z^2", "z^3", "z^2+c", "1/z", "z^z", "exp(z)", "sin(z)", "cos(z)", "log(z)",
        "sqrt(z)", "(1+z)/(1-z)", "z^2-z+i", "gamma(z)"
    };
    ImGuiStyle& st = ImGui::GetStyle();
    float wrap = ImGui::GetWindowPos().x + ImGui::GetContentRegionAvail().x;
    for (int n = 0; n < (int)(sizeof(presets) / sizeof(presets[0])); ++n) {
        float bw = ImGui::CalcTextSize(presets[n]).x + st.FramePadding.x * 2.0f;
        if (n > 0 && ImGui::GetItemRectMax().x + st.ItemSpacing.x + bw >= wrap)
            ImGui::NewLine();
        else if (n > 0)
            ImGui::SameLine();
        if (ImGui::Button(presets[n])) {
            std::snprintf(g.expr, sizeof(g.expr), "%s", presets[n]);
            g.view_dirty = true;
        }
    }
}

static void ui_controls_body() {
    if (ImGui::RadioButton("Function", g.mode == Mode::Function)) { g.mode = Mode::Function; g.view_dirty = true; }
    ImGui::SameLine();
    if (ImGui::RadioButton("Fractal", g.mode == Mode::Fractal)) { g.mode = Mode::Fractal; g.view_dirty = true; }

    const char* sys[] = { "Complex  i^2=-1", "Split  j^2=+1", "Dual  e^2=0" };
    int si = g.system < 0 ? 0 : (g.system > 0 ? 1 : 2);
    if (ImGui::Combo("Algebra", &si, sys, 3)) {
        g.system = si == 0 ? -1.0 : (si == 1 ? 1.0 : 0.0);
        g.view_dirty = true;
    }

    ImGui::Checkbox("Grid", &g.show_grid);
    ImGui::SameLine();
    ImGui::Checkbox("Axes", &g.show_axis);
    if (ImGui::Button("Reset view")) {
        g.center_x = g.center_y = 0;
        g.scale = g.mode == Mode::Fractal ? 3.5 : 4.5;
        g.view_dirty = true;
    }
    ImGui::Text("scale = %.3e", g.scale);
    if (g.scale < 1e-13)
        ImGui::TextWrapped("Double precision is exhausted; zoom out or the picture will pixelate.");

    ImGui::SeparatorText("Formula f(z)");
    if (ImGui::InputText("f(z)", g.expr, sizeof(g.expr))) g.view_dirty = true;
    ImGui::TextUnformatted("Presets:");
    ui_preset_buttons();
    if (g.mode == Mode::Fractal) {
        ImGui::TextWrapped("Iterate via scomplex M/T/B on this f. If the formula contains c, it is the full iterator (L is ignored). z^2 uses the GPU; other formulas use OpenMP.");
    }

    if (g.mode == Mode::Function) {
        ImGui::SeparatorText("Visualization");
        if (ImGui::RadioButton("Domain coloring", g.fn_view == FnView::Domain)) { g.fn_view = FnView::Domain; g.view_dirty = true; }
        if (ImGui::RadioButton("Mapping (ComplexGeom)", g.fn_view == FnView::Mapping)) { g.fn_view = FnView::Mapping; g.view_dirty = true; }
        if (ImGui::RadioButton("Polya field", g.fn_view == FnView::Polya)) { g.fn_view = FnView::Polya; g.view_dirty = true; }
        if (g.fn_view != FnView::Domain)
            if (ImGui::SliderInt("sample step", &g.map_step, 6, 40)) g.view_dirty = true;
        if (g.fn_view == FnView::Polya) {
            if (ImGui::SliderInt("Polya steps", &g.polya_omega, 16, 256)) g.view_dirty = true;
            if (ImGui::SliderFloat("Polya length", &g.polya_length, 0.5f, 8.0f, "%.2f")) g.view_dirty = true;
            ImGui::TextUnformatted("Integrator: RK4 on conj(f)/|f|.");
        }
    } else {
        ImGui::SeparatorText("M / T / B  (scomplex)");
        if (ImGui::RadioButton("M  f then +c", g.family == 0)) { g.family = 0; g.view_dirty = true; }
        if (ImGui::RadioButton("T  conj then M", g.family == 1)) { g.family = 1; g.view_dirty = true; }
        if (ImGui::RadioButton("B  abs then M", g.family == 2)) { g.family = 2; g.view_dirty = true; }
        refresh_compiled();
        if (g_fx.uses_c)
            ImGui::TextDisabled("L ignored: formula already contains c");
        else if (ImGui::Checkbox("L  (z - f(z) + c)", &g.L)) g.view_dirty = true;
        if (ImGui::Checkbox("Julia", &g.julia)) g.view_dirty = true;
        if (ImGui::InputDouble("Julia Re", &g.julia_re, 0.0, 0.0, "%.12f")) g.view_dirty = true;
        if (ImGui::InputDouble("Julia Im", &g.julia_im, 0.0, 0.0, "%.12f")) g.view_dirty = true;
        if (ImGui::SliderInt("iterations", &g.max_iter, 32, 2000)) g.view_dirty = true;
        ImGui::TextUnformatted("Wheel: zoom at cursor. Drag: pan. Right-click: Julia seed.");
        if (fractal_use_gpu())
            ImGui::TextUnformatted("Renderer: GPU (z^2).");
        else
            ImGui::TextUnformatted("Renderer: CPU OpenMP + scomplex.");
        if (!g.fractal_double)
            ImGui::TextWrapped("GPU has no double; zoom uses float (about 1e-6).");
    }

    ImGui::SeparatorText("Calculator");
    if (ImGui::InputText("expr", g.calc_in, sizeof(g.calc_in), ImGuiInputTextFlags_EnterReturnsTrue))
        run_calc();
    ImGui::InputText("z =", g.calc_z, sizeof(g.calc_z));
    if (ImGui::Button("Evaluate"))
        run_calc();
    ImGui::SameLine();
    if (ImGui::Button("Plot as f")) { std::snprintf(g.expr, sizeof(g.expr), "%s", g.calc_in); g.view_dirty = true; }
    ImGui::SameLine();
    if (ImGui::Button("Use as Julia c")) {
        try {
            cx::complex<double>::set_system(g.system);
            auto c = g.parser.eval(g.calc_in, g.parser.eval(g.calc_z, cx::complex<double>(0, 0)));
            g.julia_re = cx::real(c);
            g.julia_im = cx::imag(c);
            g.julia = true;
            g.mode = Mode::Fractal;
            g.view_dirty = true;
        } catch (const std::exception& e) { g.status = e.what(); }
    }
    ImGui::TextWrapped("%s", g.calc_out.c_str());
    if (ImGui::BeginListBox("history", ImVec2(-1, 160))) {
        for (auto& h : g.calc_hist) ImGui::Selectable(h.c_str(), false);
        ImGui::EndListBox();
    }

    if (g.has_probe) {
        ImGui::SeparatorText("Probe");
        ImGui::Text("z = %.10g %+.10g i", g.probe_x, g.probe_y);
        try {
            auto fz = eval_f(cx::complex<double>(g.probe_x, g.probe_y));
            ImGui::Text("f(z) = %.10g %+.10g i", cx::real(fz), cx::imag(fz));
        } catch (const std::exception& e) {
            ImGui::Text("f(z): %s", e.what());
        }
        if (ImGui::Button("Pin z")) g.pins.push_back({ cx::complex<double>(g.probe_x, g.probe_y) });
        ImGui::SameLine();
        if (ImGui::Button("Clear pins")) g.pins.clear();
    }
    if (!g.status.empty()) ImGui::TextColored(ImVec4(1, 0.4f, 0.3f, 1), "%s", g.status.c_str());
#ifdef _OPENMP
    ImGui::Text("OpenMP threads: %d", omp_get_max_threads());
#endif
}

static void close_controls() {
    g.control_open = false;
    g.view_dirty = true;
}

static void ui_control_screen() {
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(vp->WorkPos);
    ImGui::SetNextWindowSize(vp->WorkSize);
    ImGui::Begin("Control", nullptr,
        ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoSavedSettings);
    if (ImGui::Button("Back to Plane") || ImGui::IsKeyPressed(ImGuiKey_Escape))
        close_controls();
    ImGui::SameLine();
    ImGui::TextUnformatted("Control");
    ImGui::Separator();
    ImGui::BeginChild("control_body", ImVec2(0, 0), ImGuiChildFlags_None);
    ui_controls_body();
    ImGui::EndChild();
    ImGui::End();
}

static void ui_plane_screen() {
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(vp->WorkPos);
    ImGui::SetNextWindowSize(vp->WorkSize);
    ImGui::Begin("Plane", nullptr,
        ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoBringToFrontOnFocus |
        ImGuiWindowFlags_NoNavFocus);
    ImVec2 avail = ImGui::GetContentRegionAvail();
    if (avail.x < 16) avail.x = 16;
    if (avail.y < 16) avail.y = 16;
    ensure_fbo((int)avail.x, (int)avail.y);
    render_view();
    ImVec2 origin = ImGui::GetCursorScreenPos();
    ImGui::Image(ImTextureRef((ImTextureID)(ImU64)g.color), avail, ImVec2(0, 1), ImVec2(1, 0));
    bool hovered = ImGui::IsItemHovered();
    ImGuiIO& io = ImGui::GetIO();
    if (hovered) {
        double mx = io.MousePos.x - origin.x;
        double my = avail.y - (io.MousePos.y - origin.y);
        auto z = from_pixel(mx, my, avail.x, avail.y);
        g.probe_x = cx::real(z);
        g.probe_y = cx::imag(z);
        g.has_probe = true;
        if (io.MouseWheel != 0 && !io.KeyCtrl) {
            double factor = std::pow(1.12, -io.MouseWheel);
            g.center_x = z.real() + (g.center_x - z.real()) * factor;
            g.center_y = z.imag() + (g.center_y - z.imag()) * factor;
            g.scale *= factor;
            g.view_dirty = true;
        }
        if (ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
            double m = std::min(avail.x, avail.y);
            g.center_x -= io.MouseDelta.x / m * g.scale;
            g.center_y += io.MouseDelta.y / m * g.scale;
            g.view_dirty = true;
        }
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
            g.julia_re = g.probe_x;
            g.julia_im = g.probe_y;
        }
        if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left))
            g.pins.push_back({ z });
    }
    draw_overlay(origin, avail);
    ImGui::End();

    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + 10.0f, vp->WorkPos.y + 10.0f), ImGuiCond_Always);
    ImGui::SetNextWindowBgAlpha(0.78f);
    ImGui::Begin("##control_btn", nullptr,
        ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize |
        ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing |
        ImGuiWindowFlags_NoNav | ImGuiWindowFlags_NoMove);
    if (ImGui::Button("Control"))
        g.control_open = true;
    ImGui::End();
}

int main() {
    glfwSetErrorCallback([](int, const char* d) { std::fprintf(stderr, "%s\n", d); });
    if (!glfwInit()) return 1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GLFWwindow* win = glfwCreateWindow(1600, 900, "Complex-Geom  —  Ultra Fractal / GeoGebra", nullptr, nullptr);
    if (!win) return 1;
    glfwMakeContextCurrent(win);
    glfwSwapInterval(1);
    load_gl_pointers();

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.IniFilename = "imgui.ini";
    ImGui::StyleColorsDark();
    {
        float dpi = 1.0f;
        if (GLFWmonitor* mon = glfwGetPrimaryMonitor())
            dpi = ImGui_ImplGlfw_GetContentScaleForMonitor(mon);
        if (dpi < 1.0f) dpi = 1.0f;
        ImGuiStyle& style = ImGui::GetStyle();
        style.FontSizeBase = 16.0f;
        style.FontScaleMain = 1.0f;
        style.FontScaleDpi = dpi;
        style.FramePadding = ImVec2(8.0f, 5.0f);
        style.ItemSpacing = ImVec2(8.0f, 6.0f);
        style.ItemInnerSpacing = ImVec2(6.0f, 4.0f);
        style.WindowPadding = ImVec2(12.0f, 10.0f);
        style.GrabMinSize = 14.0f;
        style.ScrollbarSize = 16.0f;
        style.FrameRounding = 3.0f;
        style.GrabRounding = 3.0f;
        style.ScaleAllSizes(dpi);
    }
    ImGui_ImplGlfw_InitForOpenGL(win, true);
    ImGui_ImplOpenGL3_Init("#version 400");

    g.prog_domain = link_program(kVSQuad, kFSDomain);
    g.prog_fractal = link_program(kVSQuad, kFSFractal);
    g.fractal_double = g.prog_fractal != 0 && glUniform2d != nullptr;
    if (!g.prog_fractal) {
        g.prog_fractal = link_program(kVSQuad, kFSFractalF);
        g.fractal_double = false;
    }
    g.prog_points = link_program(kVSPoints, kFSPoints);
    if (!g.prog_domain || !g.prog_fractal || !g.prog_points) {
        std::fprintf(stderr, "required shaders failed to compile.\n");
        return 1;
    }
    glGenVertexArrays(1, &g.vao_quad);
    glGenVertexArrays(1, &g.vao_pts);
    glGenBuffers(1, &g.vbo_quad);
    glGenBuffers(1, &g.vbo_pts);
    float quad[8] = { -1, -1, 1, -1, -1, 1, 1, 1 };
    glBindVertexArray(g.vao_quad);
    glBindBuffer(GL_ARRAY_BUFFER, g.vbo_quad);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quad), quad, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glBindVertexArray(0);

    bool esc_down = false;

    while (!glfwWindowShouldClose(win)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        if (g.control_open)
            ui_control_screen();
        else
            ui_plane_screen();

        ImGui::Render();
        int dw, dh;
        glfwGetFramebufferSize(win, &dw, &dh);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glViewport(0, 0, dw, dh);
        glClearColor(0.06f, 0.06f, 0.07f, 1);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(win);

        bool esc = glfwGetKey(win, GLFW_KEY_ESCAPE) == GLFW_PRESS;
        if (esc && !esc_down) {
            if (g.control_open) close_controls();
            else break;
        }
        esc_down = esc;
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(win);
    glfwTerminate();
    return 0;
}
