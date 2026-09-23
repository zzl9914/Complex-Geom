#pragma once
#include "../scomplex.hpp"
#include <cctype>
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cx = std::complex::experimental::zhangzl;

struct ExprError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

enum class ExprFn : unsigned char {
    Sin, Cos, Tan, Cot, Sec, Csc, Exp, Log, Ln, Sqrt, Cbrt, Abs, Conj,
    Sinh, Cosh, Tanh, Asin, Acos, Atan, Asinh, Acosh, Atanh, Gamma, Pow, Re, Im
};

struct AstNode {
    enum Kind : unsigned char { Const, Z, C, Neg, Add, Sub, Mul, Div, Pow, Call } kind;
    unsigned char fn = 0;
    int a = -1;
    int b = -1;
    double re = 0;
    double im = 0;
};

struct CompiledExpr {
    std::vector<AstNode> n;
    int root = -1;
    bool uses_c = false;
};

// Ordinary complex (system == -1) uses principal exp/log formulas from <cmath>.
// scomplex's exp(log) path is an AGM / series that takes ~10 microseconds per
// call; a fractal frame then spends minutes inside non-escaping pixels.
// Split and dual numbers stay on scomplex, where those formulas do not apply.
namespace ordinary {
struct Z { double re, im; };
inline Z zadd(Z a, Z b) { return { a.re + b.re, a.im + b.im }; }
inline Z zsub(Z a, Z b) { return { a.re - b.re, a.im - b.im }; }
inline Z zneg(Z a) { return { -a.re, -a.im }; }
inline Z zmul(Z a, Z b) {
    return { a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re };
}
inline Z zdiv(Z a, Z b) {
    double d = b.re * b.re + b.im * b.im;
    if (!(d > 0.0)) throw ExprError("division by zero");
    return { (a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d };
}
inline double fin(double v) {
    if (!std::isfinite(v)) return 1e150;
    if (v > 1e150) return 1e150;
    if (v < -1e150) return -1e150;
    return v;
}
inline Z zpack(Z a) { return { fin(a.re), fin(a.im) }; }
inline Z zexp(Z a) {
    double e = std::exp(a.re);
    return zpack({ e * std::cos(a.im), e * std::sin(a.im) });
}
inline Z zlog(Z a) {
    double n2 = a.re * a.re + a.im * a.im;
    if (!(n2 > 0.0)) throw ExprError("log of zero");
    return { 0.5 * std::log(n2), std::atan2(a.im, a.re) };
}
inline Z zsqrt(Z a) {
    if (a.im == 0.0 && a.re >= 0.0) return { std::sqrt(a.re), 0.0 };
    return zexp(zmul(zlog(a), { 0.5, 0.0 }));
}
inline Z zpow(Z a, Z b) {
    if (b.re == 0.0 && b.im == 0.0) return { 1.0, 0.0 };
    if (a.re == 0.0 && a.im == 0.0) {
        if (b.im == 0.0 && b.re > 0.0) return { 0.0, 0.0 };
        throw ExprError("zero to a non-positive power");
    }
    if (b.im == 0.0) {
        double n = std::round(b.re);
        if (n == b.re && std::fabs(n) <= 4096.0) {
            int e = (int)n;
            bool neg = e < 0;
            unsigned u = neg ? (unsigned)(-e) : (unsigned)e;
            Z r{ 1.0, 0.0 }, base = a;
            while (u) {
                if (u & 1u) r = zmul(r, base);
                u >>= 1;
                if (u) base = zmul(base, base);
            }
            return neg ? zdiv(Z{ 1.0, 0.0 }, r) : zpack(r);
        }
    }
    return zexp(zmul(zlog(a), b));
}
inline Z zsin(Z a) {
    return zpack({ std::sin(a.re) * std::cosh(a.im), std::cos(a.re) * std::sinh(a.im) });
}
inline Z zcos(Z a) {
    return zpack({ std::cos(a.re) * std::cosh(a.im), -std::sin(a.re) * std::sinh(a.im) });
}
inline Z zsinh(Z a) {
    return zpack({ std::sinh(a.re) * std::cos(a.im), std::cosh(a.re) * std::sin(a.im) });
}
inline Z zcosh(Z a) {
    return zpack({ std::cosh(a.re) * std::cos(a.im), std::sinh(a.re) * std::sin(a.im) });
}
inline Z zasin(Z a) {
    Z w = zlog(zadd(zmul(Z{ 0.0, 1.0 }, a), zsqrt(zsub(Z{ 1.0, 0.0 }, zmul(a, a)))));
    return zpack({ w.im, -w.re });
}
inline Z zacos(Z a) { return zsub(Z{ 1.5707963267948966, 0.0 }, zasin(a)); }
inline Z zatan(Z a) {
    Z w = zlog(zdiv(zadd(Z{ 0.0, 1.0 }, a), zsub(Z{ 0.0, 1.0 }, a)));
    return zpack({ 0.5 * w.im, -0.5 * w.re });
}
inline Z zasinh(Z a) { return zlog(zadd(a, zsqrt(zadd(zmul(a, a), Z{ 1.0, 0.0 })))); }
inline Z zacosh(Z a) { return zlog(zadd(a, zsqrt(zsub(zmul(a, a), Z{ 1.0, 0.0 })))); }
inline Z zatanh(Z a) { return zmul(Z{ 0.5, 0.0 }, zlog(zdiv(zadd(Z{ 1.0, 0.0 }, a), zsub(Z{ 1.0, 0.0 }, a)))); }

inline cx::complex<double> eval(const CompiledExpr& e, cx::complex<double> z0, cx::complex<double> c0) {
    if (e.root < 0) throw ExprError("empty expression");
    auto go = [&](auto&& self, int i) -> Z {
        const AstNode& p = e.n[(size_t)i];
        switch (p.kind) {
        case AstNode::Const: return { p.re, p.im };
        case AstNode::Z: return { cx::real(z0), cx::imag(z0) };
        case AstNode::C: return { cx::real(c0), cx::imag(c0) };
        case AstNode::Neg: return zneg(self(self, p.a));
        case AstNode::Add: return zadd(self(self, p.a), self(self, p.b));
        case AstNode::Sub: return zsub(self(self, p.a), self(self, p.b));
        case AstNode::Mul: return zmul(self(self, p.a), self(self, p.b));
        case AstNode::Div: return zdiv(self(self, p.a), self(self, p.b));
        case AstNode::Pow: return zpow(self(self, p.a), self(self, p.b));
        case AstNode::Call: {
            Z x = self(self, p.a);
            Z y{ 0, 0 };
            if (p.b >= 0) y = self(self, p.b);
            switch ((ExprFn)p.fn) {
            case ExprFn::Sin: return zsin(x);
            case ExprFn::Cos: return zcos(x);
            case ExprFn::Tan: return zdiv(zsin(x), zcos(x));
            case ExprFn::Cot: return zdiv(zcos(x), zsin(x));
            case ExprFn::Sec: return zdiv(Z{ 1, 0 }, zcos(x));
            case ExprFn::Csc: return zdiv(Z{ 1, 0 }, zsin(x));
            case ExprFn::Exp: return zexp(x);
            case ExprFn::Log: return p.b >= 0 ? zdiv(zlog(x), zlog(y)) : zlog(x);
            case ExprFn::Ln: return zlog(x);
            case ExprFn::Sqrt: return zsqrt(x);
            case ExprFn::Cbrt:
                if (x.im == 0.0) return { std::cbrt(x.re), 0.0 };
                return zpow(x, Z{ 1.0 / 3.0, 0.0 });
            case ExprFn::Abs: return { std::fabs(x.re), std::fabs(x.im) };
            case ExprFn::Conj: return { x.re, -x.im };
            case ExprFn::Sinh: return zsinh(x);
            case ExprFn::Cosh: return zcosh(x);
            case ExprFn::Tanh: return zdiv(zsinh(x), zcosh(x));
            case ExprFn::Asin: return zasin(x);
            case ExprFn::Acos: return zacos(x);
            case ExprFn::Atan: return zatan(x);
            case ExprFn::Asinh: return zasinh(x);
            case ExprFn::Acosh: return zacosh(x);
            case ExprFn::Atanh: return zatanh(x);
            case ExprFn::Gamma: {
                auto g = cx::gamma(cx::complex<double>(x.re, x.im));
                return zpack({ cx::real(g), cx::imag(g) });
            }
            case ExprFn::Pow: return zpow(x, y);
            case ExprFn::Re: return { x.re, 0 };
            case ExprFn::Im: return { x.im, 0 };
            }
            throw ExprError("unknown function");
        }
        }
        throw ExprError("bad ast");
    };
    Z w = go(go, e.root);
    return cx::complex<double>(fin(w.re), fin(w.im));
}
} // namespace ordinary

inline cx::complex<double> eval_compiled(const CompiledExpr& e, cx::complex<double> z,
                                         cx::complex<double> c = cx::complex<double>()) {
    if (cx::complex<double>::get_system() == -1.0) return ordinary::eval(e, z, c);
    if (e.root < 0) throw ExprError("empty expression");
    auto go = [&](auto&& self, int i) -> cx::complex<double> {
        const AstNode& p = e.n[(size_t)i];
        switch (p.kind) {
        case AstNode::Const: return cx::complex<double>(p.re, p.im);
        case AstNode::Z: return z;
        case AstNode::C: return c;
        case AstNode::Neg: return -self(self, p.a);
        case AstNode::Add: return self(self, p.a) + self(self, p.b);
        case AstNode::Sub: return self(self, p.a) - self(self, p.b);
        case AstNode::Mul: return self(self, p.a) * self(self, p.b);
        case AstNode::Div: return self(self, p.a) / self(self, p.b);
        case AstNode::Pow: return cx::pow(self(self, p.a), self(self, p.b));
        case AstNode::Call: {
            auto x = self(self, p.a);
            cx::complex<double> y;
            if (p.b >= 0) y = self(self, p.b);
            switch ((ExprFn)p.fn) {
            case ExprFn::Sin: return cx::sin(x);
            case ExprFn::Cos: return cx::cos(x);
            case ExprFn::Tan: return cx::tan(x);
            case ExprFn::Cot: return cx::cot(x);
            case ExprFn::Sec: return cx::sec(x);
            case ExprFn::Csc: return cx::csc(x);
            case ExprFn::Exp: return cx::exp(x);
            case ExprFn::Log: return p.b >= 0 ? cx::log(x, y) : cx::log(x);
            case ExprFn::Ln: return cx::log(x);
            case ExprFn::Sqrt: return cx::sqrt(x);
            case ExprFn::Cbrt: return cx::cbrt(x);
            case ExprFn::Abs: return cx::abs(x);
            case ExprFn::Conj: return cx::conj(x);
            case ExprFn::Sinh: return cx::sinh(x);
            case ExprFn::Cosh: return cx::cosh(x);
            case ExprFn::Tanh: return cx::tanh(x);
            case ExprFn::Asin: return cx::asin(x);
            case ExprFn::Acos: return cx::acos(x);
            case ExprFn::Atan: return cx::atan(x);
            case ExprFn::Asinh: return cx::asinh(x);
            case ExprFn::Acosh: return cx::acosh(x);
            case ExprFn::Atanh: return cx::atanh(x);
            case ExprFn::Gamma: return cx::gamma(x);
            case ExprFn::Pow: return cx::pow(x, y);
            case ExprFn::Re: return cx::complex<double>(cx::real(x), 0.0);
            case ExprFn::Im: return cx::complex<double>(cx::imag(x), 0.0);
            }
            throw ExprError("unknown function");
        }
        }
        throw ExprError("bad ast");
    };
    return go(go, e.root);
}

class ExprCompiler {
public:
    CompiledExpr compile(const std::string& src) {
        s = src;
        i = 0;
        out.n.clear();
        out.root = -1;
        out.uses_c = false;
        skip();
        if (s.empty()) throw ExprError("empty expression");
        out.root = parse_add();
        skip();
        if (i < s.size()) throw ExprError("trailing input at '" + s.substr(i) + "'");
        return out;
    }

private:
    std::string s;
    size_t i = 0;
    CompiledExpr out;

    int push(AstNode n) {
        out.n.push_back(n);
        return (int)out.n.size() - 1;
    }
    void skip() {
        while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    }
    bool eat(char c) {
        skip();
        if (i < s.size() && s[i] == c) {
            ++i;
            return true;
        }
        return false;
    }
    bool starts_ident() const {
        if (i >= s.size()) return false;
        unsigned char c = static_cast<unsigned char>(s[i]);
        return std::isalpha(c) || c == '_';
    }
    std::string ident() {
        skip();
        size_t a = i;
        while (i < s.size()) {
            unsigned char c = static_cast<unsigned char>(s[i]);
            if (!std::isalnum(c) && c != '_') break;
            ++i;
        }
        return s.substr(a, i - a);
    }
    int bin(AstNode::Kind k, int a, int b) {
        AstNode n;
        n.kind = k;
        n.a = a;
        n.b = b;
        return push(n);
    }

    int parse_add() {
        int v = parse_mul();
        for (;;) {
            if (eat('+')) v = bin(AstNode::Add, v, parse_mul());
            else if (eat('-')) v = bin(AstNode::Sub, v, parse_mul());
            else break;
        }
        return v;
    }
    int parse_mul() {
        int v = parse_pow();
        for (;;) {
            skip();
            if (eat('*')) v = bin(AstNode::Mul, v, parse_pow());
            else if (eat('/')) v = bin(AstNode::Div, v, parse_pow());
            else if (i < s.size() && (s[i] == '(' || starts_ident() ||
                                      std::isdigit(static_cast<unsigned char>(s[i])) || s[i] == '.')) {
                v = bin(AstNode::Mul, v, parse_pow());
            } else break;
        }
        return v;
    }
    int parse_pow() {
        int v = parse_unary();
        skip();
        if (eat('^')) v = bin(AstNode::Pow, v, parse_pow());
        return v;
    }
    int parse_unary() {
        skip();
        if (eat('+')) return parse_unary();
        if (eat('-')) {
            AstNode n;
            n.kind = AstNode::Neg;
            n.a = parse_unary();
            return push(n);
        }
        return parse_primary();
    }
    int parse_primary() {
        skip();
        if (eat('(')) {
            int v = parse_add();
            if (!eat(')')) throw ExprError("missing ')'");
            return v;
        }
        if (i < s.size() && (std::isdigit(static_cast<unsigned char>(s[i])) || s[i] == '.')) {
            size_t a = i;
            while (i < s.size() && (std::isdigit(static_cast<unsigned char>(s[i])) || s[i] == '.' ||
                                    s[i] == 'e' || s[i] == 'E' || s[i] == '+' || s[i] == '-')) {
                if ((s[i] == '+' || s[i] == '-') && i > a && s[i - 1] != 'e' && s[i - 1] != 'E') break;
                ++i;
            }
            AstNode n;
            n.kind = AstNode::Const;
            n.re = std::stod(s.substr(a, i - a));
            return push(n);
        }
        if (starts_ident()) {
            std::string id = ident();
            if (id == "i" || id == "I" || id == "j" || id == "J") {
                AstNode n;
                n.kind = AstNode::Const;
                n.im = 1.0;
                return push(n);
            }
            if (id == "z" || id == "Z" || id == "x" || id == "X") {
                AstNode n;
                n.kind = AstNode::Z;
                return push(n);
            }
            if (id == "c" || id == "C") {
                out.uses_c = true;
                AstNode n;
                n.kind = AstNode::C;
                return push(n);
            }
            if (id == "pi" || id == "PI") {
                AstNode n;
                n.kind = AstNode::Const;
                n.re = cx::complex<double>::pi();
                return push(n);
            }
            if (id == "e") {
                AstNode n;
                n.kind = AstNode::Const;
                n.re = cx::complex<double>::e();
                return push(n);
            }
            if (!eat('(')) throw ExprError("expected '(' after " + id);
            int a = parse_add();
            int b = -1;
            if (eat(',')) b = parse_add();
            if (!eat(')')) throw ExprError("missing ')' after " + id);
            AstNode n;
            n.kind = AstNode::Call;
            n.a = a;
            n.b = b;
            n.fn = (unsigned char)fn_id(id, b >= 0);
            return push(n);
        }
        throw ExprError("unexpected token at '" + s.substr(i) + "'");
    }

    static ExprFn fn_id(const std::string& id, bool two) {
        if (id == "sin") return ExprFn::Sin;
        if (id == "cos") return ExprFn::Cos;
        if (id == "tan") return ExprFn::Tan;
        if (id == "cot") return ExprFn::Cot;
        if (id == "sec") return ExprFn::Sec;
        if (id == "csc") return ExprFn::Csc;
        if (id == "exp") return ExprFn::Exp;
        if (id == "log") return ExprFn::Log;
        if (id == "ln") return ExprFn::Ln;
        if (id == "sqrt") return ExprFn::Sqrt;
        if (id == "cbrt") return ExprFn::Cbrt;
        if (id == "abs") return ExprFn::Abs;
        if (id == "conj") return ExprFn::Conj;
        if (id == "sinh") return ExprFn::Sinh;
        if (id == "cosh") return ExprFn::Cosh;
        if (id == "tanh") return ExprFn::Tanh;
        if (id == "asin") return ExprFn::Asin;
        if (id == "acos") return ExprFn::Acos;
        if (id == "atan") return ExprFn::Atan;
        if (id == "asinh") return ExprFn::Asinh;
        if (id == "acosh") return ExprFn::Acosh;
        if (id == "atanh") return ExprFn::Atanh;
        if (id == "gamma") return ExprFn::Gamma;
        if (id == "pow") {
            if (!two) throw ExprError("pow needs two arguments");
            return ExprFn::Pow;
        }
        if (id == "re" || id == "real") return ExprFn::Re;
        if (id == "im" || id == "imag") return ExprFn::Im;
        throw ExprError("unknown function '" + id + "'");
    }
};

inline CompiledExpr compile_expr(const std::string& src) {
    return ExprCompiler().compile(src);
}

class ComplexExpr {
public:
    cx::complex<double> eval(const std::string& src, cx::complex<double> z,
                             cx::complex<double> c = cx::complex<double>()) {
        return eval_compiled(compile_expr(src), z, c);
    }
};

inline bool expr_has_gamma(const CompiledExpr& e) {
    for (const AstNode& n : e.n)
        if (n.kind == AstNode::Call && (ExprFn)n.fn == ExprFn::Gamma) return true;
    return false;
}

inline bool expr_uses_transcendental(const CompiledExpr& e) {
    for (const AstNode& n : e.n) {
        if (n.kind == AstNode::Pow) {
            const AstNode* b = (n.b >= 0 && n.b < (int)e.n.size()) ? &e.n[(size_t)n.b] : nullptr;
            double r = 0;
            bool ipow = b && b->kind == AstNode::Const && b->im == 0.0 &&
                        (r = std::round(b->re)) == b->re && std::fabs(r) <= 16.0;
            if (!ipow) return true;
        }
        if (n.kind != AstNode::Call) continue;
        switch ((ExprFn)n.fn) {
        case ExprFn::Abs: case ExprFn::Conj: case ExprFn::Re: case ExprFn::Im: break;
        default: return true;
        }
    }
    return false;
}

inline std::string glsl_number(double v, bool dbl) {
    if (!std::isfinite(v)) v = 0.0;
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    std::string s = buf;
    if (s.find('.') == std::string::npos && s.find('e') == std::string::npos && s.find('E') == std::string::npos)
        s += ".0";
    if (dbl) s += "lf";
    return s;
}

// GLSL for one ordinary-complex formula. Empty if it cannot run on the GPU
// (gamma, or an expression that did not compile).
inline std::string glsl_formula(const CompiledExpr& e, bool dbl) {
    if (e.root < 0 || expr_has_gamma(e)) return {};
    const char* T = dbl ? "dvec2" : "vec2";
    const char* S = dbl ? "double" : "float";
    auto num = [&](double v) { return glsl_number(v, dbl); };
    auto C = [&](double re, double im) {
        return std::string(T) + "(" + num(re) + "," + num(im) + ")";
    };
    std::ostringstream o;
    // Desktop GLSL 400 has float exp/log/sin, not double. Transcendentals go
    // through float; adds and multiplies stay in the shader's own type, so
    // z^2 can still iterate in double.
    auto widen = [&](const std::string& floatExpr) {
        return dbl ? "double(" + floatExpr + ")" : floatExpr;
    };
    o << "float sinh_r(float x){ return (exp(x)-exp(-x))*0.5; }\n";
    o << "float cosh_r(float x){ return (exp(x)+exp(-x))*0.5; }\n";
    o << T << " cmul(" << T << " a," << T << " b){ return " << T << "(a.x*b.x-a.y*b.y, a.x*b.y+a.y*b.x); }\n";
    o << T << " cdiv(" << T << " a," << T << " b){\n";
    o << "  " << S << " d=b.x*b.x+b.y*b.y;\n";
    o << "  if(!(d>" << num(0) << ")) return " << C(1e20, 0) << ";\n";
    o << "  return " << T << "(a.x*b.x+a.y*b.y, a.y*b.x-a.x*b.y)/d;\n}\n";
    o << T << " cexp(" << T << " a){ float y=float(a.y); float e=exp(float(a.x)); return "
      << T << "(" << widen("e*cos(y)") << "," << widen("e*sin(y)") << "); }\n";
    o << T << " clog(" << T << " a){\n";
    o << "  float n2=float(a.x)*float(a.x)+float(a.y)*float(a.y);\n";
    o << "  if(!(n2>0.0)) return " << C(-1e8, 0) << ";\n";
    o << "  return " << T << "(" << widen("0.5*log(n2)") << "," << widen("atan(float(a.y),float(a.x))") << ");\n}\n";
    o << T << " cpow(" << T << " a," << T << " b){\n";
    o << "  if(b.x==" << num(0) << " && b.y==" << num(0) << ") return " << C(1, 0) << ";\n";
    o << "  " << S << " aa=a.x*a.x+a.y*a.y;\n";
    o << "  if(!(aa>" << num(0) << ")){\n";
    o << "    if(b.y==" << num(0) << " && b.x>" << num(0) << ") return " << C(0, 0) << ";\n";
    o << "    return " << C(1e20, 0) << ";\n  }\n";
    o << "  return cexp(cmul(clog(a), b));\n}\n";
    o << T << " csqrt(" << T << " a){\n";
    o << "  if(a.y==" << num(0) << " && a.x>=" << num(0) << ") return " << T << "("
      << widen("sqrt(float(a.x))") << "," << num(0) << ");\n";
    o << "  return cpow(a," << C(0.5, 0) << ");\n}\n";
    o << T << " csin(" << T << " a){ float x=float(a.x), y=float(a.y); return " << T << "("
      << widen("sin(x)*cosh_r(y)") << "," << widen("cos(x)*sinh_r(y)") << "); }\n";
    o << T << " ccos(" << T << " a){ float x=float(a.x), y=float(a.y); return " << T << "("
      << widen("cos(x)*cosh_r(y)") << "," << widen("-sin(x)*sinh_r(y)") << "); }\n";
    o << T << " csinh(" << T << " a){ float x=float(a.x), y=float(a.y); return " << T << "("
      << widen("sinh_r(x)*cos(y)") << "," << widen("cosh_r(x)*sin(y)") << "); }\n";
    o << T << " ccosh(" << T << " a){ float x=float(a.x), y=float(a.y); return " << T << "("
      << widen("cosh_r(x)*cos(y)") << "," << widen("sinh_r(x)*sin(y)") << "); }\n";
    o << T << " casin(" << T << " a){\n";
    o << "  " << T << " w=clog(cmul(" << C(0, 1) << ",a)+csqrt(" << C(1, 0) << "-cmul(a,a)));\n";
    o << "  return " << T << "(w.y,-w.x);\n}\n";
    o << T << " cacos(" << T << " a){ return " << C(1.5707963267948966, 0) << "-casin(a); }\n";
    o << T << " catan(" << T << " a){\n";
    o << "  " << T << " w=clog(cdiv(" << C(0, 1) << "+a," << C(0, 1) << "-a));\n";
    o << "  return " << T << "(" << widen("0.5*float(w.y)") << "," << widen("-0.5*float(w.x)") << ");\n}\n";
    o << T << " casinh(" << T << " a){ return clog(a+csqrt(cmul(a,a)+" << C(1, 0) << ")); }\n";
    o << T << " cacosh(" << T << " a){ return clog(a+csqrt(cmul(a,a)-" << C(1, 0) << ")); }\n";
    o << T << " catanh(" << T << " a){ return cmul(" << C(0.5, 0) << ",clog(cdiv(" << C(1, 0) << "+a," << C(1, 0) << "-a))); }\n";
    o << T << " cf_f(" << T << " z," << T << " c){\n";
    for (int i = 0; i < (int)e.n.size(); ++i) {
        const AstNode& p = e.n[(size_t)i];
        auto vn = [](int k) { return "v" + std::to_string(k); };
        o << "  " << T << " " << vn(i) << "=";
        switch (p.kind) {
        case AstNode::Const: o << C(p.re, p.im); break;
        case AstNode::Z: o << "z"; break;
        case AstNode::C: o << "c"; break;
        case AstNode::Neg: o << T << "(-" << vn(p.a) << ".x,-" << vn(p.a) << ".y)"; break;
        case AstNode::Add: o << "(" << vn(p.a) << "+" << vn(p.b) << ")"; break;
        case AstNode::Sub: o << "(" << vn(p.a) << "-" << vn(p.b) << ")"; break;
        case AstNode::Mul: o << "cmul(" << vn(p.a) << "," << vn(p.b) << ")"; break;
        case AstNode::Div: o << "cdiv(" << vn(p.a) << "," << vn(p.b) << ")"; break;
        case AstNode::Pow: {
            const AstNode* eb = (p.b >= 0 && p.b < (int)e.n.size()) ? &e.n[(size_t)p.b] : nullptr;
            double n = 0;
            bool ipow = eb && eb->kind == AstNode::Const && eb->im == 0.0 &&
                        (n = std::round(eb->re)) == eb->re && std::fabs(n) <= 16.0;
            if (!ipow) {
                o << "cpow(" << vn(p.a) << "," << vn(p.b) << ")";
                break;
            }
            int ei = (int)n;
            if (ei == 0) { o << C(1, 0); break; }
            if (ei == 1) { o << vn(p.a); break; }
            if (ei == -1) { o << "cdiv(" << C(1, 0) << "," << vn(p.a) << ")"; break; }
            bool neg = ei < 0;
            int k = neg ? -ei : ei;
            std::string r = vn(p.a);
            for (int s = 1; s < k; ++s) r = "cmul(" + r + "," + vn(p.a) + ")";
            if (neg) r = "cdiv(" + C(1, 0) + "," + r + ")";
            o << r;
            break;
        }
        case AstNode::Call: {
            std::string a = vn(p.a);
            std::string b = p.b >= 0 ? vn(p.b) : std::string();
            switch ((ExprFn)p.fn) {
            case ExprFn::Sin: o << "csin(" << a << ")"; break;
            case ExprFn::Cos: o << "ccos(" << a << ")"; break;
            case ExprFn::Tan: o << "cdiv(csin(" << a << "),ccos(" << a << "))"; break;
            case ExprFn::Cot: o << "cdiv(ccos(" << a << "),csin(" << a << "))"; break;
            case ExprFn::Sec: o << "cdiv(" << C(1, 0) << ",ccos(" << a << "))"; break;
            case ExprFn::Csc: o << "cdiv(" << C(1, 0) << ",csin(" << a << "))"; break;
            case ExprFn::Exp: o << "cexp(" << a << ")"; break;
            case ExprFn::Log:
                if (p.b >= 0) o << "cdiv(clog(" << a << "),clog(" << b << "))";
                else o << "clog(" << a << ")";
                break;
            case ExprFn::Ln: o << "clog(" << a << ")"; break;
            case ExprFn::Sqrt: o << "csqrt(" << a << ")"; break;
            case ExprFn::Cbrt: o << "cpow(" << a << "," << C(1.0 / 3.0, 0) << ")"; break;
            case ExprFn::Abs: o << T << "(abs(" << a << ".x),abs(" << a << ".y))"; break;
            case ExprFn::Conj: o << T << "(" << a << ".x,-" << a << ".y)"; break;
            case ExprFn::Sinh: o << "csinh(" << a << ")"; break;
            case ExprFn::Cosh: o << "ccosh(" << a << ")"; break;
            case ExprFn::Tanh: o << "cdiv(csinh(" << a << "),ccosh(" << a << "))"; break;
            case ExprFn::Asin: o << "casin(" << a << ")"; break;
            case ExprFn::Acos: o << "cacos(" << a << ")"; break;
            case ExprFn::Atan: o << "catan(" << a << ")"; break;
            case ExprFn::Asinh: o << "casinh(" << a << ")"; break;
            case ExprFn::Acosh: o << "cacosh(" << a << ")"; break;
            case ExprFn::Atanh: o << "catanh(" << a << ")"; break;
            case ExprFn::Gamma: return {};
            case ExprFn::Pow: o << "cpow(" << a << "," << b << ")"; break;
            case ExprFn::Re: o << T << "(" << a << ".x," << num(0) << ")"; break;
            case ExprFn::Im: o << T << "(" << a << ".y," << num(0) << ")"; break;
            }
            break;
        }
        }
        o << ";\n";
    }
    o << "  return v" << e.root << ";\n}\n";
    return o.str();
}
