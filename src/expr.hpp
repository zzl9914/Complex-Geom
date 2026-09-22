#pragma once
#include "../scomplex.hpp"
#include <cctype>
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

inline cx::complex<double> eval_compiled(const CompiledExpr& e, cx::complex<double> z,
                                         cx::complex<double> c = cx::complex<double>()) {
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
