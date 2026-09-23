# Complex-Geom

在复平面上查看函数，并按同一套公式做分形迭代。代数不限于普通复数：还可以切到分裂复数（双曲数）和对偶数。计算用自带的 [scomplex](ref/使用说明_中文.md)，画面用 OpenGL 画，界面用 Dear ImGui。

## 界面

启动后是整幅平面。左上角 **Control** 打开控制页，**Esc** 或 **Back to Plane** 回到平面。

| 操作 | 作用 |
|---|---|
| 滚轮 | 以光标为中心缩放 |
| 左键拖动 | 平移 |
| 左键双击 | 在该点钉一枚 pin |
| 右键 | 把该点写入 Julia 种子的实部、虚部 |

控制页里可以切换：

- **Function / Fractal**：看函数，或做迭代。
- **Algebra**：`i² = −1`（普通复数）、`j² = +1`（分裂复数）、`ε² = 0`（对偶数）。
- **f(z)**：手写公式，或点预设。预设包括 `z^2`、`z^2+c`、`exp(z)`、`gamma(z)` 等。
- **Calculator**：用当前代数对表达式求值，可把结果拿去当 `f`，或写成 Julia 的 `c`。

### 函数

三种画法：

- **Domain coloring**：按 `f(z)` 的辐角着色、按模长调亮度。普通复数下，`z`、`z^2`、`z^3`、`1/z`、`exp(z)`、`sin(z)`、`cos(z)`、`log(z)`、`sqrt(z)` 走 GPU；其余公式，以及分裂复数、对偶数，走 CPU。
- **Mapping**：采样格点，同时画出 `z` 和 `f(z)`。
- **Polya field**：沿 `conj(f) / |f|` 做 RK4，画出场线。步数和长度可调。

### 分形

每一步先按族处理 `z`，再套用 `f`：

| 族 | 含义 |
|---|---|
| **M** | 直接用 `f` |
| **T** | 先取共轭，再做 M |
| **B** | 先对两个分量取绝对值，再做 M |

公式里**没有** `c` 时，迭代是 `z → f(z) + c`。勾上 **L** 则变成 `z → z − f(z) + c`。公式里**已经有** `c` 时，整式就是迭代器，L 不起作用。

**Julia** 打开后，像素是 `z`，`c` 用面板里的种子；关掉时像素是 `c`，`z` 从 0 出发。

普通复数的分形在 GPU 上迭代：当前公式会编成着色器。加减乘用双精度；`exp` / `log` / `sin` 用显卡的 float 指令，避免每个像素跑长级数把窗口堵死。`0` 的正实部复幂按 `0` 算，NaN 和 Inf 算作逃逸。放到双精度装不下像素差时，画面不再变细，但程序仍会响应。`gamma`、分裂复数和对偶数仍在 CPU 上算，并用四叉树跳过整块都不逃逸的区域。

## 公式

不区分大小写的名字：`z` 或 `x` 是自变量，`c` 是参数，`i` 或 `j` 是虚单位，`pi`、`e` 是常数。运算符是 `+`、`-`、`*`、`/`、`^` 和括号。

一元函数：`sin` `cos` `tan` `cot` `sec` `csc` `exp` `log` `ln` `sqrt` `cbrt` `abs` `conj` `sinh` `cosh` `tanh` `asin` `acos` `atan` `asinh` `acosh` `atanh` `gamma` `re`（`real`）`im`（`imag`）。

二元：`log(z, base)`、`pow(z, w)`。`^` 与 `pow` 相同。

`abs` 是对两个分量分别取绝对值，不是模长。超越函数按当前代数用库自己的乘法实现，不是 `std::sin` / `std::exp`。细节、误差和三种语言移植见 `ref/`：

- [使用说明（中文）](ref/使用说明_中文.md)
- [introduction（English）](ref/introduction_en.md)
- [Mode d'emploi（français）](ref/Mode_d'emploi_français.md)
- [使用説明（日本語）](ref/使用説明_日本語.md)

## 编译与运行

在 Windows 上，用已安装 OpenMP 的 Visual Studio 2022 或 Visual Studio 18 的 x64 开发者环境：

```bat
build.bat
bin\Complex-Geom.exe
```

`build.bat` 会调用本机的 `vcvars64.bat`，用 MSVC `/std:c++17 /O2 /openmp` 编译，并静态链接 `third_party\glfw\lib-vc2022\glfw3.lib`。需要能创建 **OpenGL 4.0 core** 上下文的显卡驱动。

产物写到 `bin\`，该目录已被忽略。运行时会在当前工作目录读写 `imgui.ini`。

## 目录

```
scomplex.hpp / scomplex.cpp    程序实际链接的 C++ 库
src/main.cpp                   窗口、绘制、分形与 Pólya
src/expr.hpp                   公式解析
src/gl_api.hpp                 所需 OpenGL 入口
ref/                           C、Haskell 移植，以及四种语言的库说明
third_party/imgui              Dear ImGui 1.93.0 WIP
third_party/glfw               GLFW 3.4 头文件与预编译库
build.bat
```

C++ 类型在 `std::complex::experimental::zhangzl`。头文件不包含 `<complex>`，因此不和 `std::complex` 撞名。`ref/` 里的 C 与 Haskell 移植不参与这个可执行文件的链接，用来和 C++ 对照。

## 第三方

- [Dear ImGui](third_party/imgui/LICENSE.txt)，MIT 许可证。
- [GLFW 3.4](third_party/glfw/LICENSE.md)，zlib/libpng 许可证。
