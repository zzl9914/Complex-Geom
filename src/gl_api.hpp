#pragma once
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#ifndef GLFW_INCLUDE_NONE
#define GLFW_INCLUDE_NONE
#endif
#include <GLFW/glfw3.h>

typedef unsigned int GLenum;
typedef unsigned int GLuint;
typedef int GLint;
typedef int GLsizei;
typedef unsigned char GLboolean;
typedef float GLfloat;
typedef double GLdouble;
typedef char GLchar;
typedef unsigned char GLubyte;
typedef std::ptrdiff_t GLsizeiptr;
typedef std::ptrdiff_t GLintptr;
typedef unsigned int GLbitfield;
typedef float GLclampf;
typedef double GLclampd;
typedef void GLvoid;

#define GL_FALSE 0
#define GL_TRUE 1
#define GL_COLOR_BUFFER_BIT 0x00004000
#define GL_POINTS 0x0000
#define GL_LINES 0x0001
#define GL_TRIANGLE_STRIP 0x0005
#define GL_ARRAY_BUFFER 0x8892
#define GL_STATIC_DRAW 0x88E4
#define GL_DYNAMIC_DRAW 0x88E8
#define GL_FLOAT 0x1406
#define GL_UNSIGNED_BYTE 0x1401
#define GL_RGB 0x1907
#define GL_RGBA 0x1908
#define GL_TEXTURE_2D 0x0DE1
#define GL_TEXTURE_MIN_FILTER 0x2801
#define GL_TEXTURE_MAG_FILTER 0x2800
#define GL_TEXTURE_WRAP_S 0x2802
#define GL_TEXTURE_WRAP_T 0x2803
#define GL_UNPACK_ALIGNMENT 0x0CF5
#define GL_LINEAR 0x2601
#define GL_CLAMP_TO_EDGE 0x812F
#define GL_RGB8 0x8051
#define GL_FRAMEBUFFER 0x8D40
#define GL_COLOR_ATTACHMENT0 0x8CE0
#define GL_FRAMEBUFFER_COMPLETE 0x8CD5
#define GL_VERTEX_SHADER 0x8B31
#define GL_FRAGMENT_SHADER 0x8B30
#define GL_COMPILE_STATUS 0x8B81
#define GL_LINK_STATUS 0x8B82
#define GL_INFO_LOG_LENGTH 0x8B84
#define GL_BLEND 0x0BE2
#define GL_SRC_ALPHA 0x0302
#define GL_ONE_MINUS_SRC_ALPHA 0x0303
#define GL_PROGRAM_POINT_SIZE 0x8642
#define GL_TEXTURE0 0x84C0
#define GL_DEPTH_TEST 0x0B71

#ifndef APIENTRY
#ifdef _WIN32
#define APIENTRY __stdcall
#else
#define APIENTRY
#endif
#endif
#ifndef APIENTRYP
#define APIENTRYP APIENTRY *
#endif

typedef GLuint(APIENTRYP PFNGLCREATESHADERPROC)(GLenum);
typedef void(APIENTRYP PFNGLSHADERSOURCEPROC)(GLuint, GLsizei, const GLchar* const*, const GLint*);
typedef void(APIENTRYP PFNGLCOMPILESHADERPROC)(GLuint);
typedef void(APIENTRYP PFNGLGETSHADERIVPROC)(GLuint, GLenum, GLint*);
typedef void(APIENTRYP PFNGLGETSHADERINFOLOGPROC)(GLuint, GLsizei, GLsizei*, GLchar*);
typedef void(APIENTRYP PFNGLDELETESHADERPROC)(GLuint);
typedef GLuint(APIENTRYP PFNGLCREATEPROGRAMPROC)(void);
typedef void(APIENTRYP PFNGLATTACHSHADERPROC)(GLuint, GLuint);
typedef void(APIENTRYP PFNGLLINKPROGRAMPROC)(GLuint);
typedef void(APIENTRYP PFNGLGETPROGRAMIVPROC)(GLuint, GLenum, GLint*);
typedef void(APIENTRYP PFNGLGETPROGRAMINFOLOGPROC)(GLuint, GLsizei, GLsizei*, GLchar*);
typedef void(APIENTRYP PFNGLUSEPROGRAMPROC)(GLuint);
typedef void(APIENTRYP PFNGLDELETEPROGRAMPROC)(GLuint);
typedef GLint(APIENTRYP PFNGLGETUNIFORMLOCATIONPROC)(GLuint, const GLchar*);
typedef void(APIENTRYP PFNGLUNIFORM1IPROC)(GLint, GLint);
typedef void(APIENTRYP PFNGLUNIFORM1FPROC)(GLint, GLfloat);
typedef void(APIENTRYP PFNGLUNIFORM2FPROC)(GLint, GLfloat, GLfloat);
typedef void(APIENTRYP PFNGLUNIFORM1DPROC)(GLint, GLdouble);
typedef void(APIENTRYP PFNGLUNIFORM2DPROC)(GLint, GLdouble, GLdouble);
typedef void(APIENTRYP PFNGLGENVERTEXARRAYSPROC)(GLsizei, GLuint*);
typedef void(APIENTRYP PFNGLBINDVERTEXARRAYPROC)(GLuint);
typedef void(APIENTRYP PFNGLDELETEVERTEXARRAYSPROC)(GLsizei, const GLuint*);
typedef void(APIENTRYP PFNGLGENBUFFERSPROC)(GLsizei, GLuint*);
typedef void(APIENTRYP PFNGLBINDBUFFERPROC)(GLenum, GLuint);
typedef void(APIENTRYP PFNGLBUFFERDATAPROC)(GLenum, GLsizeiptr, const void*, GLenum);
typedef void(APIENTRYP PFNGLDELETEBUFFERSPROC)(GLsizei, const GLuint*);
typedef void(APIENTRYP PFNGLENABLEVERTEXATTRIBARRAYPROC)(GLuint);
typedef void(APIENTRYP PFNGLVERTEXATTRIBPOINTERPROC)(GLuint, GLint, GLenum, GLboolean, GLsizei, const void*);
typedef void(APIENTRYP PFNGLGENTEXTURESPROC)(GLsizei, GLuint*);
typedef void(APIENTRYP PFNGLBINDTEXTUREPROC)(GLenum, GLuint);
typedef void(APIENTRYP PFNGLTEXIMAGE2DPROC)(GLenum, GLint, GLint, GLsizei, GLsizei, GLint, GLenum, GLenum, const void*);
typedef void(APIENTRYP PFNGLTEXPARAMETERIPROC)(GLenum, GLenum, GLint);
typedef void(APIENTRYP PFNGLACTIVETEXTUREPROC)(GLenum);
typedef void(APIENTRYP PFNGLDRAWARRAYSPROC)(GLenum, GLint, GLsizei);
typedef void(APIENTRYP PFNGLENABLEPROC)(GLenum);
typedef void(APIENTRYP PFNGLDISABLEPROC)(GLenum);
typedef void(APIENTRYP PFNGLBLENDFUNCPROC)(GLenum, GLenum);
typedef void(APIENTRYP PFNGLVIEWPORTPROC)(GLint, GLint, GLsizei, GLsizei);
typedef void(APIENTRYP PFNGLCLEARPROC)(GLbitfield);
typedef void(APIENTRYP PFNGLCLEARCOLORPROC)(GLfloat, GLfloat, GLfloat, GLfloat);
typedef void(APIENTRYP PFNGLGENFRAMEBUFFERSPROC)(GLsizei, GLuint*);
typedef void(APIENTRYP PFNGLBINDFRAMEBUFFERPROC)(GLenum, GLuint);
typedef void(APIENTRYP PFNGLFRAMEBUFFERTEXTURE2DPROC)(GLenum, GLenum, GLenum, GLuint, GLint);
typedef GLenum(APIENTRYP PFNGLCHECKFRAMEBUFFERSTATUSPROC)(GLenum);
typedef void(APIENTRYP PFNGLDELETEFRAMEBUFFERSPROC)(GLsizei, const GLuint*);
typedef void(APIENTRYP PFNGLDELETETEXTURESPROC)(GLsizei, const GLuint*);
typedef void(APIENTRYP PFNGLPIXELSTOREIPROC)(GLenum, GLint);

inline PFNGLCREATESHADERPROC glCreateShader;
inline PFNGLSHADERSOURCEPROC glShaderSource;
inline PFNGLCOMPILESHADERPROC glCompileShader;
inline PFNGLGETSHADERIVPROC glGetShaderiv;
inline PFNGLGETSHADERINFOLOGPROC glGetShaderInfoLog;
inline PFNGLDELETESHADERPROC glDeleteShader;
inline PFNGLCREATEPROGRAMPROC glCreateProgram;
inline PFNGLATTACHSHADERPROC glAttachShader;
inline PFNGLLINKPROGRAMPROC glLinkProgram;
inline PFNGLGETPROGRAMIVPROC glGetProgramiv;
inline PFNGLGETPROGRAMINFOLOGPROC glGetProgramInfoLog;
inline PFNGLUSEPROGRAMPROC glUseProgram;
inline PFNGLDELETEPROGRAMPROC glDeleteProgram;
inline PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation;
inline PFNGLUNIFORM1IPROC glUniform1i;
inline PFNGLUNIFORM1FPROC glUniform1f;
inline PFNGLUNIFORM2FPROC glUniform2f;
inline PFNGLUNIFORM1DPROC glUniform1d;
inline PFNGLUNIFORM2DPROC glUniform2d;
inline PFNGLGENVERTEXARRAYSPROC glGenVertexArrays;
inline PFNGLBINDVERTEXARRAYPROC glBindVertexArray;
inline PFNGLDELETEVERTEXARRAYSPROC glDeleteVertexArrays;
inline PFNGLGENBUFFERSPROC glGenBuffers;
inline PFNGLBINDBUFFERPROC glBindBuffer;
inline PFNGLBUFFERDATAPROC glBufferData;
inline PFNGLDELETEBUFFERSPROC glDeleteBuffers;
inline PFNGLENABLEVERTEXATTRIBARRAYPROC glEnableVertexAttribArray;
inline PFNGLVERTEXATTRIBPOINTERPROC glVertexAttribPointer;
inline PFNGLGENTEXTURESPROC glGenTextures;
inline PFNGLBINDTEXTUREPROC glBindTexture;
inline PFNGLTEXIMAGE2DPROC glTexImage2D;
inline PFNGLTEXPARAMETERIPROC glTexParameteri;
inline PFNGLACTIVETEXTUREPROC glActiveTexture;
inline PFNGLDRAWARRAYSPROC glDrawArrays;
inline PFNGLENABLEPROC glEnable;
inline PFNGLDISABLEPROC glDisable;
inline PFNGLBLENDFUNCPROC glBlendFunc;
inline PFNGLVIEWPORTPROC glViewport;
inline PFNGLCLEARPROC glClear;
inline PFNGLCLEARCOLORPROC glClearColor;
inline PFNGLGENFRAMEBUFFERSPROC glGenFramebuffers;
inline PFNGLBINDFRAMEBUFFERPROC glBindFramebuffer;
inline PFNGLFRAMEBUFFERTEXTURE2DPROC glFramebufferTexture2D;
inline PFNGLCHECKFRAMEBUFFERSTATUSPROC glCheckFramebufferStatus;
inline PFNGLDELETEFRAMEBUFFERSPROC glDeleteFramebuffers;
inline PFNGLDELETETEXTURESPROC glDeleteTextures;
inline PFNGLPIXELSTOREIPROC glPixelStorei;

#define LOAD_GL(name) name = (decltype(name))glfwGetProcAddress(#name)

inline void load_gl_pointers() {
    LOAD_GL(glCreateShader);
    LOAD_GL(glShaderSource);
    LOAD_GL(glCompileShader);
    LOAD_GL(glGetShaderiv);
    LOAD_GL(glGetShaderInfoLog);
    LOAD_GL(glDeleteShader);
    LOAD_GL(glCreateProgram);
    LOAD_GL(glAttachShader);
    LOAD_GL(glLinkProgram);
    LOAD_GL(glGetProgramiv);
    LOAD_GL(glGetProgramInfoLog);
    LOAD_GL(glUseProgram);
    LOAD_GL(glDeleteProgram);
    LOAD_GL(glGetUniformLocation);
    LOAD_GL(glUniform1i);
    LOAD_GL(glUniform1f);
    LOAD_GL(glUniform2f);
    LOAD_GL(glUniform1d);
    LOAD_GL(glUniform2d);
    LOAD_GL(glGenVertexArrays);
    LOAD_GL(glBindVertexArray);
    LOAD_GL(glDeleteVertexArrays);
    LOAD_GL(glGenBuffers);
    LOAD_GL(glBindBuffer);
    LOAD_GL(glBufferData);
    LOAD_GL(glDeleteBuffers);
    LOAD_GL(glEnableVertexAttribArray);
    LOAD_GL(glVertexAttribPointer);
    LOAD_GL(glGenTextures);
    LOAD_GL(glBindTexture);
    LOAD_GL(glTexImage2D);
    LOAD_GL(glTexParameteri);
    LOAD_GL(glActiveTexture);
    LOAD_GL(glDrawArrays);
    LOAD_GL(glEnable);
    LOAD_GL(glDisable);
    LOAD_GL(glBlendFunc);
    LOAD_GL(glViewport);
    LOAD_GL(glClear);
    LOAD_GL(glClearColor);
    LOAD_GL(glGenFramebuffers);
    LOAD_GL(glBindFramebuffer);
    LOAD_GL(glFramebufferTexture2D);
    LOAD_GL(glCheckFramebufferStatus);
    LOAD_GL(glDeleteFramebuffers);
    LOAD_GL(glDeleteTextures);
    LOAD_GL(glPixelStorei);
    if (!glCreateShader) {
        std::fprintf(stderr, "OpenGL 4.0 entry points missing.\n");
        std::exit(1);
    }
}

inline GLuint compile_shader(GLenum type, const char* src) {
    GLuint s = glCreateShader(type);
    glShaderSource(s, 1, &src, nullptr);
    glCompileShader(s);
    GLint ok = 0;
    glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[4096];
        glGetShaderInfoLog(s, 4096, nullptr, log);
        std::fprintf(stderr, "shader compile:\n%s\n", log);
        glDeleteShader(s);
        return 0;
    }
    return s;
}

inline GLuint link_program(const char* vs, const char* fs) {
    GLuint v = compile_shader(GL_VERTEX_SHADER, vs);
    GLuint f = compile_shader(GL_FRAGMENT_SHADER, fs);
    if (!v || !f) {
        if (v) glDeleteShader(v);
        if (f) glDeleteShader(f);
        return 0;
    }
    GLuint p = glCreateProgram();
    glAttachShader(p, v);
    glAttachShader(p, f);
    glLinkProgram(p);
    GLint ok = 0;
    glGetProgramiv(p, GL_LINK_STATUS, &ok);
    glDeleteShader(v);
    glDeleteShader(f);
    if (!ok) {
        char log[4096];
        glGetProgramInfoLog(p, 4096, nullptr, log);
        std::fprintf(stderr, "program link:\n%s\n", log);
        glDeleteProgram(p);
        return 0;
    }
    return p;
}
