#define GLM_ENABLE_EXPERIMENTAL
#define STB_IMAGE_IMPLEMENTATION
#include "nans.cpp"

#include "stb_image.h"
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <SDL2/SDL.h>
#include <dlfcn.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdio>
#include <unistd.h>
#include <vector>
#include <assimp/cimport.h>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <map>

#include "sdl_nans.h"
#include <ft2build.h>
#include FT_FREETYPE_H

global_variable uint64 GlobalPerfCountFrequency;
global_variable std::map<uint8, character> Characters;


global_variable char *texture_paths[] =
{
    "../res/texture/container.jpg",
    "../res/texture/earth.jpg",
    "../res/texture/checkerboard.png",
    "../res/texture/Plastic_4K_Diffuse.jpg",
    ""
};

global_variable char *vs_paths[] =
    {
        "../shaders/cube.vs",
        "../shaders/text.vs",
        "../shaders/light.vs",
        ""
    };

global_variable char *fs_paths[] =
    {
        "../shaders/cube.fs",
        "../shaders/text.fs",
        "../shaders/light.fs",
        ""
    };

static uint32
LoadTexture(char *Path)
{
    int32 Width, Height, NChannels;
    uint8 *Data = stbi_load(Path, &Width, &Height, &NChannels, 0);
    GLenum Format;
    uint32 TextureID;
    glGenTextures(1, &TextureID);

    if (Data)
    {

        if (NChannels == 1)
        {
            Format = GL_RED;
        }
        else if (NChannels == 3)
        {
            Format = GL_RGB;
        }
        else if (NChannels == 4)
        {
            Format = GL_RGBA;
        }

        glBindTexture(GL_TEXTURE_2D, TextureID);
        glTexImage2D(GL_TEXTURE_2D, 0, Format, Width, Height, 0, Format, GL_UNSIGNED_BYTE, Data);
        glGenerateMipmap(GL_TEXTURE_2D);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    }
    else
    {
        printf("ERROR::TEXTURE::Failed to load texure! Path: %s\n", Path);
    }
    stbi_image_free(Data);

    return TextureID;
}

static void
LoadTextures(sdl_render* Render)
{
    char** TexturePaths = texture_paths;
    uint32 TexIndex = 0;
    while(strcmp(*TexturePaths, ""))
    {
        Render->Textures[TexIndex++] = LoadTexture(*TexturePaths++);
    }
}

static void
CheckShaderCompilation(uint32 Shader, shader_type Type)
{
    // NOTE(Jovan): Check if shader compilation failed
    int32 success;
    GLchar infoLog[512];
    glGetShaderiv(Shader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(Shader, 512, NULL, infoLog);
        if (Type == VERTEX)
        {
            printf("ERROR::SHADER::VERTEX:COMPILATION_FAILED\n");
        }
        else
        {
            printf("ERROR::SHADER::FRAGMENT:COMPILATION_FAILED\n");
        }
        printf("%s\n", infoLog);
    }
}

static void
CheckShaderLink(uint32 Program)
{
    int32 success;
    GLchar infoLog[512];
    // NOTE(Jovan): Check for program linking failure
    glGetProgramiv(Program, GL_LINK_STATUS, &success);
    if (!success)
    {
        glGetProgramInfoLog(Program, 512, NULL, infoLog);
        printf("ERROR::SHADER_PROGRAM::LINKING_FAILED\n%s\n", infoLog);
    }
}

void LoadShader(const char *file_name, char *shader_str, int max_len)
{
    FILE *file = fopen(file_name, "r");
    if (!file)
    {
        printf("ERROR: opening file for reading: %s\n", file_name);
        return;
    }
    size_t cnt = fread(shader_str, 1, max_len - 1, file);
    if ((int)cnt >= max_len - 1)
    {
        printf("WARNING: file %s too big - truncated.\n", file_name);
    }
    if (ferror(file))
    {
        printf("ERROR: reading shader file %s\n", file_name);
        fclose(file);
        return;
    }
    shader_str[cnt] = 0;
    fclose(file);
}

static void
LoadShaders(sdl_render* Render)
{
    char** VSPaths = vs_paths;
    char** FSPaths = fs_paths;

    char VSSources[SHADER_COUNT][256 * 1024];
    char FSSources[SHADER_COUNT][256 * 1024];

    // TODO(Jovan): Rename shader loader

    uint32 ShaderIndex = 0;
    while(strcmp(*VSPaths, ""))
    {
        LoadShader(*VSPaths++, VSSources[ShaderIndex++], 256 * 1024);
    }
    ShaderIndex = 0;
    while(strcmp(*FSPaths, ""))
    {
        LoadShader(*FSPaths++, FSSources[ShaderIndex++], 256 * 1024);
    }

    for(ShaderIndex = 0;
        ShaderIndex < SHADER_COUNT;
        ++ShaderIndex)
    {
        uint32 VS = glCreateShader(GL_VERTEX_SHADER);
        uint32 FS = glCreateShader(GL_FRAGMENT_SHADER);
        const GLchar *p;
        p = (const GLchar *)VSSources[ShaderIndex];
        glShaderSource(VS, 1, &p, 0);
        p = (const GLchar *)FSSources[ShaderIndex];
        glShaderSource(FS, 1, &p, 0);
        glCompileShader(VS);
        CheckShaderCompilation(VS, VERTEX);
        glCompileShader(FS);
        CheckShaderCompilation(FS, FRAGMENT);
        uint32 ShaderProgram = glCreateProgram();
        glAttachShader(ShaderProgram, VS);
        glAttachShader(ShaderProgram, FS);
        glLinkProgram(ShaderProgram);
        CheckShaderLink(ShaderProgram);
        glDeleteShader(VS);
        glDeleteShader(FS);
        Render->Shaders[ShaderIndex] =ShaderProgram;

    }
}

static void
SDLProcessSimKeyboardButton(sdl_button_state *NewState, bool32 IsDown)
{
    Assert(NewState->EndedDown != IsDown);
    NewState->EndedDown = IsDown;
    ++NewState->HalfTransitionCount;
}

void SDLWindowResize(int32 Width, int32 Height)
{
    glViewport(0, 0, Width, Height);
}

static bool32
SDLHandleEvent(SDL_Event *Event, sdl_input *Input, bool32 *InFocus)
{
    bool32 ShouldQuit = 0;
    switch (Event->type)
    {
        // NOTE(Jovan): Quit event
    case SDL_QUIT:
    {
        ShouldQuit = 1;
    }
    break;

        // NOTE(Jovan): Window events
    case SDL_WINDOWEVENT:
    {
        switch (Event->window.event)
        {
        case SDL_WINDOWEVENT_RESIZED:
        {
            SDLWindowResize(Event->window.data1, Event->window.data2);
        }
        break;

        case SDL_WINDOWEVENT_FOCUS_LOST:
        {
            *InFocus = 0;
            printf("Lost focus, %d\n", *InFocus);
            SDL_SetRelativeMouseMode(SDL_FALSE);
        }
        break;

        case SDL_WINDOWEVENT_TAKE_FOCUS:
        {
            *InFocus = 1;
            printf("Gained focus, %d\n", *InFocus);
            SDL_SetRelativeMouseMode(SDL_TRUE);
        }
        break;
        }
    }
    break;
    case SDL_MOUSEMOTION:
    {
        if (*InFocus == 1)
        {
            Input->MouseController.X = Event->motion.x;
            Input->MouseController.Y = Event->motion.y;
            Input->MouseController.XRel = Event->motion.xrel;
            Input->MouseController.YRel = Event->motion.yrel;
        }
    }
    break;

        // NOTE(Jovan): Keyboard events
    case SDL_KEYDOWN:
    case SDL_KEYUP:
    {
        SDL_Keycode KeyCode = Event->key.keysym.sym;
        bool32 IsDown = (Event->key.state == SDL_PRESSED);

        if (Event->key.state == SDL_RELEASED)
        {
        }
        else if (Event->key.repeat)
        {
        }

        if (!(Event->key.repeat))
        {
            if (KeyCode == SDLK_w)
            {
                SDLProcessSimKeyboardButton(&Input->KeyboardController.MoveForward, IsDown);
            }
            else if (KeyCode == SDLK_a)
            {
                SDLProcessSimKeyboardButton(&Input->KeyboardController.MoveLeft, IsDown);
            }
            else if (KeyCode == SDLK_d)
            {
                SDLProcessSimKeyboardButton(&Input->KeyboardController.MoveRight, IsDown);
            }
            else if (KeyCode == SDLK_s)
            {
                SDLProcessSimKeyboardButton(&Input->KeyboardController.MoveBack, IsDown);
            }
            else if (KeyCode == SDLK_SPACE)
            {
                SDLProcessSimKeyboardButton(&Input->KeyboardController.ShootAction, IsDown);
            }
            else if (KeyCode == SDLK_UP)
            {
                SDLProcessSimKeyboardButton(&Input->KeyboardController.DebugUp, IsDown);
            }
            else if (KeyCode == SDLK_DOWN)
            {
                SDLProcessSimKeyboardButton(&Input->KeyboardController.DebugDown, IsDown);
            }
            else if (KeyCode == SDLK_LEFT)
            {
                SDLProcessSimKeyboardButton(&Input->KeyboardController.DebugLeft, IsDown);
            }
            else if (KeyCode == SDLK_RIGHT)
            {
                SDLProcessSimKeyboardButton(&Input->KeyboardController.DebugRight, IsDown);
            }
            else if (KeyCode == SDLK_f)
            {
                SDLProcessSimKeyboardButton(&Input->KeyboardController.DebugForward, IsDown);
            }
            else if (KeyCode == SDLK_b)
            {
                SDLProcessSimKeyboardButton(&Input->KeyboardController.DebugBack, IsDown);
            }
            else if (KeyCode == SDLK_r)
            {
                SDLProcessSimKeyboardButton(&Input->KeyboardController.DebugReset, IsDown);
            }
            else if (KeyCode == SDLK_p)
            {
                SDLProcessSimKeyboardButton(&Input->KeyboardController.DebugContinue, IsDown);
            }
            else if (KeyCode == SDLK_ESCAPE)
            {
                ShouldQuit = 1;
            }
        }
    }
    break;
    }

    return ShouldQuit;
}

inline static int64
SDLGetWallClock()
{
    int64 Result = SDL_GetPerformanceCounter();
    return Result;
}

inline static real32
SDLGetSecondsElapsed(int64 Start, int64 End)
{
    real32 Result = (real32)(End - Start) / (real32)GlobalPerfCountFrequency;
    return Result;
}

static inline time_t
SDLGetLastWriteTime(char *Filename)
{
    struct stat FileInfo = {};
    time_t Result = 0;
    if (stat(Filename, &FileInfo) != -1)
    {
        Result = FileInfo.st_mtim.tv_sec;
    }

    return Result;
}

static sdl_sim_code
SDLLoadSimCode(char *DynLibName)
{
    sdl_sim_code Result = {};
    Result.SimCodeDynLib = dlopen(DynLibName, RTLD_NOW | RTLD_GLOBAL);
    Result.DynLibLastWriteTime = SDLGetLastWriteTime(DynLibName);

    if (Result.SimCodeDynLib)
    {
        Result.UpdateAndRender = (sim_update_and_render *)
            dlsym(Result.SimCodeDynLib, "SimUpdateAndRender");

        Result.IsValid = (Result.UpdateAndRender != 0);
    }

    if (!Result.IsValid)
    {
        printf("WARNING::SimCode::Not loaded properly, using stub.\n");
        Result.UpdateAndRender = SimUpdateAndRenderStub;
    }

    return Result;
}

static void
SDLUnloadSimCode(sdl_sim_code *SimCode)
{
    if (dlclose(SimCode->SimCodeDynLib))
    {
        printf("ERROR::SimCode::Not unloaded properly.\n");
        char *Error = dlerror();
        printf("%s\n", Error);
    }
    SimCode->SimCodeDynLib = 0;
    SimCode->IsValid = 0;
    SimCode->UpdateAndRender = SimUpdateAndRenderStub;
}

// TODO(Jovan): Optimize and use EBO?
static void
LoadModel(const char *ModelFilepath, uint32 &ModelVAO, uint32 &ModelVBO)
{
    // TODO(Jovan): Implement
}

// TODO(jovan): Refactor!!!
static void
RenderText(sdl_render *Render, std::string Text, real32 X, real32 Y, real32 Scale, glm::vec3 Color)
{
    glUseProgram(Render->Shaders[1]);
    SetUniformF3(Render->Shaders[1], "textColor", Color.x, Color.y, Color.z);
    glActiveTexture(GL_TEXTURE0);
    glBindVertexArray(Render->VAOs[TEXTVAO]);
    std::string::const_iterator c;
    for (c = Text.begin(); c != Text.end(); ++c)
    {
        character ch = Characters[*c];
        real32 xpos = X + ch.Bearing.x * Scale;
        real32 ypos = Y - (ch.Size.y - ch.Bearing.y) * Scale;

        real32 w = ch.Size.x * Scale;
        real32 h = ch.Size.y * Scale;

        real32 vertices[6][4] = {
            {xpos, ypos + h, 0.0f, 0.0f},
            {xpos, ypos, 0.0f, 1.0f},
            {xpos + w, ypos, 1.0f, 1.0f},

            {xpos, ypos + h, 0.0f, 0.0f},
            {xpos + w, ypos, 1.0f, 1.0f},
            {xpos + w, ypos + h, 1.0f, 0.0f}};

        glBindTexture(GL_TEXTURE_2D, ch.TextureId);
        glBindBuffer(GL_ARRAY_BUFFER, Render->VBOs[TEXTVBO]);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        X += (ch.Advance >> 6) * Scale;
    }
    glBindVertexArray(0);
    glBindTexture(GL_TEXTURE_2D, 0);
}

int main()
{

    // NOTE(Jovan): Memory allocation
    memory SimMemory = {};
    SimMemory.PermanentStorageSize = Mebibytes(64);
    SimMemory.TransientStorageSize = Mebibytes(256);
    uint64 TotalStorageSize = SimMemory.PermanentStorageSize +
                              SimMemory.TransientStorageSize;
    SimMemory.PermanentStorage = mmap(0,
                                      TotalStorageSize,
                                      PROT_READ | PROT_WRITE,
                                      MAP_ANON | MAP_PRIVATE,
                                      -1,
                                      0);
    SimMemory.TransientStorage = ((uint8 *)SimMemory.PermanentStorage +
                                  SimMemory.PermanentStorageSize);
    Assert(SimMemory.PermanentStorage != (void *)-1);
    Assert(SimMemory.TransientStorage != (void *)-1);
    // NOTE(Jovan): SDL stuff
    // --------------------
    SDL_Window *Window = 0;
    SDL_GLContext GLContext = 0;
    GlobalPerfCountFrequency = SDL_GetPerformanceFrequency();
    real32 SimUpdateHz = 60.0f;
    real32 TargetSecondsPerFrame = 1.0f / SimUpdateHz;

    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 0);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

    if (SDL_Init(SDL_INIT_VIDEO) > 0)
    {
        printf("ERROR::SDL::Failed to init!\n");
        return 1;
    }

    uint32 Flags = SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE;

    Window = SDL_CreateWindow("NANS Projekat",
                              SDL_WINDOWPOS_CENTERED,
                              SDL_WINDOWPOS_CENTERED,
                              DEFAULT_WINDOW_WIDTH,
                              DEFAULT_WINDOW_HEIGHT,
                              Flags);
    GLContext = SDL_GL_CreateContext(Window);

    if (!GLContext)
    {
        printf("ERROR::SDL::Failed to create GL context!\n");
        return 1;
    }

    // NOTE(Jovan): Get window size and set viewport
    int32 Width, Height;
    SDL_GetWindowSize(Window, &Width, &Height);
    glewInit();
    glViewport(0, 0, Width, Height);

    // NOTE(Jovan): End of SDL stuff
    // -----------------------------

    // NOTE(Jovan): Font stuff
    // -----------------------

    FT_Library ft;
    if (FT_Init_FreeType(&ft))
    {
        printf("ERROR::FREETYPE: Failed to init freetype");
        return 1;
    }

    FT_Face face;
    if (FT_New_Face(ft, "../res/fonts/poppins.ttf", 0, &face))
    {
        printf("ERROR:FREETYPE: Failed to load Poppins font");
        return 1;
    }

    FT_Set_Pixel_Sizes(face, 0, 24);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    for (uint8 C = 0;
         C < 128;
         ++C)
    {
        if (FT_Load_Char(face, C, FT_LOAD_RENDER))
        {
            printf("ERROR::FREETYPE: Failed to load glyph %c\n", C);
            continue;
        }
        uint32 Texture;
        glGenTextures(1, &Texture);
        glBindTexture(GL_TEXTURE_2D, Texture);
        glTexImage2D(
            GL_TEXTURE_2D,
            0,
            GL_RED,
            face->glyph->bitmap.width,
            face->glyph->bitmap.rows,
            0,
            GL_RED,
            GL_UNSIGNED_BYTE,
            face->glyph->bitmap.buffer);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        character Character = {
            Texture,
            glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
            glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
            face->glyph->advance.x};
        Characters.insert(std::pair<uint8, character>(C, Character));
    }
    FT_Done_Face(face);
    FT_Done_FreeType(ft);

    // NOTE(Jovan): End of font stuff

    // NOTE(Jovan): GL modeling and buffering
    // --------------------------------------

    real32 CubeVertices[] =
        {
            // X  |  Y   |  Z  | Tex coords  | Normals
            -0.5f, -0.5f, -0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
            0.5f, -0.5f, -0.5f, 1.0f, 0.0f, 0.0f, 0.0f, -1.0f,
            0.5f, 0.5f, -0.5f, 1.0f, 1.0f, 0.0f, 0.0f, -1.0f,
            0.5f, 0.5f, -0.5f, 1.0f, 1.0f, 0.0f, 0.0f, -1.0f,
            -0.5f, 0.5f, -0.5f, 0.0f, 1.0f, 0.0f, 0.0f, -1.0f,
            -0.5f, -0.5f, -0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,

            -0.5f, -0.5f, 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
            0.5f, -0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f,
            0.5f, 0.5f, 0.5f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f,
            0.5f, 0.5f, 0.5f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f,
            -0.5f, 0.5f, 0.5f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f,
            -0.5f, -0.5f, 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,

            -0.5f, 0.5f, 0.5f, 1.0f, 0.0f, -1.0f, 0.0f, 0.0f,
            -0.5f, 0.5f, -0.5f, 1.0f, 1.0f, -1.0f, 0.0f, 0.0f,
            -0.5f, -0.5f, -0.5f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f,
            -0.5f, -0.5f, -0.5f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f,
            -0.5f, -0.5f, 0.5f, 0.0f, 0.0f, -1.0f, 0.0f, 0.0f,
            -0.5f, 0.5f, 0.5f, 1.0f, 0.0f, -1.0f, 0.0f, 0.0f,

            0.5f, 0.5f, 0.5f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f,
            0.5f, 0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,
            0.5f, -0.5f, -0.5f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f,
            0.5f, -0.5f, -0.5f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f,
            0.5f, -0.5f, 0.5f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
            0.5f, 0.5f, 0.5f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f,

            -0.5f, -0.5f, -0.5f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f,
            0.5f, -0.5f, -0.5f, 1.0f, 1.0f, 0.0f, -1.0f, 0.0f,
            0.5f, -0.5f, 0.5f, 1.0f, 0.0f, 0.0f, -1.0f, 0.0f,
            0.5f, -0.5f, 0.5f, 1.0f, 0.0f, 0.0f, -1.0f, 0.0f,
            -0.5f, -0.5f, 0.5f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
            -0.5f, -0.5f, -0.5f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f,

            -0.5f, 0.5f, -0.5f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f,
            0.5f, 0.5f, -0.5f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f,
            0.5f, 0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f,
            0.5f, 0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f,
            -0.5f, 0.5f, 0.5f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f,
            -0.5f, 0.5f, -0.5f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f};

    real32 LightVertices[] =
        {
            -0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f, 0.5f, -0.5f,
            0.5f, 0.5f, -0.5f,
            -0.5f, 0.5f, -0.5f,
            -0.5f, -0.5f, -0.5f,

            -0.5f, -0.5f, 0.5f,
            0.5f, -0.5f, 0.5f,
            0.5f, 0.5f, 0.5f,
            0.5f, 0.5f, 0.5f,
            -0.5f, 0.5f, 0.5f,
            -0.5f, -0.5f, 0.5f,

            -0.5f, 0.5f, 0.5f,
            -0.5f, 0.5f, -0.5f,
            -0.5f, -0.5f, -0.5f,
            -0.5f, -0.5f, -0.5f,
            -0.5f, -0.5f, 0.5f,
            -0.5f, 0.5f, 0.5f,

            0.5f, 0.5f, 0.5f,
            0.5f, 0.5f, -0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, 0.5f,
            0.5f, 0.5f, 0.5f,

            -0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, 0.5f,
            0.5f, -0.5f, 0.5f,
            -0.5f, -0.5f, 0.5f,
            -0.5f, -0.5f, -0.5f,

            -0.5f, 0.5f, -0.5f,
            0.5f, 0.5f, -0.5f,
            0.5f, 0.5f, 0.5f,
            0.5f, 0.5f, 0.5f,
            -0.5f, 0.5f, 0.5f,
            -0.5f, 0.5f, -0.5f};

    uint32 RectIndices[] =
        {
            0, 1, 3,
            0, 2, 3};

    real32 FloorVertices[] =
        {
            // X |   Y   |  Z  |  TexX | TexY
            50.0f, -0.5f, 50.0f, 2.0f, 0.0f,
            -50.0f, -0.5f, 50.0f, 0.0f, 0.0f,
            -50.0f, -0.5f, -50.0f, 0.0f, 2.0f,

            50.0f, -0.5f, 50.0f, 2.0f, 0.0f,
            -50.0f, -0.5f, -50.0f, 0.0f, 2.0f,
            50.0f, -0.5f, -50.0f, 2.0f, 2.0f};

    // TODO(Jovan): Use vertices for one line only so different lines other
    // than the coordinate system can be created / drawn
    real32 CoordinateVertices[] =
        {
            // X | Y   | Z   | R  | G  | B
            0.0f, 0.0f, 0.0f,
            1.0f, 0.0f, 0.0f};

    uint32 Stacks = 20;
    uint32 Slices = 20;

    //  real32 SphereVertices[(Stacks + 1) * (Slices + 1) * 3];
    // TODO(Jovan): Sphere texture
    real32 SphereVertices[(Stacks + 1) * (Slices + 1) * 5];
    uint32 SphereIndices[(Slices * Stacks + Slices) * 6];

    int32 Index = 0;
    for (uint32 i = 0; i <= Stacks; i++)
    {
        real32 V = (real32)i / (real32)Stacks;
        real32 Phi = V * Pi32;

        for (uint32 j = 0; j <= Slices; j++)
        {
            real32 U = (real32)j / (real32)Slices;
            real32 Theta = U * (Pi32 * 2);

            real32 x = cos(Theta) * sin(Phi);
            real32 y = cos(Phi);
            real32 z = sin(Theta) * sin(Phi);

            SphereVertices[Index++] = x;
            SphereVertices[Index++] = y;
            SphereVertices[Index++] = z;
            // // TODO(Jovan): Texture sphere
            SphereVertices[Index++] = (real32)j / (real32)Slices;
            SphereVertices[Index++] = (real32)i / (real32)Stacks;
        }
    }
    Index = 0;
    for (uint32 i = 0; i < Slices * Stacks + Slices; i++)
    {
        SphereIndices[Index++] = i;
        SphereIndices[Index++] = i + Slices + 1;
        SphereIndices[Index++] = i + Slices;

        SphereIndices[Index++] = i + Slices + 1;
        SphereIndices[Index++] = i;
        SphereIndices[Index++] = i + 1;
    }

    // NOTE(Jovan): Creating shaders and shader programs
    // TODO(Jovan): To transient storage
    sdl_render Render = {};
    LoadShaders(&Render);

    // NOTE(Jovan): VAO, EBO, VBO
    // TODO(Jovan): Gen arrays inside Render directly

    // NOTE(Jovan): Floor data
    glGenVertexArrays(VAO_COUNT, Render.VAOs);
    glGenBuffers(VBO_COUNT, Render.VBOs);

    glBindVertexArray(Render.VAOs[FLOORVAO]);
    // TODO(Jovan): Change to static draw
    glBindBuffer(GL_ARRAY_BUFFER, Render.VBOs[FLOORVBO]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(FloorVertices), FloorVertices,
                 GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(real32),
                          (void *)(0));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(real32),
                          (void *)(3 * sizeof(real32)));
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    // NOTE(Jovan): Cube data
    glBindVertexArray(Render.VAOs[CUBEVAO]);

    // TODO(Jovan): Change to static draw
    glBindBuffer(GL_ARRAY_BUFFER, Render.VBOs[CUBEVBO]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(CubeVertices), CubeVertices,
                 GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(real32),
                          (void *)0);
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(real32),
                          (void *)(3 * sizeof(real32)));
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    // NOTE(Jovan): Sphere data
    glBindVertexArray(Render.VAOs[SPHEREVAO]);

    glBindBuffer(GL_ARRAY_BUFFER, Render.VBOs[SPHEREVBO]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(SphereVertices), SphereVertices, GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(real32),
                          (void *)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(real32),
                          (void *)(3 * sizeof(real32)));
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    // NOTE(Jovan): Text data
    glBindVertexArray(Render.VAOs[TEXTVAO]);
    glBindBuffer(GL_ARRAY_BUFFER, Render.VBOs[TEXTVBO]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(real32) * 6 * 4, NULL, GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(real32), 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    // NOTE(Jovan): Light data
    glBindVertexArray(Render.VAOs[LIGHTVAO]);
    glBindBuffer(GL_ARRAY_BUFFER, Render.VBOs[CUBEVBO]);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(real32), (void *)0);
    glEnableVertexAttribArray(0);

    // NOTE(Jovan): Textures
    LoadTextures(&Render);

    // NOTE(Jovan): End of GL modeling and buffering
    // ---------------------------------------------
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // NOTE(Jovan): Main loop
    bool32 Running = 1;
    uint64 LastCounter = SDLGetWallClock();
    real32 dt = 0.0f;
    sdl_input SimInput[2];
    sdl_input *NewInput = &SimInput[0];
    sdl_input *OldInput = &SimInput[1];
    *NewInput = {};
    *OldInput = {};

    // NOTE(Jovan): Grab mouse
    SDL_SetRelativeMouseMode(SDL_TRUE);
    bool32 InFocus = 1;

    sdl_sim_code Sim = SDLLoadSimCode("nans.so");

    Render.Indices = SphereIndices;
    Render.Num = ArrayCount(SphereIndices);

    // NOTE(Jovan): Assimp model loading
    // TODO(Jovan): Move to proper location
    // -------------------------------------
    const char *ModelFilepath = "../res/model/amongus.obj";
    const struct aiScene *Scene = aiImportFile(ModelFilepath, aiProcess_Triangulate | aiProcess_FlipUVs);
    if (!Scene || Scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !Scene->mRootNode)
    {
        printf("ERROR::ASSIMP MODEL LOADING FAILED\n%s", aiGetErrorString());
    }
    printf("Loaded model with %d meshes\n", Scene->mNumMeshes);
    const struct aiMesh *Mesh = Scene->mMeshes[0];
    int32 NumMeshVertices = Mesh->mNumVertices;
    printf("Mesh with %d vertices\n", NumMeshVertices);
    real32 MeshVertices[5 * NumMeshVertices];
    uint32 VertexCount = 0;
    for (int VertexIndex = 0;
         VertexIndex < NumMeshVertices;
         ++VertexIndex)
    {
        MeshVertices[VertexCount++] = Mesh->mVertices[VertexIndex].x;
        MeshVertices[VertexCount++] = Mesh->mVertices[VertexIndex].y;
        MeshVertices[VertexCount++] = Mesh->mVertices[VertexIndex].z;
        MeshVertices[VertexCount++] = Mesh->mTextureCoords[0][VertexIndex].x;
        MeshVertices[VertexCount++] = Mesh->mTextureCoords[0][VertexIndex].y;
    }
    printf("Mesh with %d faces\n", Mesh->mNumFaces);
    int32 NumMeshFaces = Mesh->mNumFaces;

    // NOTE(Jovan): x3 because faces are triangles due to aiProcess_Triangulate flag
    uint32 MeshFaceIndices[NumMeshFaces * 3];
    uint32 IndexCount = 0;
    for (uint32 FaceIndex = 0;
         FaceIndex < Mesh->mNumFaces;
         ++FaceIndex)
    {
        MeshFaceIndices[IndexCount++] = Mesh->mFaces[FaceIndex].mIndices[0];
        MeshFaceIndices[IndexCount++] = Mesh->mFaces[FaceIndex].mIndices[1];
        MeshFaceIndices[IndexCount++] = Mesh->mFaces[FaceIndex].mIndices[2];
    }

    glBindVertexArray(Render.VAOs[MODELVAO]);

    glBindBuffer(GL_ARRAY_BUFFER, Render.VBOs[MODELVBO]);
    glBufferData(GL_ARRAY_BUFFER, 5 * NumMeshVertices * sizeof(real32), MeshVertices, GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(real32),
                          (void *)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(real32),
                          (void *)(3 * sizeof(real32)));

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    aiReleaseImport(Scene);
    //---

    Render.ModelIndices = MeshFaceIndices;
    Render.ModelNum = ArrayCount(MeshFaceIndices);
    // NOTE(Jovan): End assimp model loading

    while (Running)
    {

        time_t NewDynLibWriteTime = SDLGetLastWriteTime("nans.so");
        if ((NewDynLibWriteTime != Sim.DynLibLastWriteTime))
        {
            printf("nans.so difference: %ld\n", NewDynLibWriteTime - Sim.DynLibLastWriteTime);
            printf("Code changed!\n");
            SDLUnloadSimCode(&Sim);
            SDL_Delay(100);
            Sim = SDLLoadSimCode("nans.so");
        }

        sdl_keyboard *OldKeyboardController = &OldInput->KeyboardController;
        sdl_keyboard *NewKeyboardController = &NewInput->KeyboardController;
        *NewKeyboardController = {};

        sdl_mouse *OldMouseController = &OldInput->MouseController;
        sdl_mouse *NewMouseController = &NewInput->MouseController;
        *NewMouseController = {};
        NewMouseController->Sensitivity = 0.5f;

        for (uint32 ButtonIndex = 0;
             ButtonIndex < ArrayCount(NewKeyboardController->Buttons);
             ++ButtonIndex)
        {
            NewKeyboardController->Buttons[ButtonIndex].EndedDown =
                OldKeyboardController->Buttons[ButtonIndex].EndedDown;
        }

        NewInput->MouseController.X = OldInput->MouseController.X;
        NewInput->MouseController.Y = OldInput->MouseController.Y;
        NewInput->MouseController.XRel = OldInput->MouseController.XRel;
        NewInput->MouseController.YRel = OldInput->MouseController.YRel;

        SDL_Event Event;
        while (SDL_PollEvent(&Event))
        {
            // NOTE(Jovan): Check for exit
            if (SDLHandleEvent(&Event, NewInput, &InFocus))
            {
                Running = 0;
            }
        }

        // NOTE(Jovan): Work timing
        int64 WorkCounter = SDLGetWallClock();
        real32 WorkSecondsElapsed = SDLGetSecondsElapsed(LastCounter, WorkCounter);
        real32 SecondsElapsedForFrame = WorkSecondsElapsed;
        if (SecondsElapsedForFrame < TargetSecondsPerFrame)
        {
            while (SecondsElapsedForFrame < TargetSecondsPerFrame)
            {
                // TODO(Jovan): Check sleep granulaity on linux
                uint32 MSToSleepFor = (uint32)(1000.0f * (TargetSecondsPerFrame -
                                                          SecondsElapsedForFrame));
                SDL_Delay(MSToSleepFor);
                SecondsElapsedForFrame = SDLGetSecondsElapsed(LastCounter, SDLGetWallClock());
            }
        }
        else
        {
            // TODO(Jovan): Missed frame rate???
        }

        glClearColor(0.8f, 0.8f, 0.8f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        Sim.UpdateAndRender(&SimMemory, NewInput, &Render, dt);
        // NOTE(Jovan): Text rendering
        // ---------------------------
        glm::mat4 Projection = glm::ortho(0.0f, (real32)Width, 0.0f, (real32)Height);
        RenderText(&Render, "Jovan Ivosevic RA30/2017", 25.0f, 25.0f, 1.0f, glm::vec3(0.91f, 0.30f, 0.24f));
        SetUniformM4(Render.Shaders[TEXTSH], "projection", Projection);
        // NOTE(Jovan): End text rendering

        // NOTE(Jovan): Timing
        int64 EndCounter = SDLGetWallClock();
        real32 SPerFrame = SDLGetSecondsElapsed(LastCounter, EndCounter);
        real32 FPS = 1.0f / SPerFrame;
        real32 MSPerFrame = 1000.0f * SDLGetSecondsElapsed(LastCounter, EndCounter);
        dt = SPerFrame;
        LastCounter = EndCounter;
        char Buffer[256];
        sprintf(Buffer, "MSPerFrame = %f, FPS = %f", MSPerFrame, FPS);
        // TODO(Jovan): Move? Use DEBUG text rendering?
        RenderText(&Render, Buffer, ((real32)Width) / 2.0f, (real32)Height - 25.0f, 1.0f, glm::vec3(0.0f, 1.0f, 0.0f));
        // NOTE(Jovan): Swap buffers // TODO(Jovan): Move???
        SDL_GL_SwapWindow(Window);
        SDL_SetWindowTitle(Window, Buffer);

        sdl_input *Temp = NewInput;
        NewInput = OldInput;
        OldInput = Temp;
    }
    return 0;
}
