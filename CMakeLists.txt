cmake_minimum_required(VERSION 3.10)
project(MeltSim VERSION 0.1.0)

# Set policies
cmake_policy(SET CMP0002 NEW)
cmake_policy(SET CMP0077 NEW)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(APPLE)
    add_definitions(-DGL_SILENCE_DEPRECATION)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error=deprecated-declarations")
endif()

include(FetchContent)

# GLFW
FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG 3.3.8
)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
if(APPLE)
    set(GLFW_USE_CHDIR OFF CACHE BOOL "" FORCE)
    set(GLFW_USE_MENUBAR OFF CACHE BOOL "" FORCE)
endif()

# GLM
FetchContent_Declare(
    glm
    GIT_REPOSITORY https://github.com/g-truc/glm.git
    GIT_TAG 0.9.9.8
)
set(GLM_TEST_ENABLE OFF CACHE BOOL "" FORCE)
set(GLM_BUILD_TESTS OFF CACHE BOOL "" FORCE)

# Glad
set(GLAD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/glad")
add_library(glad "${GLAD_DIR}/src/glad.c")
target_include_directories(glad PUBLIC "${GLAD_DIR}/include")
if(APPLE)
    target_compile_definitions(glad PRIVATE GLAD_GL_IMPLEMENTATION)
endif()

# Assimp
find_package(assimp REQUIRED)
include_directories(/opt/homebrew/include)
link_directories(/opt/homebrew/lib)

# Add source files
add_executable(MeltSim
    src/main.cpp
    src/mesh.cpp
    src/mesh.h
    src/floor.cpp
    src/floor.h
    src/plane.cpp
    src/plane.h
    src/input.cpp
    src/input.h
    src/shader_utils.cpp
    src/shader_utils.h
    src/leaven/surfaceSampler.cpp
    # src/particle.cpp
    src/particle.h
    src/particleRenderer.cpp
    src/particleRenderer.h
    src/mpm.cpp
    src/mpm.h
    src/grid.cpp
    src/grid.h
    src/particleSplatter.cpp
    src/particleSplatter.h
    src/gridRenderer.cpp
    src/gridRenderer.h
    src/gridRenderer.cpp
    src/gridRenderer.h
)

# Now you can add include directories
target_include_directories(MeltSim PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen-3.4.0
)

# Executable
# add_executable(${PROJECT_NAME} ${SOURCES})

# Includes
target_include_directories(${PROJECT_NAME} PRIVATE
    src
    "${GLAD_DIR}/include"
)

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE
    glfw
    glad
    glm
    assimp
)

# macOS frameworks
if(APPLE)
    target_link_libraries(${PROJECT_NAME} PRIVATE
        "-framework OpenGL"
        "-framework Cocoa"
        "-framework IOKit"
        "-framework CoreFoundation"
    )
endif()

# Copy assets
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/assets")
    file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/assets" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
endif()

# Enable dependencies
FetchContent_MakeAvailable(glfw glm)