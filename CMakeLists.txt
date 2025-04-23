cmake_minimum_required(VERSION 3.5)
project(MeltSim VERSION 0.1.0)

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
    src/input.cpp
    src/input.h
    src/shader_utils.cpp
    src/shader_utils.h
    src/leaven/surfaceSampler.cpp
    src/particleRenderer.cpp
    src/particleRenderer.h
    src/particleSampler.cpp
    src/particleSampler.h
    # src/mpmSimulation.cpp
    src/mpmSimulator.h
    src/mpmSimulator.cpp
    
    # src/mpmSimulation.h
    src/taichi.h

    # Add more source files as needed
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