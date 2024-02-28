	cmake_minimum_required(VERSION 3.6)
project(SketchMeshVR)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


# libigl
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"           OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK              "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_LIM               "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)
option(LIBIGL_WITH_PNG               "Use PNG"            OFF)
option(LIBIGL_WITH_PYTHON            "Use Python"         OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_XML               "Use XML"            OFF)


find_package(LIBIGL REQUIRED QUIET)


# Prepare the build environment
include_directories(${LIBIGL_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/../OculusSDK/LibOVR/Include ${PROJECT_SOURCE_DIR}/../OculusSDK/LibOVR/Src ${PROJECT_SOURCE_DIR}/../OVRAvatarSDK/Include ${PROJECT_SOURCE_DIR}/../OVRPlatformSDK/Include ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/../libiglVR/external/triangle)
#add_definitions(${LIBIGL_DEFINITIONS})

# Add your project files
file(GLOB SRCFILES src/*.cpp)
add_executable(${PROJECT_NAME}_bin ${SRCFILES})
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME}_bin)
target_link_libraries(${PROJECT_NAME}_bin ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES} Winmm.lib ws2_32.lib)
target_link_libraries(${PROJECT_NAME}_bin debug ${PROJECT_SOURCE_DIR}/../OculusSDK/LibOVR/Lib/Windows/x64/Debug/VS2015/LibOVR.lib)
target_link_libraries(${PROJECT_NAME}_bin optimized ${PROJECT_SOURCE_DIR}/../OculusSDK/LibOVR/Lib/Windows/x64/Release/VS2015/LibOVR.lib)
target_link_libraries(${PROJECT_NAME}_bin ${PROJECT_SOURCE_DIR}/../OVRAvatarSDK/Windows/libovravatar.lib)
target_link_libraries(${PROJECT_NAME}_bin ${PROJECT_SOURCE_DIR}/../OVRPlatformSDK/Windows/libOVRPlatform64_1.lib)
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl_glfw igl::opengl_glfw_imgui igl::triangle)
