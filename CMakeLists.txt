cmake_minimum_required(VERSION 3.16)
project(starfish)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)

if (APPLE)
   include_directories(/usr/local/include)
   link_directories(/usr/local/lib)
endif()

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/hatchetfish/include)

add_library(${PROJECT_NAME} STATIC
    src/ArHosekSkyModel.cpp
    src/starfish-astronomy.cpp
    src/starfish-hosek-wilkie.cpp
    src/starfish-pbsky.cpp
    )

if(${CMAKE_VERSION} VERSION_GREATER "3.16.0")
    message("Starfish using precompiled headers")
    target_precompile_headers(${PROJECT_NAME} PRIVATE src/starfish_pch.cpp)
endif()

find_package(GLEW REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE GLEW::GLEW)
