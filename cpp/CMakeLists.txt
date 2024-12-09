cmake_minimum_required(VERSION 3.10)

project(MyProject)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(nlohmann_json_DIR "C:/Users/Timofey/vcpkg/installed/x64-windows/share/nlohmann_json")

find_package(nlohmann_json REQUIRED)
add_executable(MyExecutable OnlineVariance.cpp)

target_link_libraries(MyExecutable PRIVATE nlohmann_json::nlohmann_json)