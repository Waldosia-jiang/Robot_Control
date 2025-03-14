cmake_minimum_required(VERSION 3.0.2)
project(stanley_controller)

# 设置 C++17 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  std_msgs
  sensor_msgs
  geometry_msgs
)

catkin_package()

find_package(Eigen3 REQUIRED)
find_package(PythonLibs REQUIRED)

include_directories(
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include
)

add_executable(stanley_controller_node src/stanley_controller.cpp)

target_link_libraries(stanley_controller_node
  ${catkin_LIBRARIES}
  ${PYTHON_LIBRARIES}
)

