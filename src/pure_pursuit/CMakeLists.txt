cmake_minimum_required(VERSION 3.0.2)
project(pure_pursuit)

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

add_executable(pure_pursuit_node src/pure_pursuit.cpp)

target_link_libraries(pure_pursuit_node
  ${catkin_LIBRARIES}
  ${PYTHON_LIBRARIES}
)

