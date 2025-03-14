# generated from catkin/cmake/template/pkgConfig.cmake.in

# append elements to a list and remove existing duplicates from the list
# copied from catkin/cmake/list_append_deduplicate.cmake to keep pkgConfig
# self contained
macro(_list_append_deduplicate listname)
  if(NOT "${ARGN}" STREQUAL "")
    if(${listname})
      list(REMOVE_ITEM ${listname} ${ARGN})
    endif()
    list(APPEND ${listname} ${ARGN})
  endif()
endmacro()

# append elements to a list if they are not already in the list
# copied from catkin/cmake/list_append_unique.cmake to keep pkgConfig
# self contained
macro(_list_append_unique listname)
  foreach(_item ${ARGN})
    list(FIND ${listname} ${_item} _index)
    if(_index EQUAL -1)
      list(APPEND ${listname} ${_item})
    endif()
  endforeach()
endmacro()

# pack a list of libraries with optional build configuration keywords
# copied from catkin/cmake/catkin_libraries.cmake to keep pkgConfig
# self contained
macro(_pack_libraries_with_build_configuration VAR)
  set(${VAR} "")
  set(_argn ${ARGN})
  list(LENGTH _argn _count)
  set(_index 0)
  while(${_index} LESS ${_count})
    list(GET _argn ${_index} lib)
    if("${lib}" MATCHES "^(debug|optimized|general)$")
      math(EXPR _index "${_index} + 1")
      if(${_index} EQUAL ${_count})
        message(FATAL_ERROR "_pack_libraries_with_build_configuration() the list of libraries '${ARGN}' ends with '${lib}' which is a build configuration keyword and must be followed by a library")
      endif()
      list(GET _argn ${_index} library)
      list(APPEND ${VAR} "${lib}${CATKIN_BUILD_CONFIGURATION_KEYWORD_SEPARATOR}${library}")
    else()
      list(APPEND ${VAR} "${lib}")
    endif()
    math(EXPR _index "${_index} + 1")
  endwhile()
endmacro()

# unpack a list of libraries with optional build configuration keyword prefixes
# copied from catkin/cmake/catkin_libraries.cmake to keep pkgConfig
# self contained
macro(_unpack_libraries_with_build_configuration VAR)
  set(${VAR} "")
  foreach(lib ${ARGN})
    string(REGEX REPLACE "^(debug|optimized|general)${CATKIN_BUILD_CONFIGURATION_KEYWORD_SEPARATOR}(.+)$" "\\1;\\2" lib "${lib}")
    list(APPEND ${VAR} "${lib}")
  endforeach()
endmacro()


if(stanley_controller_CONFIG_INCLUDED)
  return()
endif()
set(stanley_controller_CONFIG_INCLUDED TRUE)

# set variables for source/devel/install prefixes
if("TRUE" STREQUAL "TRUE")
  set(stanley_controller_SOURCE_PREFIX /home/wheeltec/Robot_Control/src/stanley_controller)
  set(stanley_controller_DEVEL_PREFIX /home/wheeltec/Robot_Control/devel)
  set(stanley_controller_INSTALL_PREFIX "")
  set(stanley_controller_PREFIX ${stanley_controller_DEVEL_PREFIX})
else()
  set(stanley_controller_SOURCE_PREFIX "")
  set(stanley_controller_DEVEL_PREFIX "")
  set(stanley_controller_INSTALL_PREFIX /home/wheeltec/Robot_Control/install)
  set(stanley_controller_PREFIX ${stanley_controller_INSTALL_PREFIX})
endif()

# warn when using a deprecated package
if(NOT "" STREQUAL "")
  set(_msg "WARNING: package 'stanley_controller' is deprecated")
  # append custom deprecation text if available
  if(NOT "" STREQUAL "TRUE")
    set(_msg "${_msg} ()")
  endif()
  message("${_msg}")
endif()

# flag project as catkin-based to distinguish if a find_package()-ed project is a catkin project
set(stanley_controller_FOUND_CATKIN_PROJECT TRUE)

if(NOT " " STREQUAL " ")
  set(stanley_controller_INCLUDE_DIRS "")
  set(_include_dirs "")
  if(NOT " " STREQUAL " ")
    set(_report "Check the issue tracker '' and consider creating a ticket if the problem has not been reported yet.")
  elseif(NOT " " STREQUAL " ")
    set(_report "Check the website '' for information and consider reporting the problem.")
  else()
    set(_report "Report the problem to the maintainer 'wheeltec <wheeltec@todo.todo>' and request to fix the problem.")
  endif()
  foreach(idir ${_include_dirs})
    if(IS_ABSOLUTE ${idir} AND IS_DIRECTORY ${idir})
      set(include ${idir})
    elseif("${idir} " STREQUAL "include ")
      get_filename_component(include "${stanley_controller_DIR}/../../../include" ABSOLUTE)
      if(NOT IS_DIRECTORY ${include})
        message(FATAL_ERROR "Project 'stanley_controller' specifies '${idir}' as an include dir, which is not found.  It does not exist in '${include}'.  ${_report}")
      endif()
    else()
      message(FATAL_ERROR "Project 'stanley_controller' specifies '${idir}' as an include dir, which is not found.  It does neither exist as an absolute directory nor in '/home/wheeltec/Robot_Control/src/stanley_controller/${idir}'.  ${_report}")
    endif()
    _list_append_unique(stanley_controller_INCLUDE_DIRS ${include})
  endforeach()
endif()

set(libraries "")
foreach(library ${libraries})
  # keep build configuration keywords, target names and absolute libraries as-is
  if("${library}" MATCHES "^(debug|optimized|general)$")
    list(APPEND stanley_controller_LIBRARIES ${library})
  elseif(${library} MATCHES "^-l")
    list(APPEND stanley_controller_LIBRARIES ${library})
  elseif(${library} MATCHES "^-")
    # This is a linker flag/option (like -pthread)
    # There's no standard variable for these, so create an interface library to hold it
    if(NOT stanley_controller_NUM_DUMMY_TARGETS)
      set(stanley_controller_NUM_DUMMY_TARGETS 0)
    endif()
    # Make sure the target name is unique
    set(interface_target_name "catkin::stanley_controller::wrapped-linker-option${stanley_controller_NUM_DUMMY_TARGETS}")
    while(TARGET "${interface_target_name}")
      math(EXPR stanley_controller_NUM_DUMMY_TARGETS "${stanley_controller_NUM_DUMMY_TARGETS}+1")
      set(interface_target_name "catkin::stanley_controller::wrapped-linker-option${stanley_controller_NUM_DUMMY_TARGETS}")
    endwhile()
    add_library("${interface_target_name}" INTERFACE IMPORTED)
    if("${CMAKE_VERSION}" VERSION_LESS "3.13.0")
      set_property(
        TARGET
        "${interface_target_name}"
        APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES "${library}")
    else()
      target_link_options("${interface_target_name}" INTERFACE "${library}")
    endif()
    list(APPEND stanley_controller_LIBRARIES "${interface_target_name}")
  elseif(TARGET ${library})
    list(APPEND stanley_controller_LIBRARIES ${library})
  elseif(IS_ABSOLUTE ${library})
    list(APPEND stanley_controller_LIBRARIES ${library})
  else()
    set(lib_path "")
    set(lib "${library}-NOTFOUND")
    # since the path where the library is found is returned we have to iterate over the paths manually
    foreach(path /home/wheeltec/Robot_Control/devel/lib;/home/wheeltec/Robot_Control/devel/lib;/home/wheeltec/four_ws_robot/devel/lib;/home/wheeltec/PythonRobotics/devel/lib;/home/wheeltec/catkin_ws/devel/lib;/home/wheeltec/four_ws_navigation/devel/lib;/home/wheeltec/wpr_simulation/devel/lib;/home/wheeltec/demo03_ws/devel/lib;/home/wheeltec/ugv_gazebo_sim/devel/lib;/home/wheeltec/wheeltec_arm/devel/lib;/home/wheeltec/.autoware/install/ymc/lib;/home/wheeltec/.autoware/install/xsens_driver/lib;/home/wheeltec/.autoware/install/wf_simulator/lib;/home/wheeltec/.autoware/install/lattice_planner/lib;/home/wheeltec/.autoware/install/waypoint_planner/lib;/home/wheeltec/.autoware/install/waypoint_maker/lib;/home/wheeltec/.autoware/install/way_planner/lib;/home/wheeltec/.autoware/install/vlg22c_cam/lib;/home/wheeltec/.autoware/install/vision_ssd_detect/lib;/home/wheeltec/.autoware/install/vision_segment_enet_detect/lib;/home/wheeltec/.autoware/install/vel_pose_diff_checker/lib;/home/wheeltec/.autoware/install/vehicle_socket/lib;/home/wheeltec/.autoware/install/vehicle_model/lib;/home/wheeltec/.autoware/install/vehicle_gazebo_simulation_launcher/lib;/home/wheeltec/.autoware/install/vehicle_gazebo_simulation_interface/lib;/home/wheeltec/.autoware/install/vehicle_engage_panel/lib;/home/wheeltec/.autoware/install/vehicle_description/lib;/home/wheeltec/.autoware/install/op_utilities/lib;/home/wheeltec/.autoware/install/op_simulation_package/lib;/home/wheeltec/.autoware/install/op_local_planner/lib;/home/wheeltec/.autoware/install/op_global_planner/lib;/home/wheeltec/.autoware/install/lidar_kf_contour_track/lib;/home/wheeltec/.autoware/install/op_ros_helpers/lib;/home/wheeltec/.autoware/install/ff_waypoint_follower/lib;/home/wheeltec/.autoware/install/op_simu/lib;/home/wheeltec/.autoware/install/op_planner/lib;/home/wheeltec/.autoware/install/op_utility/lib;/home/wheeltec/.autoware/install/vector_map_server/lib;/home/wheeltec/.autoware/install/naive_motion_predict/lib;/home/wheeltec/.autoware/install/map_file/lib;/home/wheeltec/.autoware/install/libvectormap/lib;/home/wheeltec/.autoware/install/lane_planner/lib;/home/wheeltec/.autoware/install/imm_ukf_pda_track/lib;/home/wheeltec/.autoware/install/decision_maker/lib;/home/wheeltec/.autoware/install/vector_map/lib;/home/wheeltec/.autoware/install/vector_map_msgs/lib;/home/wheeltec/.autoware/install/vectacam/lib;/home/wheeltec/.autoware/install/udon_socket/lib;/home/wheeltec/.autoware/install/twist_generator/lib;/home/wheeltec/.autoware/install/twist_gate/lib;/home/wheeltec/.autoware/install/twist_filter/lib;/home/wheeltec/.autoware/install/twist2odom/lib;/home/wheeltec/.autoware/install/tablet_socket/lib;/home/wheeltec/.autoware/install/runtime_manager/lib;/home/wheeltec/.autoware/install/mqtt_socket/lib;/home/wheeltec/.autoware/install/tablet_socket_msgs/lib;/home/wheeltec/.autoware/install/state_machine_lib/lib;/home/wheeltec/.autoware/install/sound_player/lib;/home/wheeltec/.autoware/install/sick_lms5xx/lib;/home/wheeltec/.autoware/install/sick_ldmrs_tools/lib;/home/wheeltec/.autoware/install/sick_ldmrs_driver/lib;/home/wheeltec/.autoware/install/sick_ldmrs_msgs/lib;/home/wheeltec/.autoware/install/sick_ldmrs_description/lib;/home/wheeltec/.autoware/install/points2image/lib;/home/wheeltec/.autoware/install/rosinterface/lib;/home/wheeltec/.autoware/install/rosbag_controller/lib;/home/wheeltec/.autoware/install/pure_pursuit/lib;/home/wheeltec/.autoware/install/emergency_handler/lib;/home/wheeltec/.autoware/install/autoware_health_checker/lib;/home/wheeltec/.autoware/install/as/lib;/home/wheeltec/.autoware/install/ros_observer/lib;/home/wheeltec/.autoware/install/range_vision_fusion/lib;/home/wheeltec/.autoware/install/pos_db/lib;/home/wheeltec/.autoware/install/pixel_cloud_fusion/lib;/home/wheeltec/.autoware/install/pc2_downsampler/lib;/home/wheeltec/.autoware/install/oculus_socket/lib;/home/wheeltec/.autoware/install/obj_db/lib;/home/wheeltec/.autoware/install/nmea_navsat/lib;/home/wheeltec/.autoware/install/ndt_tku/lib;/home/wheeltec/.autoware/install/ndt_cpu/lib;/home/wheeltec/.autoware/install/multi_lidar_calibrator/lib;/home/wheeltec/.autoware/install/microstrain_driver/lib;/home/wheeltec/.autoware/install/memsic_imu/lib;/home/wheeltec/.autoware/install/marker_downsampler/lib;/home/wheeltec/.autoware/install/map_tools/lib;/home/wheeltec/.autoware/install/map_tf_generator/lib;/home/wheeltec/.autoware/install/log_tools/lib;/home/wheeltec/.autoware/install/lidar_shape_estimation/lib;/home/wheeltec/.autoware/install/lidar_point_pillars/lib;/home/wheeltec/.autoware/install/lidar_naive_l_shape_detect/lib;/home/wheeltec/.autoware/install/lidar_fake_perception/lib;/home/wheeltec/.autoware/install/lidar_apollo_cnn_seg_detect/lib;/home/wheeltec/.autoware/install/libwaypoint_follower/lib;/home/wheeltec/.autoware/install/lgsvl_simulator_bridge/lib;/home/wheeltec/.autoware/install/lanelet2_extension/lib;/home/wheeltec/.autoware/install/kvaser/lib;/home/wheeltec/.autoware/install/kitti_box_publisher/lib;/home/wheeltec/.autoware/install/javad_navsat_driver/lib;/home/wheeltec/.autoware/install/image_processor/lib;/home/wheeltec/.autoware/install/graph_tools/lib;/home/wheeltec/.autoware/install/gnss_localizer/lib;/home/wheeltec/.autoware/install/gnss/lib;/home/wheeltec/.autoware/install/glviewer/lib;/home/wheeltec/.autoware/install/gazebo_world_description/lib;/home/wheeltec/.autoware/install/gazebo_imu_description/lib;/home/wheeltec/.autoware/install/gazebo_camera_description/lib;/home/wheeltec/.autoware/install/garmin/lib;/home/wheeltec/.autoware/install/freespace_planner/lib;/home/wheeltec/.autoware/install/fastvirtualscan/lib;/home/wheeltec/.autoware/install/ekf_localizer/lib;/home/wheeltec/.autoware/install/ds4_msgs/lib;/home/wheeltec/.autoware/install/ds4_driver/lib;/home/wheeltec/.autoware/install/decision_maker_panel/lib;/home/wheeltec/.autoware/install/data_preprocessor/lib;/home/wheeltec/.autoware/install/custom_msgs/lib;/home/wheeltec/.autoware/install/calibration_publisher/lib;/home/wheeltec/.autoware/install/autoware_system_msgs/lib;/home/wheeltec/.autoware/install/autoware_rviz_plugins/lib;/home/wheeltec/.autoware/install/autoware_quickstart_examples/lib;/home/wheeltec/.autoware/install/autoware_pointgrey_drivers/lib;/home/wheeltec/.autoware/install/autoware_driveworks_interface/lib;/home/wheeltec/.autoware/install/autoware_connector/lib;/home/wheeltec/.autoware/install/astar_search/lib;/home/wheeltec/.autoware/install/amathutils_lib/lib;/home/wheeltec/.autoware/install/autoware_msgs/lib;/home/wheeltec/.autoware/install/autoware_map_msgs/lib;/home/wheeltec/.autoware/install/autoware_launcher_rviz/lib;/home/wheeltec/.autoware/install/autoware_launcher/lib;/home/wheeltec/.autoware/install/autoware_lanelet2_msgs/lib;/home/wheeltec/.autoware/install/autoware_external_msgs/lib;/home/wheeltec/.autoware/install/autoware_driveworks_gmsl_interface/lib;/home/wheeltec/.autoware/install/autoware_config_msgs/lib;/home/wheeltec/.autoware/install/autoware_can_msgs/lib;/home/wheeltec/.autoware/install/autoware_build_flags/lib;/home/wheeltec/.autoware/install/autoware_bag_tools/lib;/home/wheeltec/.autoware/install/adi_driver/lib;/home/wheeltec/wheeltec_lidar/devel/lib;/home/wheeltec/cartographer_ws/devel_isolated/cartographer_rviz/lib;/home/wheeltec/cartographer_ws/install_isolated/lib;/home/wheeltec/wheeltec_robot/devel/lib;/opt/ros/noetic/lib)
      find_library(lib ${library}
        PATHS ${path}
        NO_DEFAULT_PATH NO_CMAKE_FIND_ROOT_PATH)
      if(lib)
        set(lib_path ${path})
        break()
      endif()
    endforeach()
    if(lib)
      _list_append_unique(stanley_controller_LIBRARY_DIRS ${lib_path})
      list(APPEND stanley_controller_LIBRARIES ${lib})
    else()
      # as a fall back for non-catkin libraries try to search globally
      find_library(lib ${library})
      if(NOT lib)
        message(FATAL_ERROR "Project '${PROJECT_NAME}' tried to find library '${library}'.  The library is neither a target nor built/installed properly.  Did you compile project 'stanley_controller'?  Did you find_package() it before the subdirectory containing its code is included?")
      endif()
      list(APPEND stanley_controller_LIBRARIES ${lib})
    endif()
  endif()
endforeach()

set(stanley_controller_EXPORTED_TARGETS "")
# create dummy targets for exported code generation targets to make life of users easier
foreach(t ${stanley_controller_EXPORTED_TARGETS})
  if(NOT TARGET ${t})
    add_custom_target(${t})
  endif()
endforeach()

set(depends "")
foreach(depend ${depends})
  string(REPLACE " " ";" depend_list ${depend})
  # the package name of the dependency must be kept in a unique variable so that it is not overwritten in recursive calls
  list(GET depend_list 0 stanley_controller_dep)
  list(LENGTH depend_list count)
  if(${count} EQUAL 1)
    # simple dependencies must only be find_package()-ed once
    if(NOT ${stanley_controller_dep}_FOUND)
      find_package(${stanley_controller_dep} REQUIRED NO_MODULE)
    endif()
  else()
    # dependencies with components must be find_package()-ed again
    list(REMOVE_AT depend_list 0)
    find_package(${stanley_controller_dep} REQUIRED NO_MODULE ${depend_list})
  endif()
  _list_append_unique(stanley_controller_INCLUDE_DIRS ${${stanley_controller_dep}_INCLUDE_DIRS})

  # merge build configuration keywords with library names to correctly deduplicate
  _pack_libraries_with_build_configuration(stanley_controller_LIBRARIES ${stanley_controller_LIBRARIES})
  _pack_libraries_with_build_configuration(_libraries ${${stanley_controller_dep}_LIBRARIES})
  _list_append_deduplicate(stanley_controller_LIBRARIES ${_libraries})
  # undo build configuration keyword merging after deduplication
  _unpack_libraries_with_build_configuration(stanley_controller_LIBRARIES ${stanley_controller_LIBRARIES})

  _list_append_unique(stanley_controller_LIBRARY_DIRS ${${stanley_controller_dep}_LIBRARY_DIRS})
  _list_append_deduplicate(stanley_controller_EXPORTED_TARGETS ${${stanley_controller_dep}_EXPORTED_TARGETS})
endforeach()

set(pkg_cfg_extras "")
foreach(extra ${pkg_cfg_extras})
  if(NOT IS_ABSOLUTE ${extra})
    set(extra ${stanley_controller_DIR}/${extra})
  endif()
  include(${extra})
endforeach()
