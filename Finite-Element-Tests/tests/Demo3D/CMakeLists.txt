add_executable (Demo3D
        main.cpp
        )

# if CMake cannot use find_package for Eigen, use the explicit include path instead
find_package(Eigen3)
# set(EIGEN3_INCLUDE_DIRS /usr/include/eigen3/)

target_include_directories(Demo3D
        PUBLIC
            ${USD_INC_DIR}
            ${EIGEN3_INCLUDE_DIRS}
            ${PYTHON_INCLUDE_PATH}
        )

target_link_libraries(
        Demo3D
                ${USD_LIB_DIR}/libusd_gf.so
                ${USD_LIB_DIR}/libusd_sdf.so
                ${USD_LIB_DIR}/libusd_tf.so
                ${USD_LIB_DIR}/libusd_usd.so
                ${USD_LIB_DIR}/libusd_usdGeom.so
                ${USD_LIB_DIR}/libusd_vt.so
        )
