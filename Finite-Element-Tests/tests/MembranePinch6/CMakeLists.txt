add_executable (MembranePinch6
        main.cpp
        )

# if CMake cannot use find_package for Eigen, use the explicit include path instead
find_package(Eigen3)
# set(EIGEN3_INCLUDE_DIRS /usr/include/eigen3/)

target_include_directories(MembranePinch6
        PUBLIC
            ${USD_INC_DIR}
            ${EIGEN3_INCLUDE_DIRS}
            ${PYTHON_INCLUDE_PATH}
        )

target_link_libraries(
        MembranePinch6
                ${USD_LIB_DIR}/libgf.so
                ${USD_LIB_DIR}/libsdf.so
                ${USD_LIB_DIR}/libtf.so
                ${USD_LIB_DIR}/libusd.so
                ${USD_LIB_DIR}/libusdGeom.so
                ${USD_LIB_DIR}/libvt.so
        )
