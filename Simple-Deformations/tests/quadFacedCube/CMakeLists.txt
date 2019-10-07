add_executable (quadFacedCube
        main.cpp
        )

target_include_directories(quadFacedCube
        PUBLIC
            ${USD_INC_DIR}
        )

target_link_libraries(
        quadFacedCube
                ${USD_LIB_DIR}/libgf.so
                ${USD_LIB_DIR}/libsdf.so
                ${USD_LIB_DIR}/libtf.so
                ${USD_LIB_DIR}/libusd.so
                ${USD_LIB_DIR}/libusdGeom.so
                ${USD_LIB_DIR}/libvt.so
        )
