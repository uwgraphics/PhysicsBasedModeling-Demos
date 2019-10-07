add_executable (latticeMesh
        main.cpp
        )

target_include_directories(latticeMesh
        PUBLIC
            ${USD_INC_DIR}
        )

target_link_libraries(
        latticeMesh
                ${USD_LIB_DIR}/libgf.so
                ${USD_LIB_DIR}/libsdf.so
                ${USD_LIB_DIR}/libtf.so
                ${USD_LIB_DIR}/libusd.so
                ${USD_LIB_DIR}/libusdGeom.so
                ${USD_LIB_DIR}/libvt.so
        )
