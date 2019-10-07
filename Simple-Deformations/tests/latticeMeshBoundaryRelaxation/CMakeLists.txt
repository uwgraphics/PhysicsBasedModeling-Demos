add_executable (latticeMeshBoundaryRelaxation
        main.cpp
        )

target_include_directories(latticeMeshBoundaryRelaxation
        PUBLIC
            ${USD_INC_DIR}
        )

target_link_libraries(
        latticeMeshBoundaryRelaxation
                ${USD_LIB_DIR}/libgf.so
                ${USD_LIB_DIR}/libsdf.so
                ${USD_LIB_DIR}/libtf.so
                ${USD_LIB_DIR}/libusd.so
                ${USD_LIB_DIR}/libusdGeom.so
                ${USD_LIB_DIR}/libvt.so
        )
