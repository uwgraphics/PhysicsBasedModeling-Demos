add_executable (latticeMeshBoundaryRelaxation
        main.cpp
        )

target_include_directories(latticeMeshBoundaryRelaxation
        PUBLIC
            ${USD_INC_DIR}
        )

target_link_libraries(
        latticeMeshBoundaryRelaxation
                ${USD_LIB_DIR}/libusd_gf.so
                ${USD_LIB_DIR}/libusd_sdf.so
                ${USD_LIB_DIR}/libusd_tf.so
                ${USD_LIB_DIR}/libusd_usd.so
                ${USD_LIB_DIR}/libusd_usdGeom.so
                ${USD_LIB_DIR}/libusd_vt.so
        )
