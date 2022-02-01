add_executable (latticeMeshElasticityTest1
        main.cpp
        )

target_include_directories(latticeMeshElasticityTest1
        PUBLIC
            ${USD_INC_DIR}
        )

target_link_libraries(
        latticeMeshElasticityTest1
                ${USD_LIB_DIR}/libusd_gf.so
                ${USD_LIB_DIR}/libusd_sdf.so
                ${USD_LIB_DIR}/libusd_tf.so
                ${USD_LIB_DIR}/libusd_usd.so
                ${USD_LIB_DIR}/libusd_usdGeom.so
                ${USD_LIB_DIR}/libusd_vt.so
        )
