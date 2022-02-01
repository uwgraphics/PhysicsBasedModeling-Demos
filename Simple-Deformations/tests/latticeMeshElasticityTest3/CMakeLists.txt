add_executable (latticeMeshElasticityTest3
        main.cpp
        )

set(EIGEN3_INC_DIR /usr/include/eigen3)

target_include_directories(latticeMeshElasticityTest3
        PUBLIC
            ${USD_INC_DIR}
            ${EIGEN3_INC_DIR}
        )

target_link_libraries(
        latticeMeshElasticityTest3
                ${USD_LIB_DIR}/libusd_gf.so
                ${USD_LIB_DIR}/libusd_sdf.so
                ${USD_LIB_DIR}/libusd_tf.so
                ${USD_LIB_DIR}/libusd_usd.so
                ${USD_LIB_DIR}/libusd_usdGeom.so
                ${USD_LIB_DIR}/libusd_vt.so
        )
