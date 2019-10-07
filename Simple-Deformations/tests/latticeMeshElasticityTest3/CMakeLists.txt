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
                ${USD_LIB_DIR}/libgf.so
                ${USD_LIB_DIR}/libsdf.so
                ${USD_LIB_DIR}/libtf.so
                ${USD_LIB_DIR}/libusd.so
                ${USD_LIB_DIR}/libusdGeom.so
                ${USD_LIB_DIR}/libvt.so
        )
