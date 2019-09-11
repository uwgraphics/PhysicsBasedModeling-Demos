add_executable (latticeMeshElasticityTest2
        main.cpp
        )

find_package(PythonLibs)

target_include_directories(latticeMeshElasticityTest2
        PUBLIC
            ${USD_INC_DIR}
            ${PYTHON_INCLUDE_PATH}
        )

target_link_libraries(
        latticeMeshElasticityTest2
                ${USD_LIB_DIR}/libgf.so
                ${USD_LIB_DIR}/libsdf.so
                ${USD_LIB_DIR}/libtf.so
                ${USD_LIB_DIR}/libusd.so
                ${USD_LIB_DIR}/libusdGeom.so
                ${USD_LIB_DIR}/libvt.so
                ${USD_LIB_DIR}/libboost_python.so # todo find library
                ${PYTHON_LIBRARY}
        )
