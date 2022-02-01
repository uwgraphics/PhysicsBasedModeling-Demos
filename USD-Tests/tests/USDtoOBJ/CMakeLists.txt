add_executable (USDtoOBJ
        main.cpp
        )

target_include_directories(USDtoOBJ
        PUBLIC
            ${USD_INC_DIR}
        )

target_link_libraries(
        USDtoOBJ
                ${USD_LIB_DIR}/libusd_gf.so
                ${USD_LIB_DIR}/libusd_sdf.so
                ${USD_LIB_DIR}/libusd_tf.so
                ${USD_LIB_DIR}/libusd_usd.so
                ${USD_LIB_DIR}/libusd_usdGeom.so
                ${USD_LIB_DIR}/libusd_vt.so
        )
