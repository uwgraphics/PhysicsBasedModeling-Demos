add_executable (animateAttribute
        main.cpp
        )


target_include_directories(animateAttribute
        PUBLIC
            ${USD_INC_DIR}
        )

target_link_libraries(
        animateAttribute
                ${USD_LIB_DIR}/libusd_gf.so
                ${USD_LIB_DIR}/libusd_sdf.so
                ${USD_LIB_DIR}/libusd_tf.so
                ${USD_LIB_DIR}/libusd_usd.so
                ${USD_LIB_DIR}/libusd_usdGeom.so
                ${USD_LIB_DIR}/libusd_vt.so
        )
