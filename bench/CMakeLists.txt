
function(ADD_BENCH_TARGET TARGET POET_BENCH_LIST RT_FILES OUT_PATH)
    add_custom_target(${TARGET})

    set(bench_install_dir share/poet/${OUT_PATH})
    foreach(BENCH_FILE ${${POET_BENCH_LIST}})
        get_filename_component(BENCH_NAME ${BENCH_FILE} NAME_WE)
        set(OUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/${BENCH_NAME}.rds)

        add_custom_command(TARGET ${TARGET}
            COMMAND $<TARGET_FILE:poet_init> -o ${OUT_FILE} -s ${CMAKE_CURRENT_SOURCE_DIR}/${BENCH_FILE}
            COMMENT "Running poet_init on ${BENCH_FILE}"
            DEPENDS poet_init
            VERBATIM
        )

        install(FILES ${OUT_FILE} DESTINATION ${bench_install_dir})
    endforeach(BENCH_FILE ${${POET_BENCH_LIST}})

    # install all ADD_FILES to the same location
    install(FILES ${${RT_FILES}} DESTINATION ${bench_install_dir})
endfunction()


# define target name 
set(BENCHTARGET benchmarks)

add_custom_target(${BENCHTARGET} ALL)

add_subdirectory(barite)
add_subdirectory(dolo)
add_subdirectory(surfex)

