
function(ADD_BENCH_TARGET TARGET POET_BENCH_LIST RT_FILES OUT_PATH)
    set(bench_install_dir share/poet/${OUT_PATH})

    # create empty list 
    set(OUT_FILES_LIST "")

    foreach(BENCH_FILE ${${POET_BENCH_LIST}})
        get_filename_component(BENCH_NAME ${BENCH_FILE} NAME_WE)
        set(OUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/${BENCH_NAME}.rds)

        add_custom_command(
            OUTPUT ${OUT_FILE}
            COMMAND $<TARGET_FILE:poet_init> -o ${OUT_FILE} -s ${CMAKE_CURRENT_SOURCE_DIR}/${BENCH_FILE}
            COMMENT "Running poet_init on ${BENCH_FILE}"
            DEPENDS poet_init ${CMAKE_CURRENT_SOURCE_DIR}/${BENCH_FILE}
            VERBATIM
            COMMAND_EXPAND_LISTS
        )

        list(APPEND OUT_FILES_LIST ${OUT_FILE})

    endforeach(BENCH_FILE ${${POET_BENCH_LIST}})
        
    add_custom_target(
        ${TARGET}
        DEPENDS ${OUT_FILES_LIST})
            
    install(FILES ${OUT_FILES_LIST} DESTINATION ${bench_install_dir})

    # install all ADD_FILES to the same location
    install(FILES ${${RT_FILES}} DESTINATION ${bench_install_dir})
endfunction()


# define target name 
set(BENCHTARGET benchmarks)

add_custom_target(${BENCHTARGET} ALL)

add_subdirectory(barite)
add_subdirectory(dolo)
add_subdirectory(surfex)

