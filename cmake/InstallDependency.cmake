cmake_minimum_required(VERSION 3.16)
# Install a dependency via FetchContent at configuration time

function(INSTALL_DEPENDENCY)
    set(options CLEAN)
    set(oneValueArgs FETCHCONTENT_NAME INSTALL_DIRECTORY)
    set(multiValueArgs INSTALL_CONFIGURATIONS COMPILE_FLAGS LINK_FLAGS VARIABLE_DEFINES)
    cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}"
                        "${multiValueArgs}" ${ARGN} )

    set(_LOG_PREFIX "[InstallDependency]")
    message(STATUS "${_LOG_PREFIX} Installing ${ARGS_FETCHCONTENT_NAME}")
    if(NOT DEFINED ${ARGS_FETCHCONTENT_NAME}_POPULATED)
        FetchContent_Populate(${ARGS_FETCHCONTENT_NAME})
    endif()
    # Run CMake command
    FetchContent_GetProperties(${ARGS_FETCHCONTENT_NAME}
        SOURCE_DIR _SRC_DIR
        BINARY_DIR _BIN_DIR
    )
    if(NOT ARGS_COMPILE_FLAGS)
        set(ARGS_COMPILE_FLAGS "")
    endif()
    
    if(NOT ARGS_LINK_FLAGS)
        set(ARGS_LINK_FLAGS "")
    endif()
    if(NOT ARGS_VARIABLE_DEFINES)
        set(ARGS_VARIABLE_DEFINES "")
    else()
        list(TRANSFORM ARGS_VARIABLE_DEFINES PREPEND "-D")
        list(JOIN ARGS_VARIABLE_DEFINES " " ARGS_VARIABLE_DEFINES)
    endif()
    # Forward toolchain file if specified
    if(DEFINED CMAKE_TOOLCHAIN_FILE)
        message(STATUS "Setting toolchain for build: ${CMAKE_TOOLCHAIN_FILE}")
        set(_TOOLCHAIN -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE})    
    else()
        message(STATUS "No toolchain specified")
        set(_TOOLCHAIN "")
    endif()
    set(_EXTRA_ARGS)
    if(CMAKE_GENERATOR_PLATFORM)
        list(APPEND _EXTRA_ARGS "-A" "${CMAKE_GENERATOR_PLATFORM}")
    endif()
    #list(JOIN _EXTRA_ARGS " " _EXTRA_ARGS)


    if(ARGS_CLEAN)
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E remove_directory ${_BIN_DIR}
            COMMAND_ECHO STDOUT
            RESULT_VARIABLE _EXIT_CODE
        )
        if(_EXIT_CODE EQUAL "1")
            message( FATAL_ERROR "${_LOG_PREFIX} execute_process failed for clearing directory ${_BIN_DIR}")
        endif()
    endif()

    execute_process(
        COMMAND ${CMAKE_COMMAND} -S ${_SRC_DIR} -B ${_BIN_DIR} -G ${CMAKE_GENERATOR} -DCMAKE_INSTALL_PREFIX=${ARGS_INSTALL_DIRECTORY}
        ${ARGS_VARIABLE_DEFINES} ${_TOOLCHAIN} ${_EXTRA_ARGS}
        COMMAND_ECHO STDOUT
        RESULT_VARIABLE _EXIT_CODE
    )
    if(_EXIT_CODE EQUAL "1")
        message( FATAL_ERROR "${_LOG_PREFIX} execute_process failed for CMake config of ${ARGS_FETCHCONTENT_NAME}")
    endif()
    execute_process(
        COMMAND ${CMAKE_COMMAND} --build ${_BIN_DIR} --target install
        COMMAND_ECHO STDOUT
        RESULT_VARIABLE _EXIT_CODE
    )
    if(_EXIT_CODE EQUAL "1")
        message( FATAL_ERROR "${_LOG_PREFIX} execute_process failed for build of ${ARGS_FETCHCONTENT_NAME}")
    endif()
endfunction()
