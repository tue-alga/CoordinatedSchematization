function(CREATE_COPY_DLL_CUSTOM_COMMAND TARGET_NAME)
# Acquire all linked libraries
get_target_property(lls ${TARGET_NAME} LINK_LIBRARIES)
# Loop over them
foreach(LL ${lls})

	message(STATUS "CLI Lib: ${LL}")
    get_target_property(_isImported ${LL} IMPORTED)
    
    get_target_property(_llType ${LL} TYPE)
    # Ignore interface libraries, they cannot have dlls?
    if(${_llType} STREQUAL "INTERFACE_LIBRARY")
        continue()
    endif()

	message(STATUS "Imported: ${_isImported}")
	get_target_property(_data ${LL} IMPORTED_LOCATION_RELEASE)
	message(STATUS "-Release: ${_data}")
	get_target_property(_data ${LL} IMPORTED_LOCATION_DEBUG)
	message(STATUS "-Debug: ${_data}")
	get_target_property(_dataLL ${LL} IMPORTED_LINK_DEPENDENT_LIBRARIES)
	message(STATUS "-Dependent (dep): ${_dataLL}")
	get_target_property(_EXTRA_LIBS ${LL} INTERFACE_LINK_LIBRARIES)
	message(STATUS "-Dependent (inter): ${_EXTRA_LIBS}")
    # Determine how to copy: if the target is DLL, copy it, otherwise ignore.
    # TODO: find if PDB's are in the neighbourhood and copy those.

	add_custom_command(TARGET ${PROJ_NAME} POST_BUILD 
		COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:${LL}> $<TARGET_FILE_DIR:${PROJ_NAME}>
	)
	# TODO: move to separate function
	foreach(ELL ${_EXTRA_LIBS}) 
		get_target_property(_EXTRA_LIBS_SUB ${ELL} INTERFACE_LINK_LIBRARIES)
		message(STATUS "--- ${_EXTRA_LIBS_SUB}")
		message(STATUS "-- Extra dep: ${ELL}")
		add_custom_command(TARGET ${PROJ_NAME} POST_BUILD 
			COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:${ELL}> $<TARGET_FILE_DIR:${PROJ_NAME}>
		)
	endforeach()
	
endforeach()
endfunction()