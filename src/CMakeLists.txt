set(EXECUTABLE tb_model)

add_executable(${EXECUTABLE} 
	main.f90
  parallel.f90
  communication.f90
  inputoutput.f90
  math.f90
  constants.f90
  electronic_system.f90
  )

target_link_libraries(${EXECUTABLE} lapack)

#set_target_properties(${EXECUTABLE}
#  PROPERTIES
#  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
#  )


