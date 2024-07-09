# Copyright (c) 2019, The Science and Technology Facilities Council (STFC)
# All rights reserved.
# copy python example

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/solve.py
	${CMAKE_CURRENT_BINARY_DIR}/solve.py
 	COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/solve_weights.py
	${CMAKE_CURRENT_BINARY_DIR}/solve_weights.py
 	COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/solve_box.py
	${CMAKE_CURRENT_BINARY_DIR}/solve_box.py
 	COPYONLY)
