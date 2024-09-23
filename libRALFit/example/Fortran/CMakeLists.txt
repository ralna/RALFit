# Copyright (c) 2019, The Science and Technology Facilities Council (STFC)
# All rights reserved.
# Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
# Examples build

add_executable (nlls_example nlls_example.F90)
target_link_libraries(nlls_example ral_nlls)

add_executable (nlls_example2 nlls_example2.F90)
target_link_libraries(nlls_example2 ral_nlls)

add_executable (jacobian jacobian.F90)
target_link_libraries(jacobian ral_nlls)

add_executable (lanczos Lanczos.F90)
target_link_libraries(lanczos ral_nlls)

add_executable (LanczosBox LanczosBox.F90)
target_compile_options(LanczosBox PRIVATE "$<$<CONFIG:DEBUG>:-Wno-unused-dummy-argument>")
target_link_libraries(LanczosBox ral_nlls)

add_executable (LanczosFD LanczosFD.F90)
target_link_libraries(LanczosFD ral_nlls)

add_executable (LanczosFD2 LanczosFD2.F90)
target_link_libraries(LanczosFD2 ral_nlls)

add_executable (stp_x0sol stp_x0sol.F90)
target_compile_options(stp_x0sol PRIVATE "$<$<CONFIG:DEBUG>:-Wno-unused-dummy-argument>")
target_link_libraries(stp_x0sol ral_nlls)
