cmake_minimum_required(VERSION 3.4.1)
#file(GLOB libf2c_sources src/main/cpp/clapack/F2CLIBS/*.c)
#file(GLOB blas_sources src/main/cpp/clapack/BLAS/SRC/*.c)
#set(libf2c ${libf2c_sources})
#list(FILTER blas_sources EXCLUDE REGEX ".*Axb.c$")
#list(FILTER blas_sources EXCLUDE REGEX ".*Axb_core.c$")
#set(blas ${blas_sources})
add_library(roomshape SHARED
            # toplevel
            src/main/cpp/roomshape.cpp
            # lrs
            src/main/cpp/lrslib/lrslib.c
            src/main/cpp/lrslib/lrslong.c)
            # libf2c
            #${libf2c}
            # BLAS
            #${blas})
add_definitions(-O2
                -DLRSLONG
                -Dskip_f2c_Undefs
                -DNO_LONG_LONG)
find_library(log-lib
             log)
target_link_libraries(roomshape
                      ${log-lib})