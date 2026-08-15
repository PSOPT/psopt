cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0057 NEW)

# Verify that each QP backend plugin exports only the entry points of the plugin ABI.
#
# The plugins exist so that the ordering and factorisation code inside them cannot
# collide with another backend's. That guarantee rests on hidden visibility and on
# --exclude-libs, both of which are easy to lose to a compiler flag or a link-order
# change, and whose loss shows up as a crash in an unrelated part of the program. So it
# is checked rather than assumed.

file(GLOB plugins "${PLUGIN_DIR}/libpsopt_qp_*.so")
if(NOT plugins)
    message(FATAL_ERROR "no QP backend plugins found in ${PLUGIN_DIR}")
endif()

set(allowed "psopt_qp_solve" "psopt_qp_name" "psopt_qp_abi_version")

foreach(plugin ${plugins})
    execute_process(COMMAND ${NM} --dynamic --defined-only --format=posix ${plugin}
                    OUTPUT_VARIABLE dump RESULT_VARIABLE rc)
    if(NOT rc EQUAL 0)
        message(FATAL_ERROR "could not read the dynamic symbols of ${plugin}")
    endif()

    string(REPLACE "\n" ";" lines "${dump}")
    set(unexpected "")
    foreach(line ${lines})
        if(line MATCHES "^([^ ]+) ([A-Za-z])")
            set(symbol "${CMAKE_MATCH_1}")
            set(kind   "${CMAKE_MATCH_2}")
            # Only strong definitions matter. Weak and vague-linkage symbols (W, V)
            # are the C++ ABI's shared entities -- typeinfo, vtables, inline
            # instantiations -- which are meant to be merged across objects and would
            # break exception handling if they were hidden; the ProxQP plugin exports
            # three of them for std::bad_optional_access. What the isolation is for is
            # the strong C definitions a factorisation library provides, amd_order and
            # its relatives. The linker's own bookkeeping symbols are not definitions of
            # anything a QP library provides.
            if(kind MATCHES "^[TDBR]$"
               AND NOT symbol IN_LIST allowed
               AND NOT symbol MATCHES "^(_init|_fini|__bss_start|_edata|_end|__gmon_start__)$")
                list(APPEND unexpected "${symbol} (${kind})")
            endif()
        endif()
    endforeach()

    get_filename_component(base "${plugin}" NAME)
    if(unexpected)
        list(LENGTH unexpected count)
        string(REPLACE ";" "\n    " pretty "${unexpected}")
        message(FATAL_ERROR
            "${base} exports ${count} symbol(s) beyond the plugin ABI:\n    ${pretty}\n"
            "Hidden visibility or --exclude-libs has been lost. Exported symbols from a "
            "QP library's linear algebra can bind to another backend's and corrupt it.")
    endif()
    message(STATUS "${base}: exports only the plugin ABI")
endforeach()
