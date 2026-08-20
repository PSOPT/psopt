cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0057 NEW)

# Verify that each QP backend plugin exports only the entry points of the plugin ABI.
#
# The plugins exist so that the ordering and factorisation code inside them cannot
# collide with another backend's. That guarantee rests on hidden visibility and, on
# ELF, on --exclude-libs; on Mach-O it rests on -exported_symbols_list. All of those
# are easy to lose to a compiler flag or a link-order change, and their loss shows up
# as a crash in an unrelated part of the program. So it is checked rather than assumed.
#
# CMake variables expected: NM, PLUGIN_DIR, and PLUGIN_FORMAT ("macho" or "elf").

# CMake gives MODULE libraries the suffix .so on macOS as well as on Linux
# (CMAKE_SHARED_MODULE_SUFFIX in Platform/Darwin.cmake), so one glob serves both.
file(GLOB plugins "${PLUGIN_DIR}/libpsopt_qp_*.so")
if(NOT plugins)
    message(FATAL_ERROR "no QP backend plugins found in ${PLUGIN_DIR}")
endif()

set(allowed "psopt_qp_solve" "psopt_qp_name" "psopt_qp_abi_version"
                "psopt_qp_environment_ok")

# The first three are not merely permitted, they are required: a plugin missing one of
# them loads and then fails at the first subproblem. psopt_qp_environment_ok is
# deliberately absent from this list -- it is optional, and only the GALAHAD plugin
# provides it.
set(required "psopt_qp_abi_version" "psopt_qp_name" "psopt_qp_solve")

foreach(plugin ${plugins})
    # nm's spelling differs by object format.
    #
    # ELF: --dynamic reads the dynamic symbol table, which is exactly the set that
    # matters, and --format=posix gives "name type value size".
    #
    # Mach-O: there is no dynamic symbol table to ask for, and -U has meant opposite
    # things in Apple's classic nm and in the llvm-nm that replaced it, so neither is
    # used. Plain -g lists external symbols in BSD format, "value type name", and
    # undefined ones carry type U, which the type filter below already rejects. Mach-O
    # also prefixes every C symbol with an underscore; it is stripped so that one
    # "allowed" list serves both formats.
    if(PLUGIN_FORMAT STREQUAL "macho")
        execute_process(COMMAND ${NM} -g ${plugin}
                        OUTPUT_VARIABLE dump RESULT_VARIABLE rc)
    else()
        execute_process(COMMAND ${NM} --dynamic --defined-only --format=posix ${plugin}
                        OUTPUT_VARIABLE dump RESULT_VARIABLE rc)
    endif()
    if(NOT rc EQUAL 0)
        message(FATAL_ERROR "could not read the exported symbols of ${plugin}")
    endif()

    string(REPLACE "\n" ";" lines "${dump}")
    set(unexpected "")
    set(defined "")
    foreach(line ${lines})
        set(symbol "")
        if(PLUGIN_FORMAT STREQUAL "macho")
            # "0000000000003a10 T _psopt_qp_solve", or "                 U _malloc"
            if(line MATCHES "^[0-9a-fA-F]* *([A-Za-z]) +_?(.+)$")
                set(kind   "${CMAKE_MATCH_1}")
                set(symbol "${CMAKE_MATCH_2}")
            endif()
        else()
            if(line MATCHES "^([^ ]+) ([A-Za-z])")
                set(symbol "${CMAKE_MATCH_1}")
                set(kind   "${CMAKE_MATCH_2}")
            endif()
        endif()

        if(symbol)
            # Only strong definitions matter. Weak and vague-linkage symbols (W, V,
            # and Mach-O's S) are the C++ ABI's shared entities -- typeinfo, vtables,
            # inline instantiations -- which are meant to be merged across objects and
            # would break exception handling if they were hidden; the ProxQP plugin
            # exports three of them for std::bad_optional_access. What the isolation is
            # for is the strong C definitions a factorisation library provides,
            # amd_order and its relatives. The linker's own bookkeeping symbols are not
            # definitions of anything a QP library provides.
            if(kind MATCHES "^[TDBR]$")
                list(APPEND defined "${symbol}")
            endif()

            if(kind MATCHES "^[TDBR]$"
               AND NOT symbol IN_LIST allowed
               AND NOT symbol MATCHES "^(_init|_fini|__bss_start|_edata|_end|__gmon_start__)$"
               # Mach-O's own header symbols, as they read after the single leading
               # underscore has been stripped above -- they carry two to begin with.
               AND NOT symbol MATCHES "^(_mh_dylib_header|_mh_bundle_header|_mh_execute_header|dyld_stub_binder)$")
                list(APPEND unexpected "${symbol} (${kind})")
            endif()
        endif()
    endforeach()

    get_filename_component(base "${plugin}" NAME)
    if(unexpected)
        list(LENGTH unexpected count)
        string(REPLACE ";" "\n    " pretty "${unexpected}")
        if(PLUGIN_FORMAT STREQUAL "macho")
            set(how "-exported_symbols_list or hidden visibility has been lost")
        else()
            set(how "Hidden visibility or --exclude-libs has been lost")
        endif()
        message(FATAL_ERROR
            "${base} exports ${count} symbol(s) beyond the plugin ABI:\n    ${pretty}\n"
            "${how}. Exported symbols from a QP library's linear algebra can bind to "
            "another backend's and corrupt it.")
    endif()

    # The other half of the same question. Restricting what a plugin exports is only
    # half of keeping the export list and the plugin in step; the list can also name a
    # symbol the plugin has stopped defining, or fail to name one it still needs. On
    # Mach-O the second of those is a link error already, because an exported-symbols
    # list doubles as a set of initial undefines; on ELF nothing catches it until the
    # loader does, at the first subproblem of a user's run.
    set(missing "")
    foreach(req ${required})
        if(NOT req IN_LIST defined)
            list(APPEND missing "${req}")
        endif()
    endforeach()
    if(missing)
        string(REPLACE ";" ", " pretty "${missing}")
        message(FATAL_ERROR
            "${base} does not define the required plugin entry point(s): ${pretty}. "
            "A plugin missing one of these loads and then fails at the first QP "
            "subproblem, with a message about missing entry points.")
    endif()

    message(STATUS "${base}: exports only the plugin ABI")
endforeach()
