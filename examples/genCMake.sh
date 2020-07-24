#!/bin/sh
set -uxe
function create_cmake() {
    local f="CMakeLists-generator.txt"
    for d in $(ls -d */ | tr -d '/'); do
        cat "$f" | sed -e "s/__MYDIR__/$d/g" > "$d/CMakeList.txt"
    done
}
create_cmake
