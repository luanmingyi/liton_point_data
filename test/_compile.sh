cmake -H. -B_build_unix_debug -DCMAKE_PREFIX_PATH=$CPPLIB -DCMAKE_MODULE_PATH=$CPPLIB -DCMAKE_BUILD_TYPE=Debug
cmake --build _build_unix_debug
cmake -H. -B_build_unix_release -DCMAKE_PREFIX_PATH=$CPPLIB -DCMAKE_MODULE_PATH=$CPPLIB -DCMAKE_BUILD_TYPE=Release
cmake --build _build_unix_release
