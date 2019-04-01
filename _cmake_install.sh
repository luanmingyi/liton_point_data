cmake -H. -B_build_unix -DCMAKE_INSTALL_PREFIX=$CPPLIB/point_data
cmake --build _build_unix --config Release --target install
