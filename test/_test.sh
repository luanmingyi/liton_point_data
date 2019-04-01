cd _build_unix_debug
ctest -R Test_
cd ..

cd _build_unix_release
ctest -R Example_
cd ..
