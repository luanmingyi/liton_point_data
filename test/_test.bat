cd _build_win32
ctest -C Debug -R Test_*
ctest -C Release -R Example_*
cd ..
