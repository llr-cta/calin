cmake \
    -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
    -DPYTHON_EXECUTABLE=/opt/local/bin/python3.10 \
    -DPYTHON_INCLUDE_DIR=/opt/local/Library/Frameworks/Python.framework/Versions/3.10/include/python3.10 \
    -DPYTHON_LIBRARY=/opt/local/Library/Frameworks/Python.framework/Versions/3.10/lib/python3.10/config-3.10-darwin/libpython3.10.dylib \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCTA_CAMERASTOACT_DIR=$HOME/local \
    -DCMAKE_INSTALL_PREFIX=$HOME/local ..
