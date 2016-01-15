cmake \
    -DPYTHON_EXECUTABLE=/opt/local/bin/python3.4 \
    -DPYTHON_INCLUDE_DIR=/opt/local/Library/Frameworks/Python.framework/Versions/3.4/include/python3.4m \
    -DPYTHON_LIBRARY=/opt/local/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/config-3.4m/libpython3.4.dylib \
    -DCMAKE_BUILD_TYPE=Release \
    -DCTA_CAMERASTOACT_DIR=$HOME/local \
    -DCMAKE_INSTALL_PREFIX=$HOME/local ..
