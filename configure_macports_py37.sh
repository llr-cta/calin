cmake \ #--graphviz=calin_deps.dot \
    -DPYTHON_EXECUTABLE=/opt/local/bin/python3.7 \
    -DPYTHON_INCLUDE_DIR=/opt/local/Library/Frameworks/Python.framework/Versions/3.7/include/python3.7m \
    -DPYTHON_LIBRARY=/opt/local/Library/Frameworks/Python.framework/Versions/3.7/lib/python3.7/config-3.7m-darwin/libpython3.7.dylib \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCTA_CAMERASTOACT_DIR=$HOME/local \
    -DCMAKE_INSTALL_PREFIX=$HOME/local ..
