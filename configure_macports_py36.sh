cmake \ #--graphviz=calin_deps.dot \
    -DPYTHON_EXECUTABLE=/opt/local/bin/python3.6 \
    -DPYTHON_INCLUDE_DIR=/opt/local/Library/Frameworks/Python.framework/Versions/3.6/include/python3.6m \
    -DPYTHON_LIBRARY=/opt/local/Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/config-3.6m-darwin/libpython3.6.dylib \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCTA_CAMERASTOACT_DIR=$HOME/local \
    -DCMAKE_INSTALL_PREFIX=$HOME/local ..
