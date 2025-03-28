cmake \
    -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
    -DPython3_EXECUTABLE=/opt/local/bin/python3.10 \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCTA_ADH_APIS_DIR=$HOME/local \
    -DNO_CORSIKA8=ON \
    -DCMAKE_INSTALL_PREFIX=$HOME/local ..
