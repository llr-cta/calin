cmake3 -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DPYTHON_EXECUTABLE=/usr/bin/python3.4 \
  -DPYTHON_INCLUDE_DIR=/usr/include/python3.4m \
  -DPYTHON_LIBRARY=/usr/lib64/libpython3.so \
  -DCMAKE_PREFIX_PATH=/opt/exp_soft/vo.cta.in2p3.fr/local \
  -DCMAKE_LIBRARY_PATH=/usr/lib64 \
  -DCMAKE_INSTALL_PREFIX=$HOME/local \
  ..
