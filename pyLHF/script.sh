git clone https://github.com/yanurag994/lhf.git /lhf
mkdir /lhf/build && cd /lhf/build
cmake .. -DCMAKE_BUILDD_TYPE=Release && make -j

cp -f /lhf/build/LHFmain/libLHFlib.so /lhf/pyLHF/src/lhf/

cd /lhf/pyLHF/src/
"/opt/python/cp310-cp310/bin/pip3" wheel --no-deps -w wheelhouse/ . --no-cache

cd /lhf/pyLHF/src/wheelhouse
auditwheel show "lhf-2.0.1-py3-none-any.whl"
auditwheel repair "lhf-2.0.1-py3-none-any.whl" --plat "manylinux_2_24_x86_64"