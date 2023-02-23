target=$1
kspThreads=${2:-1}
dsThreads=${3:-1}

echo "ksp threads: ${kspThreads}"
echo "ds threads: ${dsThreads}"

buildDir="build-release"

rm -rf ${buildDir}
mkdir -p ${buildDir}
cd ${buildDir}
cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DTHREADS_KSP=${kspThreads} -DTHREADS_DS=${dsThreads} ../
cmake --build . -j 8 --target ${target}
