# cmaple


...



## Improving runtime even further:
Cmaple is optimized for speed. To make it even faster, you can swap in another high-performance allocator like [jemalloc](https://github.com/jemalloc/jemalloc.git) .
This will give another 5-10% speedup, depending on workload and hardware.

### Install jemalloc either via

#### package manager

#### ... or manually (Linux example)

```
git clone https://github.com/jemalloc/jemalloc.git
cd jemalloc
export je_build=`pwd`/build   
./autogen.sh
./configure --prefix=${je_build}
make -j20
make install
## remember this path!
echo "remember to put '${je_build}/bin' in your PATH" to make 'jemalloc-config' known
## we will do it here once, but you need to make this permanent:
export PATH="${je_build}/bin:${PATH}"
```

 ### Use jemalloc 
 
```
## run cmaple with preloaded jemalloc
LD_PRELOAD=`jemalloc-config --libdir`/libjemalloc.so.`jemalloc-config --revision` ./cmaple <more args here>
```


