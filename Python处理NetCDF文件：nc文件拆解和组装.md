# Python对NC文件的拆解组装

功能需求：对已有的NetCDF(nc)格式数据进行接续，选择特地变量并组装成一个新的nc文件，原样输出且方便作图。
实现思路：利用python中的netCDF4的包变量解析，利用numpy进行数据拷贝和赋值
特色：利用hdf5自带的压缩算法，对数据压缩即对每个chunk个用zlib压缩
如果大家有兴趣就请看下去吧!



-------------------
大象放进冰箱的三部曲。
**打开冰箱**

```
from netCDF4 import *
import numpy as np
# get NC from local inputFile
nc = Dataset(inputFileName)
# 新建一个本地文件，格式是netcdf4_classic
newnc = Dataset(outputFileName, "w", format="NETCDF4_CLASSIC")
```
pat：这里使用的格式是netcdf4_classic,自带压缩，同样netcdf4也支持，其他的不支持
**把大象放进冰箱**

```
#第一部分拷贝Dimensions
# get and copy dimensions to new file
ncdimesions = nc.dimensions
# get sample_dim
if ("sample" in ncdimesions):
    sample_dim = ncdimesions["sample"]
    newncdim_sample = newnc.createDimension(sample_dim.name, sample_dim.size)
# get ddm_dim
if ("ddm" in ncdimesions):
    ddm_dim = ncdimesions["ddm"]
    newncdim_ddm = newnc.createDimension(ddm_dim.name, ddm_dim.size)
#第二部分拷贝变量的属性，维度，和数据
# !spacecraft_num!
spacecraft_num = ncvariables["spacecraft_num"]
spacecraft_num_data = spacecraft_num[:]
# copy spacecraft_num to newnc
spacecraft_num_fill_value = spacecraft_num.getncattr("_FillValue")
#新建变量，当zlib=True,fletcher32=True,chunksizes=somlist，启动压缩，默认的zlib压缩级别是4
newncspacecraft_num = newnc.createVariable(spacecraft_num.name, spacecraft_num.dtype, spacecraft_num.dimensions,                                         fill_value=spacecraft_num_fill_value, zlib=True, fletcher32=True,
                                           chunksizes=[np.uint32(1122508), ])
# 设定newncspacecraft_num的变量属性
newncspacecraft_num.long_name = "BF-1小卫星序号"
newncspacecraft_num.coordinates = spacecraft_num.getncattr("coordinates")
newncspacecraft_num.units = spacecraft_num.getncattr("units")
newncspacecraft_num.valid_range = spacecraft_num.getncattr("valid_range")
newncspacecraft_num.comment = "飞行器编号\n\t 1= BF-1A, 2 = BF-1B, 97 = 仿真器, 98 = 其他数据"
# copy spacecraft_num to newncspacecraft_num
newncspacecraft_num[:] = spacecraft_num_data
```
**关闭冰箱**

```
# must close
nc.close()
newnc.close()
```
思路来源：
[NetCDF writing example](http://nbviewer.jupyter.org/github/Unidata/netcdf4-python/blob/master/examples/writing_netCDF.ipynb)


[1]: http://math.stackexchange.com/
[2]: https://github.com/jmcmanus/pagedown-extra "Pagedown Extra"
[3]: http://meta.math.stackexchange.com/questions/5020/mathjax-basic-tutorial-and-quick-reference
[4]: http://bramp.github.io/js-sequence-diagrams/
[5]: http://adrai.github.io/flowchart.js/
[6]: https://github.com/benweet/stackedit