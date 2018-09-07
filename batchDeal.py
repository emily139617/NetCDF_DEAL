# -*- coding: utf-8 -*-"""
"""
@Time:2018/9/6 17:03
@Author:qinjing
"""
import os
import re
import ParseNetCDF


def list_all_files(rootdir):
    import os
    _files = []
    list = os.listdir(rootdir)  # 列出文件夹下所有的目录与文件
    for i in range(0, len(list)):
        path = os.path.join(rootdir, list[i])
        if os.path.isdir(path):
            _files.extend(list_all_files(path))
        if os.path.isfile(path):
            _files.append(path)
    return _files


if __name__ == '__main__':
    listfiles = list_all_files("E:\\l2")
    listfilter = filter(lambda x: re.compile(r'.nc$').search(x), listfiles)
    # print len(listfilter)
    # print listfilter
    for i in listfilter:
        outstr = "bf1ab." + i[-56:-23] + ".l2." + "wind.v10.nc"
        outputfilename = "E:\\l2\\" + outstr;
        # ParseNetCDF.copyPartNCToLocalFileWithOutCom(i, outputfilename)
        ParseNetCDF.copyPartNCToLocalFile(i, outputfilename)
        # print i
        # print outstr
        # print i[-28:-20]
