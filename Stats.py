# -*- coding: utf-8 -*-"""
"""
@Time:2018/9/7 10:46
@Author:qinjing
"""
import os
import csv
import ParseNetCDF
import batchDeal

if __name__ == '__main__':
    timeList = []
    wysList = []
    ysList = []
    yslList = []
    listfiles = batchDeal.list_all_files("E:\\wxdatawys")
    for i in listfiles:
        timeList.append(i[20:28])
        wyssize = ParseNetCDF.get_FileSize(i)
        yssize = ParseNetCDF.get_FileSize(i.replace("wxdatawys", "wxdatays"))
        wysList.append(wyssize)
        ysList.append(yssize)
        yslList.append(round(yssize / wyssize, 2))
    # python2可以用file替代open
    with open("E:\\test.csv", "w") as csvfile:
        writer = csv.writer(csvfile)

        # 先写入columns_name
        writer.writerow(["时间", "未压缩文件大小", "压缩后文件大小", "压缩率"])
        # 写入多行用writerows
        for i in range(len(timeList)):
            writer.writerow([timeList[i],wysList[i],ysList[i],yslList[i]])
