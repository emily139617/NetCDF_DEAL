# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 09:56:11 2018

@author: Administrator
"""
# import netCDF4
import os

from netCDF4 import *
import numpy as np


# 解析函数
def parse_nc(filename):
    nc = Dataset(filename)
    #    print nc's varibles' names
    # for var in nc.variables.keys():
    #     print var
    # get dimensions for copy
    # print nc.dimensions
    # sample_dim=nc.dimensions["sample"]
    # print len(sample_dim)
    # print sample_dim.name
    # print "sample" in nc.dimensions
    # get particular variable like sp_lat to copy
    ddm_source = nc.variables["ddm_source"]
    spacecraft_num = nc.variables["spacecraft_num"]
    prn_code = nc.variables["prn_code"]
    sv_num = nc.variables["sv_num"]
    sc_lat = nc.variables["sc_lat"]
    sc_lon = nc.variables["sc_lon"]
    sc_alt = nc.variables["sc_alt"]
    fresnel_coeff = nc.variables["fresnel_coeff"]
    sample = nc.variables["sample"]
    # //fill_value需要自己制定，不能拷贝
    # print spacecraft_num.name, spacecraft_num.dimensions, spacecraft_num.dtype, spacecraft_num.shape, spacecraft_num.ncattrs(), spacecraft_num.getncattr(
    #     "_FillValue"),getlibversion()
    # print type(splatdata)
    # print splatdata
    # # <netCDF4._netCDF4.Variable>
    # print type(splat)
    print ddm_source[:]

def copyPartNCToLocalFile(inputFileName, outputFileName):
    # //get NC from local inputFile
    nc = Dataset(inputFileName)
    # 新建一个本地文件，格式是netcdf4_classic
    newnc = Dataset(outputFileName, "w", format="NETCDF4_CLASSIC")
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
    # get needed variables and copy to newnc
    ncvariables = nc.variables
    # !ddm_source!
    ddm_source = ncvariables["ddm_source"]
    ddm_source_data = ddm_source[:]
    # copy ddm_source to newnc
    ddm_source_fill_value = ddm_source.getncattr("_FillValue")
    newncddm_source = newnc.createVariable(ddm_source.name, ddm_source.dtype, ddm_source.dimensions,
                                           fill_value=ddm_source_fill_value)
    # 设定newncddm_source的变量属性
    newncddm_source.long_name = "L0数据源"
    newncddm_source.units = ddm_source.getncattr("units")
    newncddm_source.flag_values = ddm_source.getncattr("flag_values")
    newncddm_source.valid_range = ddm_source.getncattr("valid_range")
    newncddm_source.flag_meanings = "Simulator BF1 unknown"
    newncddm_source.comment = "用于风速反演的L0数据源\n\t 0 = 仿真数据 \n\t 1 = BF-1观测数据 \n\t 2 = 其他数据源\n"
    # copy ddm_source to newncddm_source
    newncddm_source[:] = ddm_source_data
    # !spacecraft_num!
    spacecraft_num = ncvariables["spacecraft_num"]
    spacecraft_num_data = spacecraft_num[:]
    # copy spacecraft_num to newnc
    spacecraft_num_fill_value = spacecraft_num.getncattr("_FillValue")
    newncspacecraft_num = newnc.createVariable(spacecraft_num.name, spacecraft_num.dtype, spacecraft_num.dimensions,
                                               fill_value=spacecraft_num_fill_value, zlib=True, fletcher32=True,
                                               chunksizes=[np.uint32(sample_dim.size), ])
    # 设定newncspacecraft_num的变量属性
    newncspacecraft_num.long_name = "BF-1小卫星序号"
    newncspacecraft_num.coordinates = spacecraft_num.getncattr("coordinates")
    newncspacecraft_num.units = spacecraft_num.getncattr("units")
    newncspacecraft_num.valid_range = spacecraft_num.getncattr("valid_range")
    newncspacecraft_num.comment = "飞行器编号\n\t 1= BF-1A, 2 = BF-1B, 97 = 仿真器, 98 = 其他数据"
    # copy spacecraft_num to newncspacecraft_num
    newncspacecraft_num[:] = spacecraft_num_data
    # !sample!
    sample = ncvariables["sample"]
    sample_data = sample[:]
    # copy sample to newnc
    newncsample = newnc.createVariable(sample.name, sample.dtype, sample.dimensions, zlib=True, fletcher32=True,
                                       chunksizes=[np.uint32(sample_dim.size/2), ])
    # 设定newncsample的变量属性
    newncsample.long_name = "采用序号"
    newncsample.units = sample.getncattr("units")
    newncsample.comment = "netCDF文件对应采样序号，从0开始"
    # copy sample to newncsample
    newncsample[:] = sample_data

    # !lat!
    lat = ncvariables["lat"]
    lat_data = lat[:]
    lat_fill_value = lat.getncattr("_FillValue")
    # copy lat to newnc
    newnclat = newnc.createVariable(lat.name, lat.dtype, lat.dimensions, zlib=True, fletcher32=True,
                                    chunksizes=[np.uint32(sample_dim.size/2), ], fill_value=lat_fill_value)
    # 设定newnclat的变量属性
    newnclat.long_name = "纬度"
    newnclat.standard_name = "latitude"
    newnclat.units = lat.getncattr("units")
    newnclat.valid_range = lat.getncattr("valid_range")
    newnclat.comment = "反演风速所用DDM数据纬度均值，北方为正"
    newnclat._CoordinateAxisType = "Lat"
    # copy lat to newnclat
    newnclat[:] = lat_data

    # !lon!
    lon = ncvariables["lon"]
    lon_data = lon[:]
    lon_fill_value = lon.getncattr("_FillValue")
    # copy lon to newnc
    newnclon = newnc.createVariable(lon.name, lon.dtype, lon.dimensions, zlib=True, fletcher32=True,
                                    chunksizes=[np.uint32(sample_dim.size/2), ], fill_value=lon_fill_value)
    # 设定newnclon的变量属性
    newnclon.long_name = "经度"
    newnclon.standard_name = "longtitude"
    newnclon.units = lon.getncattr("units")
    newnclon.valid_range = lon.getncattr("valid_range")
    newnclon.comment = "反演风速所用DDM数据经度均值，东方为正"
    newnclon._CoordinateAxisType = "lon"
    # copy lon to newnclon
    newnclon[:] = lon_data

    # !prn_code!
    prn_code = ncvariables["prn_code"]
    prn_code_data = prn_code[:]
    prn_code_fill_value = prn_code.getncattr("_FillValue")
    # copy prn_code to newnc
    newncprn_code = newnc.createVariable(prn_code.name, prn_code.dtype, prn_code.dimensions, zlib=True, fletcher32=True,
                                         chunksizes=[np.uint32(sample_dim.size), ], fill_value=prn_code_fill_value)
    # 设定newncprn_code的变量属性
    newncprn_code.long_name = "GPS PRN码"
    newncprn_code.coordinates = "sample_time lat lon"
    newncprn_code.units = prn_code.getncattr("units")
    newncprn_code.valid_range = prn_code.getncattr("valid_range")
    newncprn_code.comment = "用于风速反演的GPS PRN编码编号。取值0 到 32. \n\t 0  =  空闲.  1 – 32代表PRN码"
    # copy prn_code to newncprn_code
    newncprn_code[:] = prn_code_data

    # !sv_num!
    sv_num = ncvariables["sv_num"]
    sv_num_data = sv_num[:]
    sv_num_fill_value = sv_num.getncattr("_FillValue")
    # copy sv_num to newnc
    newncsv_num = newnc.createVariable(sv_num.name, sv_num.dtype, sv_num.dimensions, zlib=True, fletcher32=True,
                                       chunksizes=[np.uint32(sample_dim.size), ], fill_value=sv_num_fill_value)
    # 设定newncsv_num的变量属性
    newncsv_num.long_name = "GPS卫星编号"
    newncsv_num.coordinates = "sample_time lat lon"
    newncsv_num.units = sv_num.getncattr("units")
    newncsv_num.valid_range = sv_num.getncattr("valid_range")
    newncsv_num.comment = "发射对应prn_code码的GPS卫星编号"
    # copy sv_num to newncsv_num
    newncsv_num[:] = sv_num_data

    # ! antenna!
    antenna = ncvariables["antenna"]
    antenna_data = antenna[:]
    antenna_fill_value = antenna.getncattr("_FillValue")
    # copy antenna to newnc
    newncantenna = newnc.createVariable(antenna.name, antenna.dtype, antenna.dimensions, zlib=True, fletcher32=True,
                                        chunksizes=[np.uint32(sample_dim.size), ], fill_value=antenna_fill_value)
    # 设定newncantenna的变量属性
    newncantenna.long_name = "接收天线"
    newncantenna.coordinates = "sample_time lat lon"
    newncantenna.units = antenna.getncattr("units")
    newncantenna.valid_range = antenna.getncattr("valid_range")
    newncantenna.flag_values = antenna.getncattr("flag_values")
    newncantenna.flag_meanings = antenna.getncattr("flag_meanings")
    newncantenna.comment = "风速反演所用DDM信号对应的接收天线.\n\t0 = none\n\t1 = zenith (never used)\n\t2 = nadir_starboard\n\t3 = nadir_port"

    # copy antenna to newncantenna
    newncantenna[:] = antenna_data

    # !sc_lat!
    sc_lat = ncvariables["sc_lat"]
    sc_lat_data = sc_lat[:]
    sc_lat_fill_value = sc_lat.getncattr("_FillValue")
    # copy sc_lat to newnc
    newncsc_lat = newnc.createVariable(sc_lat.name, sc_lat.dtype, sc_lat.dimensions, zlib=True, fletcher32=True,
                                       chunksizes=[np.uint32(sample_dim.size/2), ], fill_value=sc_lat_fill_value)
    # 设定newncsc_lat的变量属性
    newncsc_lat.long_name = "小卫星指向，纬度"
    newncsc_lat.standard_name = "latitude"
    newncsc_lat.coordinates = "sample_time lat lon"
    newncsc_lat.units = sc_lat.getncattr("units")
    newncsc_lat.valid_range = sc_lat.getncattr("valid_range")
    newncsc_lat.comment = "用于风速反演的DDM数据对应小卫星纬度指向均值，北向为正"

    # copy sc_lat to newncsc_lat
    newncsc_lat[:] = sc_lat_data

    # !sc_lon!
    sc_lon = ncvariables["sc_lon"]
    sc_lon_data = sc_lon[:]
    sc_lon_fill_value = sc_lon.getncattr("_FillValue")
    # copy sc_lon to newnc
    newncsc_lon = newnc.createVariable(sc_lon.name, sc_lon.dtype, sc_lon.dimensions, zlib=True, fletcher32=True,
                                       chunksizes=[np.uint32(sample_dim.size/2), ], fill_value=sc_lon_fill_value)
    # 设定newncsc_lon的变量属性
    newncsc_lon.long_name = "小卫星指向：经度"
    newncsc_lon.standard_name = "longitude"
    newncsc_lon.coordinates = "sample_time lat lon"
    newncsc_lon.units = sc_lon.getncattr("units")
    newncsc_lon.valid_range = sc_lon.getncattr("valid_range")
    newncsc_lon.comment = "用于风速反演的DDM数据对应小卫星经度指向均值，东向为正"

    # copy sc_lon to newncsc_lon
    newncsc_lon[:] = sc_lon_data

    # !sc_alt!
    sc_alt = ncvariables["sc_alt"]
    sc_alt_data = sc_alt[:]
    sc_alt_fill_value = sc_alt.getncattr("_FillValue")
    # copy sc_alt to newnc
    newncsc_alt = newnc.createVariable(sc_alt.name, sc_alt.dtype, sc_alt.dimensions, zlib=True, fletcher32=True,
                                       chunksizes=[np.uint32(sample_dim.size/2), ], fill_value=sc_alt_fill_value)
    # 设定newncsc_alt的变量属性
    newncsc_alt.long_name = "小卫星高度"
    newncsc_alt.coordinates = "sample_time lat lon"
    newncsc_alt.units = sc_alt.getncattr("units")
    newncsc_alt.valid_range = sc_alt.getncattr("valid_range")
    newncsc_alt.comment = "用于风速反演的DDM数据对应小卫星在WGS-84中的高度均值"

    # copy sc_alt to newncsc_alt
    newncsc_alt[:] = sc_alt_data

    # !wind_speed!
    wind_speed = ncvariables["wind_speed"]
    wind_speed_data = wind_speed[:]
    wind_speed_fill_value = wind_speed.getncattr("_FillValue")
    # copy wind_speed to newnc
    newncwind_speed = newnc.createVariable(wind_speed.name, wind_speed.dtype, wind_speed.dimensions, zlib=True,
                                           fletcher32=True,
                                           chunksizes=[np.uint32(sample_dim.size/2), ], fill_value=wind_speed_fill_value)
    # 设定newncwind_speed的变量属性
    newncwind_speed.long_name = "NBRCS和LES风速反演结果按最小方差评估得到的风速"
    newncwind_speed.standard_name = "wind_speed"
    newncwind_speed.coordinates = "sample_time lat lon"
    newncwind_speed.units = wind_speed.getncattr("units")
    newncwind_speed.valid_range = wind_speed.getncattr("valid_range")
    newncwind_speed.comment = "以lat和lon为中心的10米表面风速"

    # copy wind_speed to newncwind_speed
    newncwind_speed[:] = wind_speed_data

    # !fds_nbrcs_wind_speed!
    fds_nbrcs_wind_speed = ncvariables["fds_nbrcs_wind_speed"]
    fds_nbrcs_wind_speed_data = fds_nbrcs_wind_speed[:]
    fds_nbrcs_wind_speed_fill_value = fds_nbrcs_wind_speed.getncattr("_FillValue")
    # copy fds_nbrcs_wind_speed to newnc
    newncfds_nbrcs_wind_speed = newnc.createVariable(fds_nbrcs_wind_speed.name, fds_nbrcs_wind_speed.dtype,
                                                     fds_nbrcs_wind_speed.dimensions, zlib=True,
                                                     fletcher32=True,
                                                     chunksizes=[np.uint32(sample_dim.size/2), ],
                                                     fill_value=fds_nbrcs_wind_speed_fill_value)
    # 设定newncfds_nbrcs_wind_speed的变量属性
    newncfds_nbrcs_wind_speed.long_name = "由nbrcs 反演的风速"
    newncfds_nbrcs_wind_speed.standard_name = "wind_speed"
    newncfds_nbrcs_wind_speed.coordinates = "sample_time lat lon"
    newncfds_nbrcs_wind_speed.units = fds_nbrcs_wind_speed.getncattr("units")
    newncfds_nbrcs_wind_speed.valid_range = fds_nbrcs_wind_speed.getncattr("valid_range")
    newncfds_nbrcs_wind_speed.comment = "由NBRCS应用FDS反演获得的10米表面风速"

    # copy fds_nbrcs_wind_speed to newncfds_nbrcs_wind_speed
    newncfds_nbrcs_wind_speed[:] = fds_nbrcs_wind_speed_data

    # !fds_les_wind_speed!
    fds_les_wind_speed = ncvariables["fds_les_wind_speed"]
    fds_les_wind_speed_data = fds_les_wind_speed[:]
    fds_les_wind_speed_fill_value = fds_les_wind_speed.getncattr("_FillValue")
    # copy fds_les_wind_speed to newnc
    newncfds_les_wind_speed = newnc.createVariable(fds_les_wind_speed.name, fds_les_wind_speed.dtype,
                                                   fds_les_wind_speed.dimensions, zlib=True,
                                                   fletcher32=True,
                                                   chunksizes=[np.uint32(sample_dim.size/2), ],
                                                   fill_value=fds_les_wind_speed_fill_value)
    # 设定newncfds_les_wind_speed的变量属性
    newncfds_les_wind_speed.long_name = "由les反演的风速"
    newncfds_les_wind_speed.standard_name = "wind_speed"
    newncfds_les_wind_speed.coordinates = "sample_time lat lon"
    newncfds_les_wind_speed.units = fds_les_wind_speed.getncattr("units")
    newncfds_les_wind_speed.valid_range = fds_les_wind_speed.getncattr("valid_range")
    newncfds_les_wind_speed.comment = "由LES应用FDS反演获得的10米表面风速"

    # copy fds_les_wind_speed to newncfds_les_wind_speed
    newncfds_les_wind_speed[:] = fds_les_wind_speed_data

    # !wind_speed_uncertainty!
    wind_speed_uncertainty = ncvariables["wind_speed_uncertainty"]
    wind_speed_uncertainty_data = wind_speed_uncertainty[:]
    wind_speed_uncertainty_fill_value = wind_speed_uncertainty.getncattr("_FillValue")
    # copy wind_speed_uncertainty to newnc
    newncwind_speed_uncertainty = newnc.createVariable(wind_speed_uncertainty.name, wind_speed_uncertainty.dtype,
                                                       wind_speed_uncertainty.dimensions, zlib=True,
                                                       fletcher32=True,
                                                       chunksizes=[np.uint32(sample_dim.size/2), ],
                                                       fill_value=wind_speed_uncertainty_fill_value)
    # 设定newncwind_speed_uncertainty的变量属性
    newncwind_speed_uncertainty.long_name = "风速反演不确定度"
    newncwind_speed_uncertainty.coordinates = "sample_time lat lon"
    newncwind_speed_uncertainty.units = wind_speed_uncertainty.getncattr("units")
    newncwind_speed_uncertainty.valid_range = wind_speed_uncertainty.getncattr("valid_range")
    newncwind_speed_uncertainty.comment = "依赖于RCG的风速反演标准偏差"

    # copy wind_speed_uncertainty to newncwind_speed_uncertainty
    newncwind_speed_uncertainty[:] = wind_speed_uncertainty_data

    # !incidence_angle!
    incidence_angle = ncvariables["incidence_angle"]
    incidence_angle_data = incidence_angle[:]
    incidence_angle_fill_value = incidence_angle.getncattr("_FillValue")
    # copy incidence_angle to newnc
    newncincidence_angle = newnc.createVariable(incidence_angle.name, incidence_angle.dtype,
                                                incidence_angle.dimensions, zlib=True,
                                                fletcher32=True,
                                                chunksizes=[np.uint32(sample_dim.size/2), ],
                                                fill_value=incidence_angle_fill_value)
    # 设定newncincidence_angle的变量属性
    newncincidence_angle.long_name = "镜面点入射角"
    newncincidence_angle.coordinates = "sample_time lat lon"
    newncincidence_angle.units = incidence_angle.getncattr("units")
    newncincidence_angle.valid_range = incidence_angle.getncattr("valid_range")
    newncincidence_angle.comment = "风速反演所用DDMs数据镜面点入射角均值"

    # copy incidence_angle to newncincidence_angle
    newncincidence_angle[:] = incidence_angle_data

    # !range_corr_gain!
    range_corr_gain = ncvariables["range_corr_gain"]
    range_corr_gain_data = range_corr_gain[:]
    range_corr_gain_fill_value = range_corr_gain.getncattr("_FillValue")
    # copy range_corr_gain to newnc
    newncrange_corr_gain = newnc.createVariable(range_corr_gain.name, range_corr_gain.dtype,
                                                range_corr_gain.dimensions, zlib=True,
                                                fletcher32=True,
                                                chunksizes=[np.uint32(sample_dim.size/2), ],
                                                fill_value=range_corr_gain_fill_value)
    # 设定newncrange_corr_gain的变量属性
    newncrange_corr_gain.long_name = "距离校正增益"
    newncrange_corr_gain.coordinates = "sample_time lat lon"
    newncrange_corr_gain.units = range_corr_gain.getncattr("units")
    newncrange_corr_gain.valid_range = range_corr_gain.getncattr("valid_range")
    newncrange_corr_gain.comment = "用于风速反演的RCG值"

    # copy range_corr_gain to newncrange_corr_gain
    newncrange_corr_gain[:] = range_corr_gain_data

    # !fresnel_coeff!
    fresnel_coeff = ncvariables["fresnel_coeff"]
    fresnel_coeff_data = fresnel_coeff[:]
    fresnel_coeff_fill_value = fresnel_coeff.getncattr("_FillValue")
    # copy fresnel_coeff to newnc
    newncfresnel_coeff = newnc.createVariable(fresnel_coeff.name, fresnel_coeff.dtype,
                                              fresnel_coeff.dimensions, zlib=True,
                                              fletcher32=True,
                                              chunksizes=[np.uint32(sample_dim.size/2), ],
                                              fill_value=fresnel_coeff_fill_value)
    # 设定newncfresnel_coeff的变量属性
    newncfresnel_coeff.long_name = "菲涅尔功率反射系数"
    newncfresnel_coeff.coordinates = "sample_time lat lon"
    newncfresnel_coeff.units = fresnel_coeff.getncattr("units")
    newncfresnel_coeff.valid_range = fresnel_coeff.getncattr("valid_range")
    newncfresnel_coeff.comment = "lat, lon对应位置平滑洋面左旋圆偏振菲涅尔电场电压反射系数的平方（@1575MHz）"

    # copy fresnel_coeff to newncfresnel_coeff
    newncfresnel_coeff[:] = fresnel_coeff_data

    # !sample_time!
    sample_time = ncvariables["sample_time"]
    sample_time_data = sample_time[:]
    # copy sample_time to newnc
    newncsample_time = newnc.createVariable(sample_time.name, sample_time.dtype,
                                            sample_time.dimensions, zlib=True,
                                            fletcher32=True,
                                            chunksizes=[np.uint32(sample_dim.size/2), ])    # 设定newncsample_time的变量属性
    newncsample_time.long_name = "采样时刻"
    newncsample_time.standard_name = "time"
    newncsample_time.calendar = "gregorian"
    newncsample_time.coordinates = "sample_time lat lon"
    newncsample_time.units = sample_time.getncattr("units")
    newncsample_time.comment = "反演风速所用DDM对应的ddm_timestamp_utc时间，相对于time_coverage_start"
    newncsample_time._CoordinateAxisType = "Time"
    # copy sample_time to newncsample_time
    newncsample_time[:] = sample_time_data
    # copy nc global attrs to new nc
    newnc.project = "BF-1"
    newnc.featureType = "trajectory"
    newnc.summary = "TEMP NONE"
    newnc.processing_level = "2"
    newnc.comment = "TEMP NONE"
    newnc.creator_type = "institution"
    newnc.creator_name = "Prototype Data Center"
    newnc.sensor = "504 DDMI"
    newnc.version_id = "1.0a"
    newnc.title = "BF-1 Level 2 Prototype Data Center Version 1.0a"
    newnc.ShortName = "BF1_L2_V1.0a"
    newnc.netcdf_version_id = getlibversion()
    newnc.date_created = nc.getncattr("date_created")
    newnc.date_issued = nc.getncattr("date_issued")
    newnc.source = "TEMP NONE"
    newnc.time_coverage_resolution = "P0DT0H0M1S"
    newnc.time_coverage_start = nc.getncattr("time_coverage_start")
    newnc.time_coverage_end = nc.getncattr("time_coverage_end")
    newnc.time_coverage_duration = "P1DT0H0M0S"
    newnc.l2_algorithm_version = "1.0a"
    newnc.time_averaging_lookup_tables_version = "1"
    newnc.nbrcs_wind_lookup_tables_version = "1"
    newnc.les_wind_lookup_tables_version = "1"
    newnc.covariance_lookup_tables_version = "1"
    newnc.standard_deviation_lookup_tables_version = "1"
    newnc.geospatial_lat_min = nc.getncattr("geospatial_lat_min")
    newnc.geospatial_lat_max = nc.getncattr("geospatial_lat_max")
    newnc.geospatial_lon_min = nc.getncattr("geospatial_lon_min")
    newnc.geospatial_lon_max = nc.getncattr("geospatial_lon_max")
    newnc.platform = "BF-1 A/B"
    newnc._CoordSysBuilder = "ucar.nc2.dataset.conv.CF1Convention"

    # must close
    nc.close()
    newnc.close()
    print "文件已经读取完毕，且关闭"

def copyPartNCToLocalFileWithOutCom(inputFileName, outputFileName):
    # //get NC from local inputFile
    nc = Dataset(inputFileName)
    # 新建一个本地文件，格式是netcdf4_classic
    newnc = Dataset(outputFileName, "w", format="NETCDF4_CLASSIC")
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
    # get needed variables and copy to newnc
    ncvariables = nc.variables
    # !ddm_source!
    ddm_source = ncvariables["ddm_source"]
    ddm_source_data = ddm_source[:]
    # copy ddm_source to newnc
    ddm_source_fill_value = ddm_source.getncattr("_FillValue")
    newncddm_source = newnc.createVariable(ddm_source.name, ddm_source.dtype, ddm_source.dimensions,
                                           fill_value=ddm_source_fill_value)
    # 设定newncddm_source的变量属性
    newncddm_source.long_name = "L0数据源"
    newncddm_source.units = ddm_source.getncattr("units")
    newncddm_source.flag_values = ddm_source.getncattr("flag_values")
    newncddm_source.valid_range = ddm_source.getncattr("valid_range")
    newncddm_source.flag_meanings = "Simulator BF1 unknown"
    newncddm_source.comment = "用于风速反演的L0数据源\n\t 0 = 仿真数据 \n\t 1 = BF-1观测数据 \n\t 2 = 其他数据源\n"
    # copy ddm_source to newncddm_source
    newncddm_source[:] = ddm_source_data
    # !spacecraft_num!
    spacecraft_num = ncvariables["spacecraft_num"]
    spacecraft_num_data = spacecraft_num[:]
    # copy spacecraft_num to newnc
    spacecraft_num_fill_value = spacecraft_num.getncattr("_FillValue")
    newncspacecraft_num = newnc.createVariable(spacecraft_num.name, spacecraft_num.dtype, spacecraft_num.dimensions,
                                               fill_value=spacecraft_num_fill_value)
    # 设定newncspacecraft_num的变量属性
    newncspacecraft_num.long_name = "BF-1小卫星序号"
    newncspacecraft_num.coordinates = spacecraft_num.getncattr("coordinates")
    newncspacecraft_num.units = spacecraft_num.getncattr("units")
    newncspacecraft_num.valid_range = spacecraft_num.getncattr("valid_range")
    newncspacecraft_num.comment = "飞行器编号\n\t 1= BF-1A, 2 = BF-1B, 97 = 仿真器, 98 = 其他数据"
    # copy spacecraft_num to newncspacecraft_num
    newncspacecraft_num[:] = spacecraft_num_data
    # !sample!
    sample = ncvariables["sample"]
    sample_data = sample[:]
    # copy sample to newnc
    newncsample = newnc.createVariable(sample.name, sample.dtype, sample.dimensions)
    # 设定newncsample的变量属性
    newncsample.long_name = "采用序号"
    newncsample.units = sample.getncattr("units")
    newncsample.comment = "netCDF文件对应采样序号，从0开始"
    # copy sample to newncsample
    newncsample[:] = sample_data

    # !lat!
    lat = ncvariables["lat"]
    lat_data = lat[:]
    lat_fill_value = lat.getncattr("_FillValue")
    # copy lat to newnc
    newnclat = newnc.createVariable(lat.name, lat.dtype, lat.dimensions,fill_value=lat_fill_value)
    # 设定newnclat的变量属性
    newnclat.long_name = "纬度"
    newnclat.standard_name = "latitude"
    newnclat.units = lat.getncattr("units")
    newnclat.valid_range = lat.getncattr("valid_range")
    newnclat.comment = "反演风速所用DDM数据纬度均值，北方为正"
    newnclat._CoordinateAxisType = "Lat"
    # copy lat to newnclat
    newnclat[:] = lat_data

    # !lon!
    lon = ncvariables["lon"]
    lon_data = lon[:]
    lon_fill_value = lon.getncattr("_FillValue")
    # copy lon to newnc
    newnclon = newnc.createVariable(lon.name, lon.dtype, lon.dimensions,fill_value=lon_fill_value)
    # 设定newnclon的变量属性
    newnclon.long_name = "经度"
    newnclon.standard_name = "longtitude"
    newnclon.units = lon.getncattr("units")
    newnclon.valid_range = lon.getncattr("valid_range")
    newnclon.comment = "反演风速所用DDM数据经度均值，东方为正"
    newnclon._CoordinateAxisType = "lon"
    # copy lon to newnclon
    newnclon[:] = lon_data

    # !prn_code!
    prn_code = ncvariables["prn_code"]
    prn_code_data = prn_code[:]
    prn_code_fill_value = prn_code.getncattr("_FillValue")
    # copy prn_code to newnc
    newncprn_code = newnc.createVariable(prn_code.name, prn_code.dtype, prn_code.dimensions, fill_value=prn_code_fill_value)
    # 设定newncprn_code的变量属性
    newncprn_code.long_name = "GPS PRN码"
    newncprn_code.coordinates = "sample_time lat lon"
    newncprn_code.units = prn_code.getncattr("units")
    newncprn_code.valid_range = prn_code.getncattr("valid_range")
    newncprn_code.comment = "用于风速反演的GPS PRN编码编号。取值0 到 32. \n\t 0  =  空闲.  1 – 32代表PRN码"
    # copy prn_code to newncprn_code
    newncprn_code[:] = prn_code_data

    # !sv_num!
    sv_num = ncvariables["sv_num"]
    sv_num_data = sv_num[:]
    sv_num_fill_value = sv_num.getncattr("_FillValue")
    # copy sv_num to newnc
    newncsv_num = newnc.createVariable(sv_num.name, sv_num.dtype, sv_num.dimensions,fill_value=sv_num_fill_value)
    # 设定newncsv_num的变量属性
    newncsv_num.long_name = "GPS卫星编号"
    newncsv_num.coordinates = "sample_time lat lon"
    newncsv_num.units = sv_num.getncattr("units")
    newncsv_num.valid_range = sv_num.getncattr("valid_range")
    newncsv_num.comment = "发射对应prn_code码的GPS卫星编号"
    # copy sv_num to newncsv_num
    newncsv_num[:] = sv_num_data

    # ! antenna!
    antenna = ncvariables["antenna"]
    antenna_data = antenna[:]
    antenna_fill_value = antenna.getncattr("_FillValue")
    # copy antenna to newnc
    newncantenna = newnc.createVariable(antenna.name, antenna.dtype, antenna.dimensions,fill_value=antenna_fill_value)
    # 设定newncantenna的变量属性
    newncantenna.long_name = "接收天线"
    newncantenna.coordinates = "sample_time lat lon"
    newncantenna.units = antenna.getncattr("units")
    newncantenna.valid_range = antenna.getncattr("valid_range")
    newncantenna.flag_values = antenna.getncattr("flag_values")
    newncantenna.flag_meanings = antenna.getncattr("flag_meanings")
    newncantenna.comment = "风速反演所用DDM信号对应的接收天线.\n\t0 = none\n\t1 = zenith (never used)\n\t2 = nadir_starboard\n\t3 = nadir_port"

    # copy antenna to newncantenna
    newncantenna[:] = antenna_data

    # !sc_lat!
    sc_lat = ncvariables["sc_lat"]
    sc_lat_data = sc_lat[:]
    sc_lat_fill_value = sc_lat.getncattr("_FillValue")
    # copy sc_lat to newnc
    newncsc_lat = newnc.createVariable(sc_lat.name, sc_lat.dtype, sc_lat.dimensions,fill_value=sc_lat_fill_value)
    # 设定newncsc_lat的变量属性
    newncsc_lat.long_name = "小卫星指向，纬度"
    newncsc_lat.standard_name = "latitude"
    newncsc_lat.coordinates = "sample_time lat lon"
    newncsc_lat.units = sc_lat.getncattr("units")
    newncsc_lat.valid_range = sc_lat.getncattr("valid_range")
    newncsc_lat.comment = "用于风速反演的DDM数据对应小卫星纬度指向均值，北向为正"

    # copy sc_lat to newncsc_lat
    newncsc_lat[:] = sc_lat_data

    # !sc_lon!
    sc_lon = ncvariables["sc_lon"]
    sc_lon_data = sc_lon[:]
    sc_lon_fill_value = sc_lon.getncattr("_FillValue")
    # copy sc_lon to newnc
    newncsc_lon = newnc.createVariable(sc_lon.name, sc_lon.dtype, sc_lon.dimensions,fill_value=sc_lon_fill_value)
    # 设定newncsc_lon的变量属性
    newncsc_lon.long_name = "小卫星指向：经度"
    newncsc_lon.standard_name = "longitude"
    newncsc_lon.coordinates = "sample_time lat lon"
    newncsc_lon.units = sc_lon.getncattr("units")
    newncsc_lon.valid_range = sc_lon.getncattr("valid_range")
    newncsc_lon.comment = "用于风速反演的DDM数据对应小卫星经度指向均值，东向为正"

    # copy sc_lon to newncsc_lon
    newncsc_lon[:] = sc_lon_data

    # !sc_alt!
    sc_alt = ncvariables["sc_alt"]
    sc_alt_data = sc_alt[:]
    sc_alt_fill_value = sc_alt.getncattr("_FillValue")
    # copy sc_alt to newnc
    newncsc_alt = newnc.createVariable(sc_alt.name, sc_alt.dtype, sc_alt.dimensions,fill_value=sc_alt_fill_value)
    # 设定newncsc_alt的变量属性
    newncsc_alt.long_name = "小卫星高度"
    newncsc_alt.coordinates = "sample_time lat lon"
    newncsc_alt.units = sc_alt.getncattr("units")
    newncsc_alt.valid_range = sc_alt.getncattr("valid_range")
    newncsc_alt.comment = "用于风速反演的DDM数据对应小卫星在WGS-84中的高度均值"

    # copy sc_alt to newncsc_alt
    newncsc_alt[:] = sc_alt_data

    # !wind_speed!
    wind_speed = ncvariables["wind_speed"]
    wind_speed_data = wind_speed[:]
    wind_speed_fill_value = wind_speed.getncattr("_FillValue")
    # copy wind_speed to newnc
    newncwind_speed = newnc.createVariable(wind_speed.name, wind_speed.dtype, wind_speed.dimensions, fill_value=wind_speed_fill_value)
    # 设定newncwind_speed的变量属性
    newncwind_speed.long_name = "NBRCS和LES风速反演结果按最小方差评估得到的风速"
    newncwind_speed.standard_name = "wind_speed"
    newncwind_speed.coordinates = "sample_time lat lon"
    newncwind_speed.units = wind_speed.getncattr("units")
    newncwind_speed.valid_range = wind_speed.getncattr("valid_range")
    newncwind_speed.comment = "以lat和lon为中心的10米表面风速"

    # copy wind_speed to newncwind_speed
    newncwind_speed[:] = wind_speed_data

    # !fds_nbrcs_wind_speed!
    fds_nbrcs_wind_speed = ncvariables["fds_nbrcs_wind_speed"]
    fds_nbrcs_wind_speed_data = fds_nbrcs_wind_speed[:]
    fds_nbrcs_wind_speed_fill_value = fds_nbrcs_wind_speed.getncattr("_FillValue")
    # copy fds_nbrcs_wind_speed to newnc
    newncfds_nbrcs_wind_speed = newnc.createVariable(fds_nbrcs_wind_speed.name, fds_nbrcs_wind_speed.dtype,
                                                     fds_nbrcs_wind_speed.dimensions,
                                                     fill_value=fds_nbrcs_wind_speed_fill_value)
    # 设定newncfds_nbrcs_wind_speed的变量属性
    newncfds_nbrcs_wind_speed.long_name = "由nbrcs 反演的风速"
    newncfds_nbrcs_wind_speed.standard_name = "wind_speed"
    newncfds_nbrcs_wind_speed.coordinates = "sample_time lat lon"
    newncfds_nbrcs_wind_speed.units = fds_nbrcs_wind_speed.getncattr("units")
    newncfds_nbrcs_wind_speed.valid_range = fds_nbrcs_wind_speed.getncattr("valid_range")
    newncfds_nbrcs_wind_speed.comment = "由NBRCS应用FDS反演获得的10米表面风速"

    # copy fds_nbrcs_wind_speed to newncfds_nbrcs_wind_speed
    newncfds_nbrcs_wind_speed[:] = fds_nbrcs_wind_speed_data

    # !fds_les_wind_speed!
    fds_les_wind_speed = ncvariables["fds_les_wind_speed"]
    fds_les_wind_speed_data = fds_les_wind_speed[:]
    fds_les_wind_speed_fill_value = fds_les_wind_speed.getncattr("_FillValue")
    # copy fds_les_wind_speed to newnc
    newncfds_les_wind_speed = newnc.createVariable(fds_les_wind_speed.name, fds_les_wind_speed.dtype,
                                                   fds_les_wind_speed.dimensions,
                                                   fill_value=fds_les_wind_speed_fill_value)
    # 设定newncfds_les_wind_speed的变量属性
    newncfds_les_wind_speed.long_name = "由les反演的风速"
    newncfds_les_wind_speed.standard_name = "wind_speed"
    newncfds_les_wind_speed.coordinates = "sample_time lat lon"
    newncfds_les_wind_speed.units = fds_les_wind_speed.getncattr("units")
    newncfds_les_wind_speed.valid_range = fds_les_wind_speed.getncattr("valid_range")
    newncfds_les_wind_speed.comment = "由LES应用FDS反演获得的10米表面风速"

    # copy fds_les_wind_speed to newncfds_les_wind_speed
    newncfds_les_wind_speed[:] = fds_les_wind_speed_data

    # !wind_speed_uncertainty!
    wind_speed_uncertainty = ncvariables["wind_speed_uncertainty"]
    wind_speed_uncertainty_data = wind_speed_uncertainty[:]
    wind_speed_uncertainty_fill_value = wind_speed_uncertainty.getncattr("_FillValue")
    # copy wind_speed_uncertainty to newnc
    newncwind_speed_uncertainty = newnc.createVariable(wind_speed_uncertainty.name, wind_speed_uncertainty.dtype,
                                                       wind_speed_uncertainty.dimensions,
                                                       fill_value=wind_speed_uncertainty_fill_value)
    # 设定newncwind_speed_uncertainty的变量属性
    newncwind_speed_uncertainty.long_name = "风速反演不确定度"
    newncwind_speed_uncertainty.coordinates = "sample_time lat lon"
    newncwind_speed_uncertainty.units = wind_speed_uncertainty.getncattr("units")
    newncwind_speed_uncertainty.valid_range = wind_speed_uncertainty.getncattr("valid_range")
    newncwind_speed_uncertainty.comment = "依赖于RCG的风速反演标准偏差"

    # copy wind_speed_uncertainty to newncwind_speed_uncertainty
    newncwind_speed_uncertainty[:] = wind_speed_uncertainty_data

    # !incidence_angle!
    incidence_angle = ncvariables["incidence_angle"]
    incidence_angle_data = incidence_angle[:]
    incidence_angle_fill_value = incidence_angle.getncattr("_FillValue")
    # copy incidence_angle to newnc
    newncincidence_angle = newnc.createVariable(incidence_angle.name, incidence_angle.dtype,
                                                incidence_angle.dimensions,
                                                fill_value=incidence_angle_fill_value)
    # 设定newncincidence_angle的变量属性
    newncincidence_angle.long_name = "镜面点入射角"
    newncincidence_angle.coordinates = "sample_time lat lon"
    newncincidence_angle.units = incidence_angle.getncattr("units")
    newncincidence_angle.valid_range = incidence_angle.getncattr("valid_range")
    newncincidence_angle.comment = "风速反演所用DDMs数据镜面点入射角均值"

    # copy incidence_angle to newncincidence_angle
    newncincidence_angle[:] = incidence_angle_data

    # !range_corr_gain!
    range_corr_gain = ncvariables["range_corr_gain"]
    range_corr_gain_data = range_corr_gain[:]
    range_corr_gain_fill_value = range_corr_gain.getncattr("_FillValue")
    # copy range_corr_gain to newnc
    newncrange_corr_gain = newnc.createVariable(range_corr_gain.name, range_corr_gain.dtype,
                                                range_corr_gain.dimensions,
                                                fill_value=range_corr_gain_fill_value)
    # 设定newncrange_corr_gain的变量属性
    newncrange_corr_gain.long_name = "距离校正增益"
    newncrange_corr_gain.coordinates = "sample_time lat lon"
    newncrange_corr_gain.units = range_corr_gain.getncattr("units")
    newncrange_corr_gain.valid_range = range_corr_gain.getncattr("valid_range")
    newncrange_corr_gain.comment = "用于风速反演的RCG值"

    # copy range_corr_gain to newncrange_corr_gain
    newncrange_corr_gain[:] = range_corr_gain_data

    # !fresnel_coeff!
    fresnel_coeff = ncvariables["fresnel_coeff"]
    fresnel_coeff_data = fresnel_coeff[:]
    fresnel_coeff_fill_value = fresnel_coeff.getncattr("_FillValue")
    # copy fresnel_coeff to newnc
    newncfresnel_coeff = newnc.createVariable(fresnel_coeff.name, fresnel_coeff.dtype,
                                              fresnel_coeff.dimensions,
                                              fill_value=fresnel_coeff_fill_value)
    # 设定newncfresnel_coeff的变量属性
    newncfresnel_coeff.long_name = "菲涅尔功率反射系数"
    newncfresnel_coeff.coordinates = "sample_time lat lon"
    newncfresnel_coeff.units = fresnel_coeff.getncattr("units")
    newncfresnel_coeff.valid_range = fresnel_coeff.getncattr("valid_range")
    newncfresnel_coeff.comment = "lat, lon对应位置平滑洋面左旋圆偏振菲涅尔电场电压反射系数的平方（@1575MHz）"

    # copy fresnel_coeff to newncfresnel_coeff
    newncfresnel_coeff[:] = fresnel_coeff_data

    # !sample_time!
    sample_time = ncvariables["sample_time"]
    sample_time_data = sample_time[:]
    # copy sample_time to newnc
    newncsample_time = newnc.createVariable(sample_time.name, sample_time.dtype,
                                            sample_time.dimensions)    # 设定newncsample_time的变量属性
    newncsample_time.long_name = "采样时刻"
    newncsample_time.standard_name = "time"
    newncsample_time.calendar = "gregorian"
    newncsample_time.coordinates = "sample_time lat lon"
    newncsample_time.units = sample_time.getncattr("units")
    newncsample_time.comment = "反演风速所用DDM对应的ddm_timestamp_utc时间，相对于time_coverage_start"
    newncsample_time._CoordinateAxisType = "Time"
    # copy sample_time to newncsample_time
    newncsample_time[:] = sample_time_data
    # copy nc global attrs to new nc
    newnc.project = "BF-1"
    newnc.featureType = "trajectory"
    newnc.summary = "TEMP NONE"
    newnc.processing_level = "2"
    newnc.comment = "TEMP NONE"
    newnc.creator_type = "institution"
    newnc.creator_name = "Prototype Data Center"
    newnc.sensor = "504 DDMI"
    newnc.version_id = "1.0a"
    newnc.title = "BF-1 Level 2 Prototype Data Center Version 1.0a"
    newnc.ShortName = "BF1_L2_V1.0a"
    newnc.netcdf_version_id = getlibversion()
    newnc.date_created = nc.getncattr("date_created")
    newnc.date_issued = nc.getncattr("date_issued")
    newnc.source = "TEMP NONE"
    newnc.time_coverage_resolution = "P0DT0H0M1S"
    newnc.time_coverage_start = nc.getncattr("time_coverage_start")
    newnc.time_coverage_end = nc.getncattr("time_coverage_end")
    newnc.time_coverage_duration = "P1DT0H0M0S"
    newnc.l2_algorithm_version = "1.0a"
    newnc.time_averaging_lookup_tables_version = "1"
    newnc.nbrcs_wind_lookup_tables_version = "1"
    newnc.les_wind_lookup_tables_version = "1"
    newnc.covariance_lookup_tables_version = "1"
    newnc.standard_deviation_lookup_tables_version = "1"
    newnc.geospatial_lat_min = nc.getncattr("geospatial_lat_min")
    newnc.geospatial_lat_max = nc.getncattr("geospatial_lat_max")
    newnc.geospatial_lon_min = nc.getncattr("geospatial_lon_min")
    newnc.geospatial_lon_max = nc.getncattr("geospatial_lon_max")
    newnc.platform = "BF-1 A/B"
    newnc._CoordSysBuilder = "ucar.nc2.dataset.conv.CF1Convention"

    # must close
    nc.close()
    newnc.close()
    print "文件已经读取完毕，且关闭"

# 读取文件大小，做对比分析
def get_FileSize(filePath):
    filePath = unicode(filePath,"utf8")
    fsize = os.path.getsize(filePath)
    fsize = fsize / float(1024*1024)
    return round(fsize,2)

if __name__ == '__main__':
    # parse_nc("F:\\fsdata\\cygl2.nc")
    print get_FileSize("F:\\fsdata\\cyg03.nc")
    #17.8.1
    # copyPartNCToLocalFile("E:\\l2\\213-20170801\\cyg.ddmi.s20170801-000000-e20170801-235959.l2.wind-mss.a20.d20.nc", "E:\\l2\\bf1ab.s20170801-000000-e20170801-235959.l2.wind.v10.nc")
    # 17.8.17
    # copyPartNCToLocalFile("E:\\l2\\229\\cyg.ddmi.s20170817-000000-e20170817-235959.l2.wind-mss.a20.d20.nc","E:\\l2\\bf1ab.s20170817-000000-e20170817-235959.l2.wind.v10.nc")
    #17.8.23
    #copyPartNCToLocalFile("E:\\l2\\235\\cyg.ddmi.s20170823-000000-e20170823-235959.l2.wind-mss.a20.d20.nc","E:\\l2\\bf1ab.s20170823-000000-e20170823-235959.l2.wind.v10.nc")
    #17.9.1
    # copyPartNCToLocalFile("E:\\l2\\244-20170901\\cyg.ddmi.s20170901-000000-e20170901-235959.l2.wind-mss.a20.d20.nc", "E:\\l2\\bf1ab.s20170901-000000-e20170901-235959.l2.wind.v10.nc")
    #17.9.17
    #copyPartNCToLocalFile("E:\\l2\\260\\cyg.ddmi.s20170917-000000-e20170917-235959.l2.wind-mss.a20.d20.nc", "E:\\l2\\bf1ab.s20170917-000000-e20170917-235959.l2.wind.v10.nc")
    #17.9.20
    # copyPartNCToLocalFile("E:\\l2\\263\\cyg.ddmi.s20170920-000000-e20170920-235959.l2.wind-mss.a20.d20.nc", "E:\\l2\\bf1ab.s20170920-000000-e20170920-235959.l2.wind.v10.nc")
