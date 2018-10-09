# -*- coding: utf-8 -*-"""
"""
@Time:2018/9/17 9:00
@Author:qinjing
"""
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

if __name__ == '__main__':
    # 点图实现color渐变
    nc = Dataset("F:\\wxdata\\wxdatays\\bf1ab.s20170801-000000-e20170801-235959.l2.wind.v10.nc")
    lats = nc.variables["lat"][:]
    longs = nc.variables["lon"][:]
    longs = np.where(longs > 180, longs - 360, longs)
    # create Basemap instance for Robinson projection.
    map = Basemap()
    map.drawcoastlines()
    wind = nc.variables["wind_speed"]
    wind_speed = wind[:]
    max_wind_speed = max(wind_speed)
    min_wind_speed = min(wind_speed)
    title = wind.getncattr("long_name")
    units = wind.getncattr("units")
    # 为了颜色
    # alpha是透明度
    map.scatter(longs, lats, s=15, c=wind_speed, alpha=50, cmap=plt.cm.jet)
    colorbar = map.colorbar(location="bottom")
    colorbar.set_label("min:"+str(min_wind_speed)+","+"max:"+str(max_wind_speed)+" "+"units:"+units)
    # 去除刻度
    plt.xticks(())
    plt.yticks(())
    # 添加标题
    plt.title(title, fontproperties='SimHei', fontsize=15)
    plt.show()

    # 椭球型地图上画点
    # names = []
    # pops = []
    # lats = []
    # lons = []
    # countries = []
    # for line in file("major_city"):
    #     info = line.split()
    #     names.append(info[0])
    #     pops.append(float(info[1]))
    #     lat = float(info[2][:-1])
    #     if info[2][-1] == 'S': lat = -lat
    #     lats.append(lat)
    #     lon = float(info[3][:-1])
    #     if info[3][-1] == 'W': lon = -lon + 360.0
    #     lons.append(lon)
    #     country = info[4]
    #     countries.append(country)
    #     # ============================================
    #     # set up map projection with
    #     # use low resolution coastlines.
    #     map = Basemap(projection='ortho', lat_0=35, lon_0=120, resolution='l')
    #     # draw coastlines, country boundaries, fill continents.
    #     map.drawcoastlines(linewidth=0.25)
    #     map.drawcountries(linewidth=0.25)
    #     # draw the edge of the map projection region (the projection limb)
    #     map.drawmapboundary(fill_color='#689CD2')
    #     # draw lat/lon grid lines every 30 degrees.
    #     map.drawmeridians(np.arange(0, 360, 30))
    #     map.drawparallels(np.arange(-90, 90, 30))
    #     # Fill continent wit a different color
    #     map.fillcontinents(color='#BF9E30', lake_color='#689CD2', zorder=0)
    #     # compute native map projection coordinates of lat/lon grid.
    #     x, y = map(lons, lats)
    #     max_pop = max(pops)
    #     # Plot each city in a loop.
    #     # Set some parameters
    #     size_factor = 80.0
    #     y_offset = 15.0
    #     rotation = 30
    #     for i, j, k, name in zip(x, y, pops, names):
    #         size = size_factor * k / max_pop
    #         cs = map.scatter(i, j, s=size, marker='o', color='#FF5600')
    #         plt.text(i, j + y_offset, name, rotation=rotation, fontsize=10)
    # plt.title('Major Cities in Asia & Population')
    # plt.show()

    # countouf still failure becauese outof memory
    #
    # # read data from wxsj
    # nc = Dataset("F:\\wxdata\\wxdatays\\bf1ab.s20170801-000000-e20170801-235959.l2.wind.v10.nc")
    # lats = nc.variables["lat"][:][0:100].data
    # lons = nc.variables["lon"][:][0:100].data
    # wind_speed = nc.variables["wind_speed"][:][0:100]
    # print type(lats), type(lons), type(wind_speed)
    # # # print type(lats)
    # # # print type(lons)
    # # # print lats
    # # # lons = [113.14, 113.4, 121.29, 116.24]
    # # # lats = [23.08, 34.46, 31.14, 39.55]
    # Lons, Lats = np.meshgrid(lons, lats)
    #
    # # print Lons
    # # map.scatter(x,y,marker="D",color="m")
    # map = Basemap()
    # # Fill the continents with the land color
    # # map.fillcontinents(color='coral', lake_color='aqua')
    # map.drawcoastlines()
    # cs = map.contourf(Lons, Lats, wind_speed, cmap=plt.cm.jet)
    # cbar = map.colorbar()
    # map.drawcoastlines(linewidths=1.25)
    # plt.show()
