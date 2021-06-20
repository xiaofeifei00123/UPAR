#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读temp
读pressure
并插值到站点上
垂直插值+水平插值

600, 590, 580, 570, 560, 550, 540, 530, 520, 510, 500,
475, 450, 425, 400, 350, 300, 200, 100
-----------------------------------------
Time             :2021/06/17 18:59:18
Author           :Forxd
Version          :2.0
'''

import xarray as xr
import os
import numpy as np
import pandas as pd
import wrf
import netCDF4 as nc
from netCDF4 import Dataset
from wrf import getvar, vinterp, interplevel
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,FormatStrFormatter
import cmaps 
class Regrid():

    def get_pressure_lev(self, ):
        #### 得到各层的pressure值
        ## 不同时刻各层气压值，差别可以忽略不计, 
        # 后面还要对气压层进行插值, 这里不对它做过高精度要求
        path = '/mnt/zfm_18T/Asses_PBL/wrfout_data/'
        flnm_pressure = os.path.join(path, 'pressure_Jul_YSU_latlon')
        ds_pressure = xr.open_dataset(flnm_pressure)
        pr = ds_pressure.pressure
        prb = pr.sel(time='2016-07-01 13:00')
        prc = prb.sel(lat=32.13, lon=92.5, method='nearest')
        return prc

    def regrid(self, da_temp, station, pressure_lev):
        """[summary]

        Args:
            da_temp ([DataArray]): [需要插值的变量]

        Returns:
            [type]: [description]
        """

        #### 将bottom_top坐标换成气压坐标 
        time_coord = da_temp.time.values
        lat_coord = da_temp.lat.values
        lon_coord = da_temp.lon.values

        prc = self.get_pressure_lev()
        pressure_coord = prc.values
        # print(pressure_coord)

        da_temp_reset = da_temp.values
        # print(da_temp_reset)
        da_temp = xr.DataArray(da_temp_reset, 
                coords=[time_coord, pressure_coord, lat_coord, lon_coord],
                dims=['time', 'pressure', 'lat', 'lon']
                )

        #### 对换成气压坐标的temp进行插值
        ## 水平插值
        # da_temp = da_temp.sel(lat=33.5, lon=97.5, method='nearest')
        da_temp = da_temp.sel(lat=station['lat'], lon=station['lon'], method='nearest')
        ## 垂直插值
        da_temp_return = da_temp.interp(pressure=pressure_lev)
        # print(da_temp_return)
        return da_temp_return


class Draw():

    def __init__(self, da_var, name):
        self.var = da_var
        self.name = name
        self.station = station
        pass

    def draw_main(self,):
        fig = plt.figure(figsize=(7,4.5), dpi=400)
        ax = fig.add_axes([0.12,0.01,0.85,0.9])
        self.draw_contourf(ax, fig)
        # self.draw_barb(ax,fig)
        fig.savefig(self.name['fig_name'])
        pass

    def draw_contourf(self, ax, fig):

        val = self.var
        # print(val)
        x = val.time  # 各行的名称
        y = val.pressure.values


        aa = x.values
        bb = pd.Series(aa)
        x = bb

        val = val.values.swapaxes(0,1)  # 转置矩阵
        color_map = cmaps.BlueDarkRed18
        # color_map = cmaps.GMT_panoply

        level = np.arange(-20, 21, 1)
        CS = ax.contourf(x,y,val,levels=level, cmap=color_map, extend='both')
        cb = fig.colorbar(CS, orientation='horizontal', shrink=0.8, pad=0.14, fraction=0.14) # 这里的cs是画填色图返回的对象
    #    ## 设置标签大小
        ax.set_xticks(x[::24])  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax.set_xticklabels(x[::24].dt.strftime('%d'))

        # ax.gca().invert_yaxix()
        plt.gca().invert_yaxis()

        ax.set_yticks(y[::5])
        ax.set_yticklabels(y[::5])
        ax.set_title(self.name['title'])

        # ## 设置标签名称
        ax.set_xlabel("Time(UTC, day)", fontsize=14)
        ax.set_ylabel("Pressure (hPa)", fontsize=14)

        fig.savefig(fig_name,bbox_inches = 'tight')

class Main():
    pass

if __name__ == '__main__':


    station={'lat':32.3, 'lon':84.0, 'name':'GaiZe'}
    # station={'lat':30.56, 'lon':88.42, 'name':'ShenZha'}

    pressure_level = np.arange(300, 550+1, 5)
    var = 'temp'
    month = 'Jul'
    model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
    model = model_list[0]

    path = '/mnt/zfm_18T/Asses_PBL/wrfout_data/'
    # flnm_temp = os.path.join(path, 'temp_Jul_YSU_latlon')
    file_name = str(var)+"_"+str(month)+"_"+str(model)+"_latlon"
    flnm_var = os.path.join(path, file_name)

    ds_var = xr.open_dataset(flnm_var)
    da_var = ds_var[var]


    Re = Regrid()
    da_var = Re.regrid(da_var, station, pressure_level)


    # print(da_var)
    fig_path = '/home/fengxiang/Project/Asses_PBL/Draw/UPAR/'
    fig_name = os.path.join(fig_path, str(model)+"_"+station['name'])
    title = str(model)+'_'+str(station['name'])
    name = {'fig_name':fig_name, 'title':title}

    Dr = Draw(da_var, name)
    Dr.draw_main()