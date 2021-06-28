#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
画探空的廓线图
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
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import cmaps
from get_cmap import get_cmap_temp
import datetime

## 探空图
from metpy.plots import Hodograph, SkewT
from metpy.units import units
from data_process import TransferData
# from temp import *
import metpy


class Draw_skewt():
    def __init__(self, station):
        pass
        self.station = station

    def draw_main(self, ):

        time_index = ['00', '06', '12']
        # time_index = ['00']

        ## 循环出一个时次一个站点
        for time_select in time_index:
            ### 获得数据
            # gd = Get_data()
            tr = TransferData(self.station)
            # model_dic = gd.get_data_main(var, station)
            # print("yes")
            model_dic_t = tr.transfer_data('temp')
            model_dic_td = tr.transfer_data('td')
            # model_dic_t_td = tr.transfer_data('t_td')
            # model_dic_wind = tr.transfer_data('wind')
            model_dic_z = tr.transfer_data('z')
            # model_dic_z = tr.transfer_data('z')
            # print(model_dic_t)
            # model_dic_td = gd.get_data_main('td', self.station)
            # model_dic_t = gd.get_data_main('temp', self.station)
            # model_dic_td = gd.get_data_main('td', self.station)
            ## 画图
            title = {'time': time_select, 'station_name': self.station['name']}
            # cc = self.draw_skewt(model_dic_t, model_dic_td, title)
            # self.draw_upar(
            #     model_dic_t, 
            #     model_dic_td,
            #     model_dic_t_td,
            #     model_dic_wind,
            #     title,
            #     )
            self.caculate_theta_v(model_dic_t, model_dic_td, model_dic_z, title)

    def draw_upar(self, model_dic_t, model_dic_td, model_dic_t_td,
                  model_dic_wind, title_dic):
        """画探空的廓线
        """
        pass
        # print()
        fig = plt.figure(figsize=(8, 7), dpi=400)
        # ax1 = fig.add_axes([0.76, 0.1, 0.23, 0.85])  # 左右起始位置，宽度和高度
        ax1 = fig.add_axes([0.1, 0.1, 0.22, 0.8])  # 左右起始位置，宽度和高度
        ax2 = fig.add_axes([0.42, 0.1, 0.22, 0.8])  # 左右起始位置，宽度和高度
        ax3 = fig.add_axes([0.75, 0.1, 0.22, 0.8])  # 左右起始位置，宽度和高度

        # model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF', 'obs']
        model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'obs']
        # model_list = ['TEMF', 'obs']
        model_t = {}  # 存储每个试验的某个时次的平均值
        model_td = {}  # 存储每个试验的某个时次的平均值
        model_t_td = {}
        model_wind = {}

        time_index = model_dic_t['obs'].time.sel(time=datetime.time(
            int(title_dic['time'])))  # 00还是06,12时间筛选在这

        #### 取特殊时次的平均值
        ## 删除7月1日10时和7月1日12时两个时次的值的数据
        ## 所有的时间序列都转换成numpy的时间序列形式，pandas和xarray 都支持
        time_index = time_index.values  # 转换成numpy数组
        t1 = pd.date_range('20160701 00', '20160701 1200', freq='6H')
        ## numpy筛选两个array的不同部分
        time_index = np.setdiff1d(time_index, t1.values)
        for i in range(len(model_list)):
            model_t[model_list[i]] = model_dic_t[model_list[i]].sel(
                time=time_index).mean(dim='time', skipna=True)

            model_td[model_list[i]] = model_dic_td[model_list[i]].sel(
                time=time_index).mean(dim='time', skipna=True)

            model_t_td[model_list[i]] = model_dic_t_td[model_list[i]].sel(
                time=time_index).mean(dim='time', skipna=True)

            model_wind[model_list[i]] = model_dic_wind[model_list[i]].sel(
                time=time_index).mean(dim='time', skipna=True)

        # color_list = ['orange', 'red', 'cyan', 'blue', 'green', 'black']
        color_list = ['orange', 'red', 'cyan', 'blue', 'black']

        theta_v = self.caculate_theta_v(model_t, model_td)

        model_t = theta_v

        aa = model_t['YSU'].to_series()
        aa1 = aa.dropna()
        y_bottom = aa1.index[0]

        bb = model_t_td['YSU'].to_series()
        bb1 = bb.dropna()
        y2_bottom = bb1.index[0]

        cc = model_wind['YSU'].to_series()
        cc1 = cc.dropna()
        y3_bottom = cc1.index[0]

        for [model, i] in zip(model_list, color_list):

            y = model_t[model].pressure
            x = model_t[model].values
            # y_bottom = model_t['YSU'].pressure[0]
            # y_bottom = bottom1
            y = y_bottom - y
            # print(model_t['TEMF'])

            # ax1.invert_yaxis()
            ax1.plot(x, y, color=i, label=model)
            ax1.set_ylim(0, 500)
            # plt.savefig('testttt')

            y2 = model_t_td[model].pressure
            x2 = model_t_td[model].values
            # y2_bottom = model_t_td['YSU'].pressure[0]
            y2 = y2_bottom - y2

            # ax1.invert_yaxis()
            ax2.plot(x2, y2, color=i, label=model)
            ax2.set_ylim(0, 500)

            y3 = model_wind[model].pressure
            x3 = model_wind[model].values
            # y3_bottom = model_wind['YSU'].pressure[0]
            y3 = y3_bottom - y3

            # ax1.invert_yaxis()
            ax3.plot(x3, y3, color=i, label=model)
            ax3.set_ylim(0, 500)
            ax3.set_xlim(0, 25)
            # ax3.set_xlim(0,25)

        fts = 14
        ax1.legend(loc='lower right')
        ax1.set_ylabel('Pressure off the ground (hPa)', fontsize=fts)
        ax1.set_xlabel('theta_v', fontsize=fts)
        # ax1.invert_yaxis()
        # ax1.set_xticks(np.arange(-75, 26, 25))  # 这个是选择哪几个坐标画上来的了,都有只是显不显示

        ax2.set_xlabel('T-Td', fontsize=fts)
        # ax2.set_xticks(x2[::5])  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax2.set_xticks(np.arange(0, 26, 5))  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax3.set_xlabel('Wind_speed', fontsize=fts)
        ax3.set_xticks(np.arange(0, 26, 5))  # 这个是选择哪几个坐标画上来的了,都有只是显不显示

        # ax.set_xticks(x)  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        # xlabel = x.dt.strftime('%d').values
        # ax.set_xticklabels(xlabel, rotation=45)

        # ax2.legend()
        # ax3.legend()
        fig.suptitle(title_dic['station_name'] + "_" + title_dic['time'],
                     fontsize=fts * 1.3)

        fig_name = title_dic['station_name'] + "_" + title_dic['time']
        path = '/home/fengxiang/Project/Asses_PBL/Draw/UPAR/'
        fgnm = os.path.join(path, fig_name)

        fig.savefig(fgnm)

    def caculate_theta_v(self, model_dic_t, model_dic_td, model_dic_z, title_dic):
        # -------------获取数据-------------------

        ### TEMF方案温度是缺省值是怎么回事呢????

        # model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'obs']
        # model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF', 'obs']
        model_list = ['TEMF']
        model_t = {}  # 存储每个试验的某个时次的平均值
        model_td = {}
        model_z = {}
        q = {}  # specific humidity
        w = {}  # mixing ratio
        theta_v = {}  # virtual potential temperature

        time_index = model_dic_t['obs'].time.sel(
            time=datetime.time(int(title_dic['time'])))   # 00还是06,12时间筛选在这

        time_index = time_index.values  # 转换成numpy数组
        t1 = pd.date_range('20160701 00', '20160701 1200', freq='6H')
        ## numpy筛选两个array的不同部分
        time_index = np.setdiff1d(time_index, t1.values)
        # print(time_index)

        for i in range(len(model_list)):
            model_t[model_list[i]] = model_dic_t[model_list[i]].sel(
                time=time_index).mean(dim='time', skipna=True)
            model_td[model_list[i]] = model_dic_td[model_list[i]].sel(
                time=time_index).mean(dim='time', skipna=True)
            model_z[model_list[i]] = model_dic_z[model_list[i]].sel(
                time=time_index).mean(dim='time', skipna=True)

        dic_return = {}

        for model in model_list:
            pass
            td = model_td[model].dropna(dim='pressure')
            prc = td.pressure.values
            height = model_z[model].sel(pressure=prc).values - self.station['height']
            t = model_t[model].sel(pressure=prc).dropna(dim='pressure')
            if len(t)==0 or len(td)==0:
                pass
                da = model_t[model]
                dic_return[model] = da
            else: 
                pressure = units.Quantity(td.pressure.values,"hPa")

                temperature = units.Quantity(t.values,"degC")
                dew_point = units.Quantity(td.values,"degC")

                # if temperature[0]

                # q[model] = metpy.calc.specific_humidity_from_dewpoint(model_td[model].pressure.values, model_td[model].values)
                q[model] = metpy.calc.specific_humidity_from_dewpoint(pressure, dew_point)
                # print(q[model])
                w[model] = metpy.calc.mixing_ratio_from_specific_humidity(q[model])
                theta_v[model] = metpy.calc.virtual_potential_temperature(pressure, temperature, w[model])
                # da = xr.DataArray(theta_v[model], coords=prc, dim=pressure)
                # da = xr.DataArray(theta_v[model], coords=[prc], dims=['pressure'])
                da = xr.DataArray(theta_v[model], coords=[height], dims=['height'])
                # print(da)
                dic_return[model] = da
                # print(theta_v[model])

        return dic_return

    # def caculate_theta_v(self, model_t, model_td):
    #     # -------------获取数据-------------------

    #     q = {}
    #     w = {}
    #     theta_v = {}

    #     dic_return = {}

    #     model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'obs']
    #     for model in model_list:
    #         pass
    #         td = model_td[model].dropna(dim='pressure')
    #         prc = td.pressure.values
    #         t = model_t[model].sel(pressure=prc).dropna(dim='pressure')

    #         pressure = units.Quantity(td.pressure.values, "hPa")
    #         temperature = units.Quantity(t.values, "degC")
    #         dew_point = units.Quantity(td.values, "degC")

    #         # q[model] = metpy.calc.specific_humidity_from_dewpoint(model_td[model].pressure.values, model_td[model].values)
    #         q[model] = metpy.calc.specific_humidity_from_dewpoint(
    #             pressure, dew_point)
    #         # print(q[model])
    #         w[model] = metpy.calc.mixing_ratio_from_specific_humidity(q[model])
    #         theta_v[model] = metpy.calc.virtual_potential_temperature(
    #             pressure, temperature, w[model])
    #         # da = xr.DataArray(theta_v[model], coords=prc, dim=pressure)
    #         da = xr.DataArray(theta_v[model], coords=[prc], dims=['pressure'])
    #         print(da)
    #         dic_return[model] = da
    #         # print(theta_v[model])

    #     return dic_return


if __name__ == '__main__':

    pass
    station_dic = {
        # 'GaiZe': {
            # 'lat': 32.3,
            # 'lon': 84.0,
            # 'name': 'GaiZe',
            # 'number': '55248',
            # 'height': 4400,
        # },
        # 'ShenZha': {
        #     'lat': 30.9,
        #     'lon': 88.7,
        #     'name': 'ShenZha',
        #     'number': '55472',
        #     'height': 4672
        # },
        'ShiQuanhe': {
            'lat': 32.4,
            'lon': 80.1,
            'name': 'ShiQuanhe',
            'number': '55228',
            'height': 4280
        },
    }
    #### 最终画图
    var_list = [
        'temp', 't_td', 'wind', 'temp_grads', 't_td_grads', 'wind_grads'
    ]
    # var = var_list[4]
    station = station_dic['GaiZe']
    Dr = Draw_skewt(station)
    Dr.draw_main(station_dic)
