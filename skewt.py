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

# 探空图
from metpy.plots import Hodograph, SkewT
from metpy.units import units
from data_process import TransferData
# from temp import *
import metpy


class Draw_skewt():
    def __init__(self, station):
        pass
        self.station = station
        if station['name'] == 'ShenZha':
            self.model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'obs']
            self.color_list = ['orange', 'red', 'cyan', 'blue', 'black']
        else:
            self.model_list = [
                'ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF', 'obs'
            ]
            # self.color_list = [
            #     'orange', 'red', 'cyan', 'blue', 'green', 'black'
            # ]
            self.color_list = [
                'red', 'red', 'blue', 'blue', 'green', 'black'
            ]
        # self.model_list = ['obs']
        # self.model_list = ['TEMF', 'obs']
        self.line_style_list = ['solid', 'dashed','solid', 'dashed','solid', 'solid', 'solid']

    def draw_main(self, ):

        time_index = ['00', '06', '12']
        # time_index = ['12']

        # 循环出一个时次一个站点
        for time_select in time_index:
            # 获得数据
            # gd = Get_data()
            tr = TransferData(self.station)
            # print(self.station)
            # model_dic = gd.get_data_main(var, station)
            # print("yes")
            # model_dic_t = tr.transfer_data('temp')
            # model_dic_td = tr.transfer_data('td')
            model_dic_theta_v = tr.transfer_data('theta_v')
            # model_dic_t_td = tr.transfer_data('t_td')
            model_dic_q = tr.transfer_data('rh')
            model_dic_wind = tr.transfer_data('wind_s')
            # model_dic_height = tr.transfer_data('z')
            # model_dic_z = tr.transfer_data('z')
            # print(model_dic_t)
            # model_dic_td = gd.get_data_main('td', self.station)
            # model_dic_t = gd.get_data_main('temp', self.station)
            # model_dic_td = gd.get_data_main('td', self.station)
            # 画图
            title = {'time': time_select, 'station_name': self.station['name']}
            # cc = self.draw_skewt(model_dic_t, model_dic_td, title)
            self.draw_upar(
                # model_dic_t,
                model_dic_theta_v,
                model_dic_q,
                model_dic_wind,
                title,
            )
            # self.caculate_theta_v(model_dic_t, model_dic_td, model_dic_z, title)

    def draw_upar(self, model_dic_theta_v, model_dic_q,
                  model_dic_wind, title_dic):
        """画探空的廓线
        """
        pass
        # print()
        fig = plt.figure(figsize=(9, 7), dpi=400)
        ax1 = fig.add_axes([0.1, 0.15, 0.22, 0.8])  # 左下起始位置，宽度和高度
        ax2 = fig.add_axes([0.42, 0.15, 0.22, 0.8])
        ax3 = fig.add_axes([0.75, 0.15, 0.22, 0.8])

        model_list = self.model_list
        # model_t = {}  # 存储每个试验的某个时次的平均值
        # model_td = {}  # 存储每个试验的某个时次的平均值
        # model_t_td = {}
        model_q = {}
        model_wind = {}
        # model_height = {}
        model_theta_v = {}

        time_index1 = model_dic_theta_v['obs'].time.sel(
            time=datetime.time(int(title_dic['time'])))  # 00还是06,12时间筛选在这

        time_index2 = model_dic_theta_v['TEMF'].time.sel(
            time=datetime.time(int(title_dic['time'])))  # 00还是06,12时间筛选在这

                                #   time_index2.values)
        ## 求两个时间的交集
        time_index = np.intersect1d(time_index1.values,
                                    time_index2.values)
        # print(time_index)
        # 取特殊时次的平均值 ## 删除7月1日10时和7月1日12时两个时次的值的数据
        # 所有的时间序列都转换成numpy的时间序列形式，pandas和xarray 都支持
        # time_index = time_index.values  # 转换成numpy数组
        # t1 = pd.date_range('20160701 00', '20160701 1200', freq='6H')
        # # numpy筛选两个array的不同部分
        # time_index = np.setdiff1d(time_index, t1.values)
        for i in range(len(model_list)):
            # model_t[model_list[i]] = model_dic_t[model_list[i]].sel(
            #     time=time_index).mean(dim='time', skipna=True)

            model_theta_v[model_list[i]] = \
                model_dic_theta_v[model_list[i]].sel(
                    time=time_index).mean(dim='time', skipna=True)

            # model_height[model_list[i]] = \
            #     model_dic_height[ model_list[i]].sel(
            #     time=time_index).mean(dim='time', skipna=True)

            # model_t_td[model_list[i]] = \
            #     model_dic_t_td[model_list[i]].sel(
            #         time=time_index).mean(dim='time', skipna=True)
            model_q[model_list[i]] = \
                model_dic_q[model_list[i]].sel(
                    time=time_index).mean(dim='time', skipna=True)

            model_q['fnl'] = \
                model_dic_q['fnl'].sel(
                    time=time_index).mean(dim='time', skipna=True)
            model_wind[model_list[i]] = \
                model_dic_wind[model_list[i]].sel(
                    time=time_index).mean(dim='time', skipna=True)

        # color_list = ['orange', 'red', 'cyan', 'blue', 'green', 'black']
        color_list = self.color_list
        # color_list = ['orange', 'red', 'cyan', 'blue', 'black']

        # theta_v = self.caculate_theta_v_height(model_t, model_td, model_height)

        # model_t = theta_v
        model_t = model_theta_v
        model_t_td = model_q

        # aa = model_t['YSU'].to_series()
        # aa1 = aa.dropna()
        # y_bottom = aa1.index[0]

        bb = model_t_td['YSU'].to_series()
        bb1 = bb.dropna()
        y2_bottom = bb1.index[0]

        cc = model_wind['YSU'].to_series()
        cc1 = cc.dropna()
        y3_bottom = cc1.index[0]

        for [model, i, j] in zip(model_list, color_list, self.line_style_list):

            # y = model_t[model].height_coord
            y = model_t[model].pressure
            x = model_t[model].values
            ax1.plot(x, y, color=i, label=model, linestyle=j)

            y2 = model_t_td[model].pressure
            # y2 = model_t_td[model].height_coord
            x2 = model_t_td[model].values
            # y2 = y2_bottom - y2
            ax2.plot(x2, y2, color=i, label=model, linestyle=j)

            # ax2.plot(x2, y2, color=i, label=model)

            y3 = model_wind[model].pressure
            # y3 = model_wind[model].height_coord
            x3 = model_wind[model].values
            # y3 = y3_bottom - y3
            ax3.plot(x3, y3, color=i, label=model, linestyle=j)
            ax3.set_xlim(0, 25)

        fts = 14
        # ax1.set_ylabel('Height (m)', fontsize=fts)
        ax1.set_ylabel('Pressure/(hPa)', fontsize=fts)
        ax1.set_xlabel('Theta_v/(K)', fontsize=fts)
        # ax1.set_ylim(0, 5000)
        # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax1.set_xticks(np.arange(320, 355, 10))
        ax1.set_xlim(320, 360)

        # ax2.set_xlabel('q/(g/kg)', fontsize=fts)
        ax2.set_xlabel('rh', fontsize=fts)
        # ax2.set_xticks(np.arange(0, 26, 5))
        # ax2.set_ylim(0, 5000)

        x2_fnl = model_t_td['fnl'].values*0.01
        y2_fnl = model_t_td['fnl'].pressure
        ax2.plot(x2_fnl, y2_fnl, color='cyan', label='fnl', linestyle='solid')
        ax2.legend(loc='lower center',
                   bbox_to_anchor=(0.25, -0.18, 0.2, 0.2),
                   ncol=4)

        ax3.set_xlabel('Wind_speed (m/s)', fontsize=fts)
        ax3.set_xticks(np.arange(0, 26, 5))
        # ax3.set_ylim(0, 5000)

        ax1.invert_yaxis()
        ax2.invert_yaxis()
        ax3.invert_yaxis()
        ax1.set_ylim(570, 300)
        ax2.set_ylim(570, 300)
        ax3.set_ylim(570, 300)

        fig.suptitle(
            title_dic['station_name'] + "_" + title_dic['time'],
            fontsize=fts * 1.3
        )

        fig_name = title_dic['station_name'] + "_" + title_dic['time']
        path = '/home/fengxiang/Project/Asses_PBL/Draw/UPAR/'
        fgnm = os.path.join(path, fig_name)
        fig.savefig(fgnm)


if __name__ == '__main__':

    pass
    station_dic = {
        # 'GaiZe': {
        #     'lat': 32.3,
        #     'lon': 84.0,
        #     'name': 'GaiZe',
        #     'number': '55248',
        #     'height': 4400,
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
    # 最终画图
    var_list = [
        'temp', 't_td', 'wind', 'temp_grads', 't_td_grads', 'wind_grads'
    ]
    # var = var_list[4]
    # station = station_dic['ShiQuanhe']
    #    station = station_dic['GaiZe']

    for key in station_dic:
        station = station_dic[key]
        Dr = Draw_skewt(station)
        Dr.draw_main()
