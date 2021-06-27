#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
单变量的插值和绘图
读temp, td, u, v
读pressure
并插值到站点上
垂直插值+水平插值
画一个站点的图
每天所有时次的
原始气压坐标
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
    def __init__(self, ):
        pass

    def draw_main(self, station_dic):

        time_index = ['00', '06', '12']
        # time_index = ['00']

        ## 一个时次一个站点
        for time_select in time_index:
            for key in station_dic:
                pass
                station = station_dic[key]

                ### 获得数据
                # gd = Get_data()
                tr = TransferData(station)
                # model_dic = gd.get_data_main(var, station)
                # print("yes")
                model_dic_t = tr.transfer_data('temp')
                model_dic_t_td = tr.transfer_data('t_td')
                model_dic_wind = tr.transfer_data('wind')
                # print(model_dic_t)
                # model_dic_td = gd.get_data_main('td', station)
                # model_dic_t = gd.get_data_main('temp', station)
                # model_dic_td = gd.get_data_main('td', station)
                ## 画图
                title = {'time': time_select, 'station_name': station['name']}
                # cc = self.draw_skewt(model_dic_t, model_dic_td, title)
                self.draw_upar(model_dic_t, model_dic_t_td, model_dic_wind,
                               title)

    def draw_upar(self, model_dic_t, model_dic_t_td, model_dic_wind,
                  title_dic):
        """画探空的廓线
        """
        pass
        # print()
        fig = plt.figure(figsize=(8, 7), dpi=400)
        # ax1 = fig.add_axes([0.76, 0.1, 0.23, 0.85])  # 左右起始位置，宽度和高度
        ax1 = fig.add_axes([0.1, 0.1, 0.22, 0.8])  # 左右起始位置，宽度和高度
        ax2 = fig.add_axes([0.42, 0.1, 0.22, 0.8])  # 左右起始位置，宽度和高度
        ax3 = fig.add_axes([0.75, 0.1, 0.22, 0.8])  # 左右起始位置，宽度和高度

        model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF', 'obs']
        # model_list = ['TEMF', 'obs']
        model_t = {}  # 存储每个试验的某个时次的平均值
        model_t_td = {}
        model_wind = {}

        time_index = model_dic_t['obs'].time.sel(
            time=datetime.time(int(title_dic['time'])))

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

            model_t_td[model_list[i]] = model_dic_t_td[model_list[i]].sel(
                time=time_index).mean(dim='time', skipna=True)

            model_wind[model_list[i]] = model_dic_wind[model_list[i]].sel(
                time=time_index).mean(dim='time', skipna=True)

        color_list = ['orange', 'red', 'cyan', 'blue', 'green', 'black']

        
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
        ax1.legend(loc='lower left')
        ax1.set_ylabel('Pressure off the ground (hPa)', fontsize=fts)
        ax1.set_xlabel('T', fontsize=fts)
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

    def caculate_theta_v(self, model_dic_t, model_dic_td, title_dic):
        # -------------获取数据-------------------
        model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF', 'obs']
        model_t = {}  # 存储每个试验的某个时次的平均值
        model_td = {}

        time_index = model_dic_t['obs'].time.sel(
            time=datetime.time(int(title_dic['time'])))

        print(time_index)
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

        metpy.calc.virtual_potential_temperature()

    def draw_skewt(slef, model_dic_t, model_dic_td, title_dic):
        """画skewt图

        Args:
            slef ([type]): [description]
            model_dic_t ([type]): 温度文
            model_dic_td ([type]): [description]
            title_dic ([type]): [description]
        """

        # -------------获取数据-------------------
        model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF', 'obs']
        model_t = {}
        model_td = {}

        time_index = model_dic_t['obs'].time.sel(
            time=datetime.time(int(title_dic['time'])))

        print(time_index)
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

        fig = plt.figure(figsize=(8, 4), dpi=400)

        ###  T_lnP & Td_lnP
        skew = SkewT(fig, rotation=30, rect=[0.1, 0.1, float(8 / 17), 0.85])
        skew1 = SkewT(
            fig,
            rotation=30,
            rect=[float(8 / 17 + 0.1 + 0.1), 0.1,
                  float(5 / 17), 0.85])

        # 画T_lnP图
        num = len(model_t)
        colors = [
            'cyan',
            'red',
            'green',
            'blue',
            'orange',
            'black',
        ]

        print("画T_lnP")
        for i in range(num):  # 模式个数+观测
            # print(i)
            # line = skew.plot(module.keys[i])
            key = list(model_t.keys())
            val = list(model_t.values())
            module_value = val[i]
            module_name = key[i]
            # print(module_name)
            # val = modlue.values()[i]
            # print(key[i])
            # print(val)
            # line = skew.plot_barbs(val[''])
            line = skew.plot(module_value['pressure'],
                             module_value.values,
                             linewidth=2.0,
                             label=module_name,
                             color=colors[i])

        print("画Td_lnP")
        for i in range(num):  # 模式个数+观测
            # print(i)
            # line = skew.plot(module.keys[i])
            key = list(model_td.keys())
            val = list(model_td.values())
            module_value = val[i]
            module_name = key[i]
            # print(module_name)
            line = skew.plot(module_value['pressure'],
                             module_value.values,
                             linestyle='--',
                             linewidth=2.0,
                             color=colors[i])

        for i in range(num):  # 模式个数+观测
            key = list(model_t.keys())
            val_t = list(model_t.values())
            val_td = list(model_td.values())
            module_value_t = val_t[i]
            module_value_td = val_td[i]
            module_name = key[i]
            # print(module_name)
            line = skew1.plot(
                module_value_t['pressure'],
                module_value_t.values - module_value_td.values,
                linestyle='-.',
                linewidth=2.,
                alpha=1.,
                # marker='*',
                color=colors[i])
            # line.set_linewidht(2.0)
        skew.ax.set_xlabel('T,Td (℃)', fontsize=14)
        skew.ax.set_ylabel('Pressure (hPa)', fontsize=14)
        skew.ax.set_ylim(700, 100)
        skew.ax.set_xlim(-50, 30)
        skew.plot_dry_adiabats(alpha=0.3)  # 干绝热线
        skew.plot_moist_adiabats(alpha=0.3)  # 湿绝热线
        skew.ax.legend(edgecolor='white')

        skew.ax.set_title(title_dic['station_name'], loc='right', fontsize=14)
        skew.ax.set_title(title_dic['time'], loc='left', fontsize=14)
        skew.ax.xaxis.set_tick_params(labelsize=12)
        skew.ax.yaxis.set_tick_params(labelsize=12)

        skew1.ax.set_xlabel('T-Td (℃)', fontsize=14)
        skew1.ax.set_ylim(700, 100)
        skew1.ax.set_xlim(0, 50)
        # skew1.plot_dry_adiabats(alpha=0.3)  # 干绝热线
        # skew1.plot_moist_adiabats(alpha=0.3)  # 湿绝热线
        # fig_name = fig_name+keys[0]
        path = '/home/fengxiang/Project/Asses_PBL/Draw/UPAR/picture/'
        fig_name = title_dic['station_name'] + "_" + str(title_dic['time'])
        fig_name = os.path.join(path, fig_name)
        fig.savefig(fig_name)


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
        'ShenZha': {
            'lat': 30.9,
            'lon': 88.7,
            'name': 'ShenZha',
            'number': '55472',
            'height': 4672
        },
        # 'ShiQuanhe': {
        #     'lat': 32.4,
        #     'lon': 80.1,
        #     'name': 'ShiQuanhe',
        #     'number': '55228',
        #     'height': 4280
        # },
    }
    #### 最终画图
    var_list = [
        'temp', 't_td', 'wind', 'temp_grads', 't_td_grads', 'wind_grads'
    ]
    # var = var_list[4]

    Dr = Draw_skewt()
    Dr.draw_main(station_dic)
