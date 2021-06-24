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
from read_data import Get_data
# from temp import *


class Draw_skewt():
    def __init__(self, ):
        pass

    def draw_main(self, station_dic):

        time_index = ['00', '06', '12']

        for time_select in time_index:
            for key in station_dic:
                pass
                station = station_dic[key]

                ### 获得数据
                gd = Get_data()
                # model_dic = gd.get_data_main(var, station)
                # print("yes")
                model_dic_t = gd.get_data_main('temp', station)
                model_dic_td = gd.get_data_main('td', station)
                ## 画图
                title = {'time':time_select, 'station_name':station['name']}
                cc = self.draw_skewt(model_dic_t, model_dic_td, 
                                        title)


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
        time_index = time_index.values   # 转换成numpy数组
        t1 = pd.date_range('20160701 00', '20160701 1200', freq='6H')
        time_index = np.setdiff1d(time_index, t1.values)   # numpy 筛选两个array的不同部分

        for i in range(len(model_list)):
            # time_index = model_dic_t[model_list[i]].time.sel(
            #     time=datetime.time(int(key)))
            model_t[model_list[i]] = model_dic_t[model_list[i]].sel(time=time_index).mean(dim='time', skipna=True)

            model_td[model_list[i]] = model_dic_td[model_list[i]].sel(time=time_index).mean(dim='time', skipna=True)


        fig = plt.figure(figsize=(8, 4), dpi=400)

        # 画探空曲线
        skew = SkewT(fig, rotation=30, rect=[0.1, 0.1, float(8 / 17), 0.85])
        skew1 = SkewT(
            fig,
            rotation=30,
            rect=[float(8 / 17 + 0.1 + 0.1), 0.1,
                  float(5 / 17), 0.85])

        # 画T_lnP图
        num = len(model_t)
        colors = ['cyan','red', 'green', 'blue', 'orange','black', ]

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
        skew1.plot_dry_adiabats(alpha=0.3)  # 干绝热线
        skew1.plot_moist_adiabats(alpha=0.3)  # 湿绝热线
        # fig_name = fig_name+keys[0]
        path = '/home/fengxiang/Project/Asses_PBL/Draw/UPAR/picture/'
        fig_name = title_dic['station_name'] + "_" +str(title_dic['time'])
        fig_name = os.path.join(path, fig_name)
        fig.savefig(fig_name)


if __name__ == '__main__':

    pass
    station_dic = {
        'GaiZe': {
            'lat': 32.3,
            'lon': 84.0,
            'name': 'GaiZe',
            'number': '55248'
        },
        'ShenZha': {
            'lat': 30.9,
            'lon': 88.7,
            'name': 'ShenZha',
            'number': '55472'
        },
        'ShiQuanhe': {
            'lat': 32.4,
            'lon': 80.1,
            'name': 'ShiQuanhe',
            'number': '55228'
        },
    }

    #### 最终画图
    var_list = [
        'temp', 't_td', 'wind', 'temp_grads', 't_td_grads', 'wind_grads'
    ]
    # var = var_list[4]

    Dr = Draw_skewt()
    Dr.draw_main(station_dic)
