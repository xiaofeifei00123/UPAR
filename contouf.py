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
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import cmaps
from get_cmap import get_cmap_temp
import datetime
from read_data import Get_data


class Draw():
    def __init__(self, ):
        # self.da_var = da_var
        # self.name = name
        # self.station = station
        pass

    def draw_main(self, station_dic, var):

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
        # var = 'temp'
        # var = 't_td'
        # var = 'wind'
        for key in station_dic:
            pass
            station = station_dic[key]

            ### 获得数据
            gd = Get_data()
            model_dic = gd.get_data_main(var, station)
            # print("yes")

            ## 画图
            bb = self.combine_fig(var, model_dic, station['name'])

    def draw_contourf_single(self, var, ax, val, title, time_index):
        """[summary]

        Args:
            val (DataArray): 需要换图的变量
        """

        # val = self.var
        # print(val)
        x = val.time  # 各行的名称
        y = val.pressure.values

        # time_index1 = pd.date_range(start='20160702 00', end='20160731 00', freq='24H')
        # time_index2 = pd.date_range(start='20160702 00', end='20160731 00', freq='24H')
        # time_index_dic = {'00':time_index1, '12':time_index2}
        # for key in time_index_dic:
        x = time_index
        # time_index = time_index['data']
        val = val.sel(time=time_index)

        val = val.values.swapaxes(0, 1)  # 转置矩阵
        color_map = get_cmap_temp()

        if var == 'temp':
            level = np.arange(-20, 30, 2.5)
        elif var == 't_td':
            level = np.arange(0, 30, 1)
        elif var == 'wind':
            level = np.arange(0, 25, 2.5)
        elif var == 'temp_grads':
            level = np.arange(-0.3, 0.3, 0.02)
        elif var == 't_td_grads':
            # level = np.arange(0, 0.25, 0.01)
            level = np.arange(-0.25, 0.26, 0.01)
        elif var == 'wind_grads':
            level = np.arange(-0.25, 0.26, 0.05)

        # elif var in ['temp_grads', 't_td_grads', 'wind_grads']:
        #     level = np.arange(0, 0.5, 0.01)
        #     pass
        # print(len(x))
        # print(len(y))
        # print(val.shape)
        # print(level)

        CS = ax.contourf(x,
                         y,
                         val,
                         levels=level,
                         cmap=color_map,
                         extend='both')
        # cb = fig.colorbar(CS, orientation='horizontal', shrink=0.8, pad=0.14, fraction=0.14) # 这里的cs是画填色图返回的对象
        ## 设置标签大小
        # ax.set_xticks(x[::2])  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        # xlabel = x[::2].dt.strftime('%m%d').values
        ax.set_xticks(x)  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        xlabel = x.dt.strftime('%d').values
        ax.set_xticklabels(xlabel, rotation=45)

        plt.gca().invert_yaxis()
        # ax.set_ylim(0, 300)
        ax.set_title(title, fontsize=14)
        # ## 设置标签名称
        ax.set_xlabel("Time(UTC, day)", fontsize=14)
        ax.set_ylabel("Pressure (hPa)", fontsize=14)
        # fig.savefig(fig_name,bbox_inches = 'tight')
        return CS

    def combine_fig(self, var, model_dic, station_name):
        """[summary]

        Args:
            model_dic ([type]): 各试验数据的字典
        """
        # fig = plt.figure(figsize=(8, 5), dpi=400)  # 创建页面
        # ax = fig.add_axes([0.1, 0.13, 0.85, 0.8])  # 重新生成一个新的坐标图
        # print(area_dic['all'])
        fig = plt.figure(figsize=(10, 10), dpi=200)  # 创建页面
        grid = plt.GridSpec(3,
                            2,
                            figure=fig,
                            left=0.10,
                            right=0.98,
                            bottom=0.18,
                            top=0.95,
                            wspace=0.3,
                            hspace=0.35)

        axes = [None] * 6  # 设置一个维度为8的空列表
        axes[0] = fig.add_subplot(grid[0, 0:1])
        axes[1] = fig.add_subplot(grid[0, 1:2])
        axes[2] = fig.add_subplot(grid[1, 0:1])
        axes[3] = fig.add_subplot(grid[1, 1:2])
        axes[4] = fig.add_subplot(grid[2, 0:1])
        axes[5] = fig.add_subplot(grid[2, 1:2])

        model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF', 'obs']

        ax6 = fig.add_axes([0.18, 0.06, 0.7, 0.02])  # 重新生成一个新的坐标图
        keys = ['00', '12']
        # for model in model_list:
        for key in keys:
            for i in range(len(model_list)):
                time_index = model_dic[model_list[i]].time.sel(
                    time=datetime.time(int(key)))
                # time_index = time_index_dic[key]
                CS = None
                # time_index = model_dic[model_list[i]].time.sel(time=datetime.time(12))
                title = str(station_name) + "_" + model_list[i] + "_" + str(
                    key)
                CS = self.draw_contourf_single(var, axes[i],
                                               model_dic[model_list[i]], title,
                                               time_index)

            cb = fig.colorbar(CS,
                              cax=ax6,
                              orientation='horizontal',
                              shrink=0.8,
                              pad=0.14,
                              fraction=0.14)  # 这里的cs是画填色图返回的对象

            path = '/home/fengxiang/Project/Asses_PBL/Draw/UPAR/picture/'
            fig_name = os.path.join(
                path,
                str(var) + "_" + str(station_name) + "_" + str(key))
            # fig.savefig('/home/fengxiang/Project/Asses_PBL/Draw/tt.png')
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
    station = station_dic['GaiZe']
    # station = station_dic['ShiQuanhe']
    # gd = Get_data()
    # aa = gd.get_data_main('temp', station)

    # var = 'temp'
    # var = 't_td'

    ## 测试wrfout插值
    # re = Regrid()
    # aa = re.get_pressure_lev(station)
    # print(aa.values)

    # ## 测试obs插值
    # gb = Get_obs()
    # aa = gb.read_obs(station)
    # print(aa['wind_s'])
    # print(aa)

    #### 最终画图
    var_list = ['temp', 't_td', 'wind', 'temp_grads', 't_td_grads', 'wind_grads']
    var = var_list[4]

    Dr = Draw()
    Dr.draw_main(station_dic, var)
