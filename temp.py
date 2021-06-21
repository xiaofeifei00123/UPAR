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


class Upar():
    pass

    def __init__(var):
        pass

class Get_obs():

    def read_single(self, flnm):
        pass
        """这里出现的这些问题，都是由于自己对pandas库不熟悉所导致的，
        我希望自己的工作，不仅仅是打工者的角色, 要做到更多更好一点才好
        """
        col_names = ['pressure', 'height', 'temp', 'td', 'wind_d', 'wind_s']
        ## 按列数据, 要哪几列的数据, 再命名
        df = pd.read_table(flnm,
                            sep=' ',
                        #  skiprows=0,
                        usecols=[26, 27, 29, 30, 32, 33],
                        names=col_names)

        df1 = df.where(df<9999,np.nan)  # 将缺省值赋值为NaN
        df2 = df1.dropna(axis=0, subset=['pressure'])  # 将含有缺省值的行删掉
        df2 = df2.dropna(axis=0, subset=['temp'])
        ##　这里是导致自己花了大量时间的原因
        df3 = df2.drop_duplicates('pressure', 'first', inplace=True)  # 将preuusre这一列中，含有相同值的取第一个，其他删掉
        df3 = df2.set_index(['pressure'], inplace=True)  # 将pressure这一列设为index
        df3 = df2.sort_values('pressure')   # 按照某一列排序
        ds = xr.Dataset.from_dataframe(df3)
        pressure_level = [
            600, 575, 550, 525, 500, 450, 400, 350, 300, 250, 200, 150, 100
        ]
        cc = ds.interp(pressure=pressure_level)
        return cc

    def read_obs(self, station):
        number = station['number']
        # path = "/mnt/zfm_18T/Asses_PBL/GPS_Upar_2016/SCEX_TIPEX3_UPAR_GPS_MUL_55228-201607/"
        path1 = "/mnt/zfm_18T/Asses_PBL/GPS_Upar_2016/SCEX_TIPEX3_UPAR_GPS_MUL_"
        path = path1+str(number)+"-201607"

        aa = os.listdir(path)  # 文件名列表
        # print(type(aa))
        aa.sort()  # 排一下顺序，这个是对列表本身进行操作

        ds_time = []  # 每个时次变量的列表
        ttt = []   # 时间序列列表
        for flnm in aa:

            fl_time = flnm[-14:-4]
            tt = pd.to_datetime(fl_time, format='%Y%m%d%H')
            ttt.append(tt)
            ## 这时间是不规则的
            flnm = os.path.join(path, flnm)

            aa = self.read_single(flnm)
            ds_time.append(aa)

        ds = xr.concat(ds_time,dim='time')
        ds.coords['time'] = ttt
        # print(ds['temp'])
        return ds

class Get_data():
    """
    读取，处理数据
    """
    def __init__(self, ):
        pass

    def get_data_single(self,
                        station={
                            'lat': 32.4,
                            'lon': 80.1,
                            'name': 'ShiQuanhe',
                            'number': '55228'
                        }):
        """读取单个变量,5次wrf试验的数据
        """
        # station={'lat':32.3, 'lon':84.0, 'name':'GaiZe'}
        #### 循环出需要进行处理的单个试验
        pressure_level = [
            600, 575, 550, 525, 500, 450, 400, 350, 300, 250, 200, 150, 100
        ]
        ## 站点
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
        var = 'temp'  # 先就处理温度, 后面的数据是需要处理得到的
        month = 'Jul'  # 后面进行循环即可
        model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        # model = model_list[0]

        path = '/mnt/zfm_18T/Asses_PBL/wrfout_data/'
        # station = station_dic['GaiZe']

        model_dic = {}
        for model in model_list:
            # flnm_temp = os.path.join(path, 'temp_Jul_YSU_latlon')
            file_name = str(var) + "_" + str(month) + "_" + str(
                model) + "_latlon"
            flnm_var = os.path.join(path, file_name)

            ds_var = xr.open_dataset(flnm_var)
            da_var = ds_var[var]

            Re = Regrid()
            da_var = Re.regrid(da_var, station, pressure_level)
            model_dic[model] = da_var
        return model_dic

    def get_data_single_once(self, var, flnm_var, station):
        '''单个变量，单个试验的数据读取'''
        pass
        pressure_level = [
            600, 575, 550, 525, 500, 450, 400, 350, 300, 250, 200, 150, 100
        ]

        ds_var = xr.open_dataset(flnm_var)
        da_var = ds_var[var]

        Re = Regrid()
        da_var = Re.regrid(da_var, station, pressure_level)
        return da_var
        # model_dic[model] = da_var
        # return model_dic

    def get_data_t_td(self, station):
        pass
        model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        var1 = 'temp'  # 先就处理温度, 后面的数据是需要处理得到的
        var2 = 'td'  # 先就处理温度, 后面的数据是需要处理得到的
        month = 'Jul'
        model_dic = {}
        for model in model_list:
            pass
            path = '/mnt/zfm_18T/Asses_PBL/wrfout_data/'
            file_name1 = str(var1) + "_" + str(month) + "_" + str(
                model) + "_latlon"
            flnm_var1 = os.path.join(path, file_name1)
            file_name2 = str(var2) + "_" + str(month) + "_" + str(
                model) + "_latlon"
            flnm_var2 = os.path.join(path, file_name2)
            temp = self.get_data_single_once(var1, flnm_var1, station)
            td = self.get_data_single_once(var2, flnm_var2, station)
            model_dic[model] = temp - td
        return model_dic


    def get_data_wind(self, station):
        pass
        model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        var1 = 'U'  # 先就处理温度, 后面的数据是需要处理得到的
        var2 = 'V'  # 先就处理温度, 后面的数据是需要处理得到的
        month = 'Jul'
        model_dic = {}
        for model in model_list:
            pass
            path = '/mnt/zfm_18T/Asses_PBL/wrfout_data/'
            file_name1 = str(var1) + "_" + str(month) + "_" + str(
                model) + "_latlon"
            flnm_var1 = os.path.join(path, file_name1)
            file_name2 = str(var2) + "_" + str(month) + "_" + str(
                model) + "_latlon"
            flnm_var2 = os.path.join(path, file_name2)
            U = self.get_data_single_once(var1, flnm_var1, station)
            V = self.get_data_single_once(var2, flnm_var2, station)
            model_dic[model] = xr.ufuncs.sqrt(U**2 + V**2)
        return model_dic

    def get_data_temp(self, station):
        pass
        model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        var = 'temp'  # 先就处理温度, 后面的数据是需要处理得到的
        # var2 = 'V'  # 先就处理温度, 后面的数据是需要处理得到的
        month = 'Jul'
        model_dic = {}
        for model in model_list:
            pass
            path = '/mnt/zfm_18T/Asses_PBL/wrfout_data/'
            file_name = str(var) + "_" + str(month) + "_" + str(
                model) + "_latlon"
            flnm_var = os.path.join(path, file_name)
            temp = self.get_data_single_once(var, flnm_var, station)
            model_dic[model] = temp
        return model_dic

    def grads_data(self, model_dic):
        """根据输入的垂直文件，求它做梯度后的文件

        Args:
            model_dic ([type]): [description]
        """
        pass
        # print("yes")
        model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF', 'obs']
        dic_return = {}
        for model in model_list:
            ser = []
            da = model_dic[model]
            pressure = da.pressure.values
            ## 第0层是600hPa
            for i in range(len(pressure)):
                if i == 0:
                    sr = (da[:,i+1]-da[:,i])/(pressure[i+1]-pressure[i])  # 边界上使用前差或后差
                    ser.append(sr)
                elif i>0 and i<len(pressure)-1:
                    sr = (da[:,i+1]-da[:,i-1])/(pressure[i]-pressure[i-1])/2   # 中间的使用中央差分
                    ser.append(sr)
                elif i==len(pressure)-1:
                    pass
                    sr = (da[:,i]-da[:,i-1])/(pressure[i]-pressure[i-1])
                    ser.append(sr)
            da = xr.concat(ser, 'pressure')
            da.coords['pressure'] = pressure
            dic_return[model] = da
        return dic_return


    def get_data_main(self, var, station):
        pass
        gb = Get_obs()
        ds_obs = gb.read_obs(station)  # 观测数据
        if var == 'temp':
            model_dic = self.get_data_temp(station)
            model_dic['obs'] = ds_obs['temp']
            print(model_dic)
        elif var == 't_td':
            pass
            model_dic = self.get_data_t_td(station)
            model_dic['obs'] = ds_obs['temp'] - ds_obs['td']
            print(model_dic)
        elif var == 'wind':
            model_dic = self.get_data_wind(station)
            model_dic['obs'] = ds_obs['wind_s']
            print(model_dic)
            pass
        elif var == 'temp_grads':  ## 梯度
            model_dic = self.get_data_temp(station)
            model_dic['obs'] = ds_obs['temp']

            model_dic = self.grads_data(model_dic)
            print(model_dic)
            pass

        elif var == 't_td_grads':  ## 梯度
            model_dic = self.get_data_t_td(station)
            model_dic['obs'] = ds_obs['temp'] - ds_obs['td']

            model_dic = self.grads_data(model_dic)
            print(model_dic)
            pass

        elif var == 'wind_grads':  ## 梯度
            model_dic = self.get_data_wind(station)
            model_dic['obs'] = ds_obs['wind_s']

            model_dic = self.grads_data(model_dic)
            print(model_dic)
            pass

        return model_dic


class Regrid():
    """对数据进行插值, 
    水平插值和垂直插值
    """
    def get_pressure_lev(self, station):
        #### 得到各层的pressure值
        ## 不同时刻各层气压值，差别可以忽略不计,
        # 后面还要对气压层进行插值, 这里不对它做过高精度要求
        path = '/mnt/zfm_18T/Asses_PBL/wrfout_data/'
        flnm_pressure = os.path.join(path, 'pressure_Jul_YSU_latlon')
        ds_pressure = xr.open_dataset(flnm_pressure)
        pr = ds_pressure.pressure
        prb = pr.sel(time='2016-07-01 13:00')
        lat = station['lat']
        lon = station['lon']
        # prc = prb.sel(lat=32.13, lon=92.5, method='nearest')
        prc = prb.sel(lat=lat, lon=lon, method='nearest')
        return prc

    def regrid(self, da_temp, station, pressure_lev):
        """对wrfout数据进行插值的
        需要水平插值和垂直插值两项
        Args:
            da_temp ([DataArray]): [需要插值的变量]

        Returns:
            [type]: [description]
        """

        #### 将bottom_top坐标换成气压坐标
        time_coord = da_temp.time.values
        lat_coord = da_temp.lat.values
        lon_coord = da_temp.lon.values
        # loc = {'lat':station['lat'], 'lon'}
        prc = self.get_pressure_lev(station)
        pressure_coord = prc.values
        # print(pressure_coord)

        da_temp_reset = da_temp.values
        # print(da_temp_reset)
        da_temp = xr.DataArray(
            da_temp_reset,
            coords=[time_coord, pressure_coord, lat_coord, lon_coord],
            dims=['time', 'pressure', 'lat', 'lon'])

        #### 对换成气压坐标的temp进行插值
        ## 水平插值
        # da_temp = da_temp.sel(lat=33.5, lon=97.5, method='nearest')
        da_temp = da_temp.sel(lat=station['lat'],
                              lon=station['lon'],
                              method='nearest')
        ## 垂直插值
        da_temp_return = da_temp.interp(pressure=pressure_lev)
        # print(da_temp_return)
        return da_temp_return


class Draw():
    def __init__(self, ):
        # self.da_var = da_var
        # self.name = name
        # self.station = station
        pass

    def draw_main(self, ):

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
        var = 'wind'
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
            level = np.arange(-19, 20, 3)
        elif var == 't_td':
            level = np.arange(0, 15, 1)
        elif var == 'wind':
            level = np.arange(0, 15, 1)

        CS = ax.contourf(x,
                         y,
                         val,
                         levels=level,
                         cmap=color_map,
                         extend='both')
        # cb = fig.colorbar(CS, orientation='horizontal', shrink=0.8, pad=0.14, fraction=0.14) # 这里的cs是画填色图返回的对象
        ## 设置标签大小
        ax.set_xticks(x[::2])  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax.set_xticklabels(x[::2].strftime('%m%d'))

        plt.gca().invert_yaxis()
        ax.set_ylim(600, 300)
        ax.set_title(title)
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
                            hspace=0.3)

        axes = [None] * 6  # 设置一个维度为8的空列表
        axes[0] = fig.add_subplot(grid[0, 0:1])
        axes[1] = fig.add_subplot(grid[0, 1:2])
        axes[2] = fig.add_subplot(grid[1, 0:1])
        axes[3] = fig.add_subplot(grid[1, 1:2])
        axes[4] = fig.add_subplot(grid[2, 0:1])

        model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']

        time_index1 = pd.date_range(start='20160702 00',
                                    end='20160731 00',
                                    freq='24H')
        time_index2 = pd.date_range(start='20160702 12',
                                    end='20160731 12',
                                    freq='24H')
        # time_index_dic = {'00':{'data':time_index1, 'name':'00'},
        #  '12':{'data':time_index2, 'name':'12'}}

        # time_index_dic = {'00':time_index1, '12':time_index2}
        time_index_dic = {'00': time_index1}
        # time_index_dic = {'12':time_index2}

        # for model in model_list:
        for key in time_index_dic:
            time_index = time_index_dic[key]
            CS = None
            for i in range(len(model_list)):
                title = str(station_name) + model_list[i] + str(key)
                CS = self.draw_contourf_single(var, axes[i],
                                               model_dic[model_list[i]], title,
                                               time_index)

            cb = fig.colorbar(CS,
                              orientation='horizontal',
                              shrink=0.8,
                              pad=0.14,
                              fraction=0.14)  # 这里的cs是画填色图返回的对象

            path = '/home/fengxiang/Project/Asses_PBL/Draw/UPAR/'
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
    gd = Get_data()
    aa = gd.get_data_main('temp_grads', station)
    # Dr = Draw()
    # Dr.draw_main()