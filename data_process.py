#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
准备不同垂直坐标下:
1. 气压
2. 离地气压
3. 离地高度
多时次：
最低时次循环(观测6h, 模式1h)
的temp, td, wind;

包括:垂直插值+水平插值
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


class GetObs():
    def __init__(self, station) -> None:
        self.station = station

    def add_sureface(self, df):
        """预处理GPS探空数据
        主要是增加地面层的高度值

        Args:
            df ([type]): [description]
        """
        pass
        df = df.where(df < 9999, np.nan)  # 将缺省值赋值为NaN
        df = df.sort_values('pressure', ascending=False)  # 按照某一列排序
        max_pressure = df['pressure'].max()
        # df = df.where(df[''] < 9999, np.nan)  # 将缺省值赋值为NaN
        # if df['pressure'] == max_pressure
        df.loc[df['pressure'] == max_pressure,
               'height'] = self.station['height']
        # print(df)
        return df

    def read_single(self, flnm):
        pass
        """这里出现的这些问题，都是由于自己对pandas库不熟悉所导致的，
        我希望自己的工作，不仅仅是打工者的角色, 要做到更多更好一点才好
        """
        # pressure_level = np.arange(600, 100, -25)
        pressure_level = np.arange(600, 100, -5)
        # pressure_level = np.arange(0, 601, 25)
        col_names = [
            'pressure', 'height', 'temp', 'td', 't_td', 'wind_d', 'wind_s'
        ]
        ## 按列数据, 要哪几列的数据, 再命名
        df = pd.read_table(
            flnm,
            sep=' ',
            #  skiprows=0,
            usecols=[26, 27, 29, 30, 31, 32, 33],
            names=col_names,
        )

        df = self.add_sureface(df)

        df1 = df.dropna(axis=0, subset=['pressure'])  # 将含有缺省值的行删掉

        da_list = []
        var_list = ['temp', 'td', 't_td', 'wind_s', 'height']
        for var in var_list:
            pass
            df2 = df1  # 每一次重新传入这个df1,尽可能多的保留数据
            df2 = df2.dropna(axis=0, subset=[var])

            # 坐标处理
            df2.drop_duplicates('pressure', 'first', inplace=True)

            # 改变pressure坐标为相对的还是绝对的
            # df2['pressure'] = df2['pressure'].max() - df2['pressure']

            df2.set_index(['pressure'], inplace=True)  # 将pressure这一列设为index
            df3 = df2.sort_values('pressure')  # 按照某一列排序

            if var in ['td', 't_td']:
                if df3['t_td'].size == 0:
                    da = xr.DataArray([np.nan, np.nan],
                                      coords=[[0, 1]],
                                      dims='pressure')
                    df3['t_td'] = da.to_series()
                elif df3['td'].size == 0:
                    da = xr.DataArray([np.nan, np.nan],
                                      coords=[[0, 1]],
                                      dims='pressure')
                    df3['td'] = da.to_series()

            da = xr.DataArray.from_series(df3[var])
            dda = da.interp(pressure=pressure_level)
            # dda = dda.dropna(dim='pressure')
            # print(dda)
            da_list.append(dda)
        da_return = xr.concat(da_list, dim=var_list)
        # da_return = da_return.dropna(dim='pressure')
        # print(da_return.sel(concat_dim='height'))
        return da_return

    def read_obs(self, ):
        number = self.station['number']
        path1 = "/mnt/zfm_18T/Asses_PBL/GPS_Upar_2016/SCEX_TIPEX3_UPAR_GPS_MUL_"
        path = path1 + str(number) + "-201607"

        aa = os.listdir(path)  # 文件名列表
        aa.sort()  # 排一下顺序，这个是对列表本身进行操作

        ds_time = []  # 每个时次变量的列表
        ttt = []  # 时间序列列表
        for flnm in aa:

            fl_time = flnm[-14:-4]
            tt = pd.to_datetime(fl_time, format='%Y%m%d%H')
            ttt.append(tt)
            # 这时间是不规则的
            flnm = os.path.join(path, flnm)

            aa = self.read_single(flnm)
            ds_time.append(aa)

        ds = xr.concat(ds_time, dim='time')
        ds.coords['time'] = ttt
        return ds


class GetWrfout():
    """
    读取，处理数据
    """
    def __init__(self, station):
        self.station = station

    def get_data_single_once(
        self,
        var,
        flnm_var,
    ):
        '''单个变量，单个试验的数据读取'''
        pass
        # pressure_level = [
        #     600, 575, 550, 525, 500, 450, 400, 350, 300, 250, 200, 150, 100
        # ]
        # pressure_level = np.arange(0, 601, 25)
        # pressure_level = np.arange(600, 100, -25)
        pressure_level = np.arange(600, 100, -5)

        ds_var = xr.open_dataset(flnm_var)
        da_var = ds_var[var]

        # Re = Regrid()
        da_var = self.regrid_wrfout(da_var, pressure_level)
        return da_var
        # model_dic[model] = da_var
        # return model_dic

    def get_data_t_td(self, ):
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
            temp = self.get_data_single_once(var1, flnm_var1)
            td = self.get_data_single_once(var2, flnm_var2)
            model_dic[model] = temp - td
        return model_dic

    def get_data_wind(self, ):
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
            U = self.get_data_single_once(var1, flnm_var1)
            V = self.get_data_single_once(var2, flnm_var2)
            model_dic[model] = xr.ufuncs.sqrt(U**2 + V**2)
        return model_dic

    def get_data_temp(self, ):
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
            temp = self.get_data_single_once(var, flnm_var)
            if model == 'TEMF':
                print(temp)

            model_dic[model] = temp
        # print(model_dic['TEMF'].values)
        return model_dic

    def get_data_td(self, ):
        pass
        model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        var = 'td'  # 先就处理温度, 后面的数据是需要处理得到的
        # var2 = 'V'  # 先就处理温度, 后面的数据是需要处理得到的
        month = 'Jul'
        model_dic = {}
        for model in model_list:
            pass
            path = '/mnt/zfm_18T/Asses_PBL/wrfout_data/'
            file_name = str(var) + "_" + str(month) + \
                "_" + str(model) + "_latlon"
            flnm_var = os.path.join(path, file_name)
            temp = self.get_data_single_once(var, flnm_var)
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
            # 第0层是600hPa
            for i in range(len(pressure)):
                if i == 0:
                    # 边界上使用前差或后差
                    sr = (da[:, i + 1] - da[:, i]) / \
                        (pressure[i + 1] - pressure[i])
                    ser.append(sr)
                # 中间的使用中央差分
                elif i > 0 and i < len(pressure) - 1:
                    sr = (da[:, i + 1] - da[:, i - 1]) / \
                        (pressure[i] - pressure[i - 1]) / 2

                    ser.append(sr)
                elif i == len(pressure) - 1:
                    pass
                    sr = (da[:, i] - da[:, i - 1]) / \
                        (pressure[i] - pressure[i - 1])
                    ser.append(sr)
            da = xr.concat(ser, 'pressure')
            da.coords['pressure'] = pressure

            da = da.transpose()

            dic_return[model] = da
            # print(dic_return)
        return dic_return

    def regrid_wrfout(self, da_temp, pressure_lev):
        """对数据进行插值, 
        水平插值和垂直插值
        这是一个需要调用的类，也可以看成一个方法
        """

        # da_temp = da_temp

        def get_pressure_lev():
            # 得到各层的pressure值
            # 不同时刻各层气压值，差别可以忽略不计,
            # 后面还要对气压层进行插值, 这里不对它做过高精度要求
            path = '/mnt/zfm_18T/Asses_PBL/wrfout_data/'
            # 不同的方案和时间，在同一站点，各层气压值相差小于1度
            # 故不作分开考虑
            flnm_pressure = os.path.join(path, 'pressure_Jul_YSU_latlon')
            ds_pressure = xr.open_dataset(flnm_pressure)
            pr = ds_pressure.pressure
            prb = pr.sel(time='2016-07-01 13:00')
            lat = self.station['lat']
            lon = self.station['lon']
            # prc = prb.sel(lat=32.13, lon=92.5, method='nearest')
            prc = prb.sel(lat=lat, lon=lon, method='nearest')
            # prc = prc[0] - prc  # 离地气压高度
            return prc

        def get_height_lev():
            # 得到各层的pressure值
            # 不同时刻各层气压值，差别可以忽略不计,
            # 后面还要对气压层进行插值, 这里不对它做过高精度要求
            path = '/mnt/zfm_18T/Asses_PBL/wrfout_data/'
            # 不同的方案和时间，在同一站点，各层气压值相差小于1度
            # 故不作分开考虑
            flnm_pressure = os.path.join(path, 'z_Jul_YSU_latlon')
            ds = xr.open_dataset(flnm_pressure)
            pr = ds.z
            prb = pr.sel(time='2016-07-01 13:00')
            lat = self.station['lat']
            lon = self.station['lon']
            # prc = prb.sel(lat=32.13, lon=92.5, method='nearest')
            prc = prb.sel(lat=lat, lon=lon, method='nearest')
            # print(self.station['height'])
            # prc = prc - self.station['height']
            # prc = prc[0] - prc  # 离地气压高度
            return prc

        def regrid():
            """对wrfout数据进行插值的
            需要水平插值和垂直插值两项
            Args:
                da_temp (DataArray): 需要插值的变量

            Returns:
                DataArray: 插值后的变量
            """
            # 将bottom_top坐标换成气压坐标
            time_coord = da_temp.time.values
            lat_coord = da_temp.lat.values
            lon_coord = da_temp.lon.values
            # loc = {'lat':station['lat'], 'lon'}
            prc = get_pressure_lev()
            pressure_coord = prc.values
            # print(pressure_coord)

            da_temp_reset = da_temp.values
            # print(da_temp_reset)
            da = xr.DataArray(
                da_temp_reset,
                coords=[time_coord, pressure_coord, lat_coord, lon_coord],
                dims=['time', 'pressure', 'lat', 'lon'])

            # 对换成气压坐标的temp进行插值
            # 水平插值
            # da_temp = da_temp.sel(lat=33.5, lon=97.5, method='nearest')
            da = da.sel(lat=self.station['lat'],
                        lon=self.station['lon'],
                        method='nearest')
            # 垂直插值
            da_return = da.interp(pressure=pressure_lev)
            # da_return = da_return.dropna(dim='pressure')
            # print(da_return)
            return da_return

        # regrid()

        # return get_pressure_lev()
        # return get_height_lev()
        return regrid()


class TransferData(GetObs, GetWrfout):
    pass

    ## 当子类赋值属性和父类一样时， 子类实例化之后
    ## 父类拥有和子类相同的属性值
    def __init__(self, station):
        self.station = station  # 公有属性的处理

    def transfer_data(self, var):
        # var = self.var
        # pass
        ds_obs = self.read_obs()  # 观测数据

        if var == 'temp':
            pass
            # model_dic = None
            model_dic = self.get_data_temp()
            model_dic['obs'] = ds_obs.sel(concat_dim='temp')
        elif var == 'td':
            model_dic = self.get_data_td()
            model_dic['obs'] = ds_obs.sel(concat_dim='td')
        elif var == 't_td':
            pass
            model_dic = self.get_data_t_td()
            model_dic['obs'] = ds_obs.sel(concat_dim='t_td')
        elif var == 'wind':
            model_dic = self.get_data_wind()
            model_dic['obs'] = ds_obs.sel(concat_dim='wind_s')
            pass
        elif var == 'temp_grads':  # 梯度
            model_dic = self.get_data_temp()
            # model_dic['obs'] = ds_obs['temp']
            model_dic['obs'] = ds_obs.sel(concat_dim='temp')

            model_dic = self.grads_data(model_dic)
            # print(model_dic)
            pass
        elif var == 't_td_grads':  # 梯度
            model_dic = self.get_data_t_td()
            # model_dic['obs'] = ds_obs['temp'] - ds_obs['td']
            model_dic['obs'] = ds_obs.sel(concat_dim='t_td')

            model_dic = self.grads_data(model_dic)
            # print(model_dic)
            pass

        elif var == 'wind_grads':  # 梯度
            model_dic = self.get_data_wind()
            # model_dic['obs'] = ds_obs['wind_s']
            model_dic['obs'] = ds_obs.sel(concat_dim='wind_s')
            model_dic = self.grads_data(model_dic)
            # print(model_dic)
            pass
        else:
            model_dic = None
            print("错误")

        return model_dic


if __name__ == '__main__':

    pass
    station_dic = {
        'GaiZe': {
            'lat': 32.3,
            'lon': 84.0,
            'name': 'GaiZe',
            'number': '55248',
            'height': 4400,
        },
        'ShenZha': {
            'lat': 30.9,
            'lon': 88.7,
            'name': 'ShenZha',
            'number': '55472',
            'height': 4672
        },
        'ShiQuanhe': {
            'lat': 32.4,
            'lon': 80.1,
            'name': 'ShiQuanhe',
            'number': '55228',
            'height': 4280
        },
    }
    station = station_dic['GaiZe']
    # station = station_dic['ShenZha']

    # ## 传数据
    # tr = TransferData(station)
    # aa = tr.transfer_data('temp')
    # # print(aa)

    # 测试wrfout插值
    gw = GetWrfout(station)
    # aa = gw.get_pressure_lev(station)
    # aa = gw.regrid_wrfout('temp', [1,2,3])
    gw.get_data_temp()
    # print(aa)
    # # print(aa.values)

    # # 测试obs插值
    # flnm = '/mnt/zfm_18T/Asses_PBL/GPS_Upar_2016/SCEX_TIPEX3_UPAR_GPS_MUL_55228-201607/upar_G_55228_2016070206.txt'
    # gb = GetObs(station)
    # bb = gb.read_single(flnm)
    # # aa = gb.read_obs(station)
    # # print(aa.sel(concat_dim='td'))
