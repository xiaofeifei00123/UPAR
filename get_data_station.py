#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:

-----------------------------------------
Time             :2021/06/21 18:02:34
Author           :Forxd
Version          :1.0
'''


import xarray as xr
import pandas as pd
import numpy as np
import os
import datetime
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




if __name__ == "__main__":
    pass

    station = {
                'lat': 32.3,
                'lon': 84.0,
                'name': 'GaiZe',
                'number': '55248'
            }

    gb = Get_obs()
    ds = gb.read_obs(station)
    # print(ds.time.sel(time='2016-07* 12'))
    print(ds.time.sel(time=datetime.time(12)))
    # path = "/mnt/zfm_18T/Asses_PBL/GPS_Upar_2016/SCEX_TIPEX3_UPAR_GPS_MUL_55228-201607/"
    # aa = os.listdir(path)
    # # print(type(aa))
    # aa.sort()

    # bb = []
    # ttt = []
    # for flnm in aa:

    #     fl_time = flnm[-14:-4]
    #     tt = pd.to_datetime(fl_time, format='%Y%m%d%H')
    #     ttt.append(tt)
    #     ## 这数是不规则的
    #     # print(flnm)
    #     # print(tt)
    #     flnm = os.path.join(path, flnm)

    #     aa = read_single(flnm)
    #     bb.append(aa)
    

    # dd = xr.concat(bb,dim='time')
    # dd.coords['time'] = ttt
    # print(dd['temp'])



