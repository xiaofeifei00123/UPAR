import xarray as xr
import pandas as pd
import numpy as np
import temp as tp

ta = pd.date_range('20160701 1200', '20160731 1200', freq='1H')


# ds = xr.Data
# ds = xr.Dataset.from_dataframe(ta)
# print(ds)
da = xr.DataArray(ta.values)

# ttt = aa.to_series()
ttt = da.values
# print(ttt)
# print(type(ttt))

# ta = da.values
# # ta = pd.Series(ta)  # 变成numpy


# # t1 = pd.to_datetime('20160701 1200').values
t1 = pd.date_range('20160701 00', '20160701 1200', freq='12H')
# t2 = pd.to_datetime('20160731 1200')


cc = np.setdiff1d(ttt, t1.values)
print(cc)



# # print(ta)
# # print(ta.values)

# for i in t1:
#     # print(i)
#     if i in ta:
#         # print(i)
#         # ta = ta.drop(i)
#         dd = np.delete(ta, i)
# # print(ta)
# print(dd)



# print(tb)
# tc = ta.where(ta!=tb, np.nan)
# if t1 in ta:
#     aa = ta.drop(t1)
# print(aa)






# xlabel = aa.dt.strftime('%d').values

# print(ta)
# print(tc.dropna())







# gd = tp.Get_data()

# path = '/mnt/zfm_18T/Asses_PBL/wrfout_data/pressure_Jul_TEMF_latlon'
# path = '/mnt/zfm_18T/Asses_PBL/wrfout_d02_2016-07-07_03:00:00'
# ds = xr.open_dataset(path)
# pass
# print(ds.T.values)


# da = np.arange(1,200,5)
# cor = np.arange(0.1, 20, 0.5)

# da = xr.DataArray(da, coords=[cor,], dims=['cor'])

# print(da.sel(cor=0.9, method='nearest'))
# print(da.interp(cor=[1,2,3]))