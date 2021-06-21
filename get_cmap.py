
import matplotlib as mpl
import cmaps
import numpy as np
import meteva.base as mb
import matplotlib.pyplot as plt

def get_cmap_rain():
    ccc = cmaps.precip3_16lev_r
    ccc,clev = mb.def_cmap_clevs(mb.cmaps.rain_1h)
    colors = mpl.cm.get_cmap(ccc)
    col = colors(np.linspace(0, 1, 18))
    # print(col)
    cccc = mpl.colors.ListedColormap([
        col[0],
        col[1],
        col[2],
        col[3], 
        col[4], 
        # col[5], 
        col[6], 
        col[7], 
        col[8], 
        col[9], 
        col[10], 
        col[11], 
        # col[12], 
        col[13], 
        # col[14], 
        # (231 / 250, 177 / 250, 22 / 250),
        # col[4],
        # col[6],
        # '#85f485',
        # '#16c516',
        # 'white',
    ])
    cmap = cccc
    return cmap

def get_cmap_temp():
    # ccc = mb.cmaps.rain_1h
    # ccc,clev = mb.def_cmap_clevs(mb.cmaps.temp_2m, vmin=-30, vmax=35)
    # print(clev)

    # cccc, clev = mb.def_cmap_clevs(cmap=ccc, clevs=np.arange(300))

    ccc = cmaps.GMT_panoply
    # cmap = cmaps.temp_19lev
    # ccc = cmaps.circular_1
    # ccc = cmaps.precip3_16lev
    # print(ccc)
    colors = mpl.cm.get_cmap(ccc)
    col = colors(np.linspace(0, 1, 15))
    # # print(col)
    # col = colors(np.arange(0,1,11))
    # # # print(col)
    cccc = mpl.colors.ListedColormap([
        col[0],
        col[1],
        col[2],
        col[3],
        col[4],
        col[5],
        # col[6],
        # col[7],
        'white',
        col[8],
        col[9],
        col[10],
        col[11],
        col[12],
        col[13],
        col[14],
    ])

    # ## 在原来的基础上，多搞几个细分的颜色
    # cccc, clev = mb.def_cmap_clevs(cmap=cccc, clevs=np.arange(43))

    # # # print(len(cmap.values))
    cmap = cccc
    return cmap

def draw():
    fig, ax = plt.subplots(figsize=(6, 1))
    fig.subplots_adjust(bottom=0.5)
    # cmap = get_cmap_rain()
    cmap = get_cmap_temp()
    # norm = mpl.colors.Normalize(vmin=5, vmax=10)
    fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap),
             cax=ax, orientation='horizontal', label='Some Units')
    fig.savefig('test')

if __name__ == '__main__':
        
    # aa = get_cmap_rain()
    # print(aa)
    draw()