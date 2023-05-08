import matplotlib.pyplot as plt
import numpy as np

import mne

def show_topo(data ,vmin=0,vmax=1,contours = 6 , cmap = 'jet'):
    # 后面的操作不接受一维变量
    #plt.figure(figsize = (10 , 20) , dpi = 3000)
    data = data.reshape(62 , 1)
    # 创建信息
    locs_info_path='SEED_GRAPH.loc'
    montage=mne.channels.read_custom_montage(locs_info_path)
    ch_names=montage.ch_names
    #print(ch_names)
    info = mne.create_info(ch_names = ch_names , sfreq = 1000 , ch_types = 'eeg')
    # 创建evoke对象
    evoked = mne.EvokedArray(data , info)
    # 设置电极
    evoked.set_montage(montage)
    #evoked.plot_sensors(show_names=True)
    # 显示拓扑图
    ax , countour = mne.viz.plot_topomap(evoked.data[: , 0] ,evoked.info , 
                                        ch_type='eeg', sensors=True,names=None,
                                        outlines='head',extrapolate='head',
                                        image_interp='cubic',
                                        vlim=[vmin,vmax],cmap=cmap,
    				                    contours = contours,
                                        show = True) 
    #plt.show()
    #plt.savefig("1.jpg",dpi=1000, bbox_inches='tight')
    return ax , countour

data=np.random.rand(62)
show_topo(data)

