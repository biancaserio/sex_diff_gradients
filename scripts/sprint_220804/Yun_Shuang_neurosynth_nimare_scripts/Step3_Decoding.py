
# step1_read dataset
from neurosynth.base.dataset import Dataset
dataset = Dataset('/data/pt_02542/social_networks/results/Decoding/Datasets/v3/database.txt')
dataset.add_features('/data/pt_02542/social_networks/results/Decoding/Datasets/v3/features.txt')
dataset.get_feature_names()
#dataset.save('/data/pt_02542/social_networks/results/Decoding/Datasets/v3/Dataset.pkl')


# step2_generating ROI mask in shell





# step3_decoding
from neurosynth.base.dataset import Dataset
from neurosynth.analysis import decode
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['svg.fonttype'] = 'none'
# rank topics according to volumeNO.
def getOrder(d, thr):
    dh = []
    for i in range(0,len(d)):
        di = d[i]
        dh.append(np.average(np.array(range(0,len(d[i]))) + 1, weights=di))
    heatmapOrder = np.argsort(dh)
    return heatmapOrder

# Import neurosynth database:
pickled_dataset = '/data/pt_02542/social_networks/results/Decoding/Datasets/v3/dataset.pkl'
dataset = Dataset.load(pickled_dataset)


# Analysis with 24 terms:
features = pd.read_csv('/data/pt_02542/social_networks/results/Decoding/Datasets/v3/topic50.txt', sep='\t', index_col=0)

#34 topics
topics_to_keep=[0,1,6,7,8,9,11,13,14,15,16,17,19,20,23,25,26,27,28,
29,30,32,33,35,37,38,40,41,42,44,45,46,47,49]

features = features.iloc[:, topics_to_keep]

features=features.rename(columns = lambda x: x.replace('LDA50_abstract_weight__', ''))

dataset.add_features(features, append=False)

# Gradient 1

decoder = decode.Decoder(dataset, method='roi')

# Set threshold:
thr = 0
vmin = 0
vmax = 4

tot = 5
data = decoder.decode([str('/mnt/hgfs/Share/neurosynth/G2/volume_%02d_%02d.nii.gz' % (i * tot, (i * tot) + tot)) for i in range(0,int(100/tot))])
#data = decoder.decode([str('/mnt/hgfs/Share/neurosynth/CPcAllmasks/volume_%02d_%02d.nii.gz' % (i * tot, (i * tot) + tot)) for i in range(0,int(100/tot))])

df = []
df = data.copy()
newnames = []
[newnames.append(('%s-%s' % (str(i * tot), str((i*tot) + tot)))) for i in range(0,len(df.columns))]
df.columns = newnames


dfnew = data
df[df<thr] = 0 
# 9 topics are noise, leaving 25 topics
test=dfnew[~(df==0).all(1)]
test[test<thr] = 0
heatmapOrder = getOrder(np.array(test), thr)

test.columns = newnames

sns.set(context="paper", font="sans-serif", font_scale=2)
f, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(30, 20), sharey=True)
plotData = test.reindex(test.index[heatmapOrder])
cax = sns.heatmap(plotData, linewidths=1, square=True, cmap='Greys', robust=False, 
            ax=ax1, vmin=vmin, vmax=vmax, mask=plotData == 0)

sns.axlabel('Percentile along gradient', 'NeuroSynth topics terms')
cbar = cax.collections[0].colorbar
cbar.set_label('z-stat', rotation=270)
cbar.set_ticks(ticks=[thr,vmax])
cbar.set_ticklabels(ticklabels=[thr,vmax])
cbar.outline.set_edgecolor('black')
cbar.outline.set_linewidth(0.5)

plt.draw()
plt.show()
f.savefig('/mnt/hgfs/Share/neurosynth/newest/33_32cognitive_G1ROI_thr0.svg', format='svg')

