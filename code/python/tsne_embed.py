from config import config
import pickle
from scipy.io import loadmat
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

import pandas as pd
from ggplot import *


data_dir, cache_dir, plot_save_dir = config()
lfp_file = data_dir.joinpath('lfp_sections.mat')
label_file = data_dir.joinpath('labels.mat')
tsne_out_file = cache_dir.joinpath('tsne_out_cache.pickle')
tsne_plot_file = plot_save_dir.joinpath('tsne_plot.png')

lfp_sections = loadmat(lfp_file.as_posix(), chars_as_strings=True, appendmat=False)['lfp_sections']
labels = loadmat(label_file.as_posix(), chars_as_strings=True, appendmat=False)['labels'][0]

if not tsne_out_file.exists():
    print('Running PCA')
    pca_mod = PCA(n_components=50)
    pca_reduced = pca_mod.fit_transform(lfp_sections)
    print('Running t-SNE')
    tsne_mod = TSNE()
    tsne_out = tsne_mod.fit_transform(pca_reduced)
    with tsne_out_file.open('wb') as f:
        pickle.dump(tsne_out, f)
else:
    with tsne_out_file.open('rb') as f:
        tsne_out = pickle.load(f)

all_data = [[x[0], x[1], labels[i][0], (i % 89)+1] for i,x in enumerate(tsne_out)]
df = pd.DataFrame(all_data , columns=['x', 'y', 'labels', 'electrode_num'])
df.to_csv(data_dir.resolve().parent.joinpath('R').joinpath('tsne.csv').as_posix())

embedded_plot = ggplot(df, aes('x', 'y', color='labels', label='electrode_num')) + \
                geom_point() + \
                geom_text() + \
                theme_bw()
embedded_plot.save(tsne_plot_file.as_posix())