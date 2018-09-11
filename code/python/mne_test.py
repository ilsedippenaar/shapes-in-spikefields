import numpy as np
import mne
import pandas as pd
import matplotlib.pyplot as plt

from scipy.io import loadmat
from pathlib import Path

# %%

data_dir = Path('../..', 'data/python').resolve()
cache_dir = Path('../..', 'cache/python').resolve()

tbl = loadmat(str(data_dir / 'tbl'))['tbl']

a = {k: [x[0] for x in tbl[k]] for k in tbl.dtype.names}
df = pd.DataFrame(a)
lfps = np.stack([x.T for x in df.lfps]) # n_epochs, n_channels, n_samples

info = mne.create_info(ch_names=[f'e{i}' for i in range(53)],
                       ch_types='ecog',
                       sfreq=1000,
                       line_freq=60)
events = np.stack([i, 1, df['result'][i] == 'true_positive'] for i in range(len(df)))
epochs = mne.EpochsArray(data=lfps[:,:,0:256], info=info, tmin=-0.256,
                         events=events, event_id={'success': 1, 'failure': 0})

_ = epochs['failure'].average().plot()

con = mne.connectivity.spectral_connectivity(epochs, method='ppc', fmax=100, mt_bandwidth=10,
                                             indices=(np.zeros(52), np.arange(1,53)), verbose=False)

plt.plot(con[1], con[0].T[:,0])
plt.show()