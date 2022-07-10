import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm
from matplotlib.ticker import LinearLocator

import seaborn as sns
import os.path
import json

import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("file_path", type=Path)
# parser.add_argument('--angle_abbrv', type=str, required=False)

args = parser.parse_args()

assert args.file_path.exists() , "The given directory doesn't exist"

output_dir = str(args.file_path)

df = pd.read_csv(output_dir+'/marginals/trc_df.csv', index_col=False)


sns.set()
sns.set_context("talk", font_scale=3)
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{amssymb}')


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})


# Make data.
X = df['theta'].to_numpy()
Y = df['rho'].to_numpy()
Z = df['avgcost'].to_numpy()


gridshape = (10,10)  # tuple

# # restore X,Y shapes for plotting
X = X.reshape(gridshape)
Y = Y.reshape(gridshape)
Z = Z.reshape(gridshape)


# Plot the surface: grayscale
# surf = ax.plot_surface(X, Y, Z, cmap=cm.gray, linewidth=0, antialiased=False)

# Plot the surface: color scale
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                    linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=10)

# Customize the z axis.
ax.set_xlabel(r'$\theta$', fontsize=20)
ax.set_ylabel(r'$\rho$', fontsize=20)
ax.zaxis.set_rotate_label(False) 
ax.set_zlabel(r'Average final cost', fontsize=20, rotation = 90)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

ax.set_zlim(np.min(Z), np.max(Z))
# change fontsize
for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(16)


# ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
# ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar, which maps values to colors.
# fig.colorbar(surf, shrink=0.5, aspect=5)

plt.savefig('viz_2d_avg_costs.png',format='png', bbox_inches='tight')
plt.show()
