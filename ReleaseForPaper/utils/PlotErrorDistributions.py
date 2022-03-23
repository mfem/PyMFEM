from math import sqrt
import numpy as np
import pandas as pd
from pylab import *
import seaborn as sns
from seaborn.palettes import color_palette
from prob_envs.StationaryProblem import StationaryProblem

sns.set()
# sns.set_context("notebook")
# sns.set_context("paper")
sns.set_context("paper", font_scale=1.5)
# sns.set_context("talk")
sns.set_style("ticks")

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{amssymb}')

#### PARAMETERS
fig = 'a' # 'a' or 'b'
num_refs = 6
num_unif_ref = 1
recompute = True
save_fig = True
####

if fig == 'a':
      file_name = "./RLModule/out/H10errors_sinsin.csv"
      prob_config = {
            'estimator_type'    : 'exact',
            'problem_type'      : 'sinsin',
            'mesh_name'         : 'inline-quad_2by2.mesh',
            'num_unif_ref'      : num_unif_ref,
            'order'             : 1,
            'error_file_name'   : file_name
      }
else:
      file_name = "./RLModule/out/H10errors_LShaped.csv"
      prob_config = {
            'estimator_type'    : 'exact',
            'problem_type'      : 'lshaped',
            'mesh_name'         : 'l-shape-benchmark.mesh',
            'num_unif_ref'      : num_unif_ref,
            'order'             : 1,
            'error_file_name'   : file_name
      }

if recompute:
      env = StationaryProblem(**prob_config)
      env.reset(save_errors=True)
      for _ in range(num_refs):
            env.step(0.0)
            print(np.sum(env.errors**2))

dofs = []
df = pd.read_csv(file_name)
for i, col in enumerate(df.columns):
      dofs.append(float(col))
      num_dofs = float(col)
      num_non_zeros = len(df[col]) - df[col].isna().sum()
      # df[col] *= df[col]
      df[col] *= num_non_zeros
      # df[col] *= num_non_zeros
      # df[col] = - np.log(df[col]) / np.log(num_non_zeros**2)
      # df[col] = np.log(df[col]) / np.log(num_dofs)
      df.rename(columns={col:str(i)}, inplace=True)

dofs = np.array(dofs)
means = df.mean().to_numpy()
medians = df.median().to_numpy()
mins = df.min().to_numpy()
maxes = df.max().to_numpy()

proxy_df = pd.DataFrame(columns = df.columns)
for i, col in enumerate(proxy_df.columns):
      proxy_df[col] = [means[i]]

# ax = sns.stripplot(data=proxy_df, zorder=10,  color="white", linewidth=1, jitter=False, edgecolor="black")
ax = sns.boxenplot(data=df, width=.6,
                  # palette="gray"
                  palette="RdBu_r" ## !!
                  # palette="PuBu_r"
                  # palette="coolwarm"
                  # palette="Spectral"
                  )
if fig != 'a':
      ax.set_yscale('log')
ax.set_ylabel(r'Local element errors (normalized)')
ax.set_xlabel(r'Refinement')

## Convert x tick labels to LaTeX
xticklabels = ax.get_xticklabels()
for i in range(len(xticklabels)):
      label = xticklabels[i].get_text()
      xticklabels[i].set_text('$\\mathdefault{%i}$' % int(label))
ax.set_xticklabels(xticklabels)

print('Means   = ', means)
print('Medians = ', medians)
print('Mins    = ', mins)
print('Maxes   = ', maxes)
print('conv. rate of means   = ', 2 * np.log(means[-1]/means[-2])/np.log(dofs[-2]/dofs[-1]))
print('conv. rate of medians = ', 2 * np.log(medians[-1]/medians[-2])/np.log(dofs[-2]/dofs[-1]))
print('conv. rate of mins    = ', 2 * np.log(mins[-2]/mins[-3])/np.log(dofs[-3]/dofs[-2]))
print('conv. rate of maxes   = ', 2 * np.log(maxes[-1]/maxes[-2])/np.log(dofs[-2]/dofs[-1]))

if save_fig:
      if fig == 'a':
            plt.savefig('./RLModule/figures/fig1a_H10Errors_sinsin.pdf')
      else:
            plt.savefig('./RLModule/figures/fig1b_H10Errors_LShaped.pdf')
plt.show()

# print(np.log(env.global_errors[1][-1]/env.global_errors[1][-5])/np.log(env.global_errors[0][-5]/env.global_errors[0][-1]))
# print(np.log(env.global_errors[2][-1]/env.global_errors[2][-5])/np.log(env.global_errors[0][-5]/env.global_errors[0][-1]))
# print(np.log(env.global_errors[3][-1]/env.global_errors[3][-5])/np.log(env.global_errors[0][-5]/env.global_errors[0][-1]))


# plt.title('Errors')
# plt.loglog(env.global_errors[0], env.global_errors[1], 'r-', label='|| grad(u - u_h)||')
# plt.loglog(env.global_errors[0], env.global_errors[2], 'b-', label='|| y - grad(u_h)||')
# plt.loglog(env.global_errors[0], env.global_errors[3], 'k-', label='|| y_ZZ - grad(u_h)||')
# plt.legend()