from distutils import errors
import mfem.ser as mfem
from math import sqrt
import numpy as np
import pandas as pd
from pylab import *
import seaborn as sns
from seaborn.palettes import color_palette
from prob_envs.Poisson import Poisson

def SaveErrorsToFile(env, df_ErrorHistory, error_file_name):
    num_dofs = env.fespace.GetTrueVSize()
    df_tmp = pd.DataFrame({str(num_dofs):env.errors})
    df_ErrorHistory = pd.concat([df_ErrorHistory,df_tmp], axis=1)
    df_ErrorHistory.to_csv(error_file_name, index=False)
    return df_ErrorHistory

def GetElementVertices(env, k):
    Tr = env.mesh.GetElementTransformation(k)
    physical_pts = np.zeros((4,2))
    reference_pt = mfem.IntegrationPoint()
    for i in range(2):
        for j in range(2):
            reference_pt.Set(float(i),float(j),0.0,0.0)
            physical_pts[i+2*j,:] = Tr.Transform(reference_pt)
    return physical_pts

def SaveDistsToFile(env, df_dist, dist_file_name):
    dist = []
    for i in range(env.mesh.GetNE()):
        vertex_coords = GetElementVertices(env, i)
        min_dist = 1e8
        for j in range(4):
            x = vertex_coords[j,0]
            y = vertex_coords[j,1]
            min_dist = min(min_dist, sqrt(x**2 + y**2))
        dist.append(min_dist)
    df_tmp = pd.DataFrame(dist)
    df_dist = pd.concat([df_dist,df_tmp], axis=1)
    df_dist.to_csv(dist_file_name, index=False)
    return df_dist

sns.set()
sns.set_context("paper", font_scale=1.5)
sns.set_style("ticks")

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{amssymb}')

#### PARAMETERS
num_refs = 6
num_unif_ref = 1
order = 1
recompute = False
save_fig = True
theta = 0.0
####

error_file_name = "./ReleaseForPaper/out/H10errors_LShaped.csv"
dist_file_name = "./ReleaseForPaper/out/Dists_LShaped.csv"
prob_config = {
    'problem_type'      : 'lshaped',
    'mesh_name'         : 'l-shape-benchmark.mesh',
    'num_unif_ref'      : num_unif_ref,
    'order'             : order,
}

if recompute:
    env = Poisson(**prob_config)
    env.reset()
    df_ErrorHistory = pd.DataFrame()
    df_DistHistory = pd.DataFrame()
    for _ in range(num_refs):
        env.step(theta)
        df_ErrorHistory = SaveErrorsToFile(env, df_ErrorHistory, error_file_name)
        df_DistHistory = SaveDistsToFile(env, df_DistHistory, dist_file_name)
        print(np.sum(env.errors**2))

dofs = []
df_zeta = pd.read_csv(error_file_name)
df_dist = pd.read_csv(dist_file_name)
print(df_dist.iloc[:,0])

for i, col in enumerate(df_zeta.columns):
      dofs.append(float(col))
      num_dofs = float(col)
      num_non_zeros = len(df_zeta[col]) - df_zeta[col].isna().sum()
      df_zeta[col] *= num_non_zeros**(1/2) * num_dofs**(order/2)
      df_zeta.rename(columns={col:str(i)}, inplace=True)



# for dists, errors in zip(df_dist.column,df_zeta.columns):
for i in range(num_refs-1):
    df_tmp1 = pd.DataFrame({"dist":np.sqrt(df_dist.iloc[:,i])})
    # df_tmp1 = pd.DataFrame({"dist":df_dist.iloc[:,i]})
    df_tmp2 = pd.DataFrame({"error":np.log(df_zeta.iloc[:,i])})
    # df_tmp2 = pd.DataFrame({"error":np.log(1+df_zeta.iloc[:,i])})
    # df_tmp2 = pd.DataFrame({"error":df_zeta.iloc[:,i]})
    df = pd.concat([df_tmp1,df_tmp2], axis=1)
    ax = sns.regplot(data=df, x="dist", y="error")
    # ax = sns.regplot(data=df[(df["error"] >= theta*df["error"].max())], x="dist", y="error")
# ax = sns.scatterplot(x=df_zeta, y=df_dist, palette="RdBu_r")
ax.set_ylabel(r'$\{\overline{\eta}_T \colon T\in\mathcal{T}_k\}$')
ax.set_xlabel(r'Distance from singularity')
ax.legend()

## Convert x tick labels to LaTeX
# xticklabels = ax.get_xticklabels()
# for i in range(len(xticklabels)):
#       label = xticklabels[i].get_text()
#       xticklabels[i].set_text('$\\mathdefault{%i}$' % int(label))
# ax.set_xticklabels(xticklabels)

if save_fig:
    plt.savefig('./ReleaseForPaper/figures/Correlation.pdf')
plt.show()