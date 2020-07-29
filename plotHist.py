# Sicheng He hschsc@umich.edu

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib as matplotlib

# color
my_blue = '#4C72B0'
my_red = '#C54E52'
my_green = '#56A968' 
my_brown = '#b4943e'
my_purple = '#684c6b'
my_orange = '#cc5500'

# font
matplotlib.rcParams.update({'font.size': 20})
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

# area distribution
area_file_name = "./solution/old_opt_sol/areas.dat"
areas = np.loadtxt(area_file_name)
print(min(areas))
print(max(areas))

n_bins = 100

fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
axs.hist(np.log10(areas), bins=n_bins)

axs.set_xlabel(r'$\log_{10}(Area)$', fontsize=20)
axs.set_ylabel(r'number', fontsize=20, rotation=0)

axs.spines['right'].set_visible(False)
axs.spines['top'].set_visible(False)

fig.savefig("./solution/old_opt_sol/area.pdf", bbox_inches='tight')


stress_file_name_list = ["./solution/new_case_1/sol_stress.txt",
    "./solution/new_case_2/sol_stress.txt",
    "./solution/new_case_3/sol_stress.txt",
    "./solution/old_case_1/sol_stress.txt",
    "./solution/old_case_2/sol_stress.txt",
    "./solution/old_case_3/sol_stress.txt"]
buckling_file_name_list = ["./solution/new_case_1/sol_buckling.txt",
    "./solution/new_case_2/sol_buckling.txt",
    "./solution/new_case_3/sol_buckling.txt",
    "./solution/old_case_1/sol_buckling.txt",
    "./solution/old_case_2/sol_buckling.txt",
    "./solution/old_case_3/sol_buckling.txt"]

case_name_list = ["new_case_1",
"new_case_2",
"new_case_3",
"old_case_1",
"old_case_2",
"old_case_3",]

for i in range(len(case_name_list)):

    stress = np.loadtxt(stress_file_name_list[i])
    buckling = np.loadtxt(buckling_file_name_list[i])

    n_bins = 30

    # stress
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    N, bins, patches = axs.hist(stress, bins=n_bins)

    fracs = N / N.max()

    coolwarm_cmap = matplotlib.cm.get_cmap('coolwarm')
    norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)
    j = 0
    for thisfrac, thispatch in zip(fracs, patches):
        sigma = patches[j].xy[0] + patches[j]._width / 2
        color = coolwarm_cmap(norm(-sigma))
        thispatch.set_facecolor(color)

        j += 1

    axs.set_xlabel(r'$\sigma / \sigma_Y$', fontsize=20)
    axs.set_ylabel(r'number', fontsize=20, rotation=0)

    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)

    height = N.max() / 2
    axs.text(1.0, height, r'Yield stress', color = 'k', fontsize=20, alpha=0.3)
    axs.plot([-1, -1], [0, height], '-', color='k', alpha=0.3)
    axs.plot([1, 1], [0, height], '-', color='k', alpha=0.3)

    file_name = "./solution/" + case_name_list[i] + "/" + "stress_hist.pdf"
    fig.savefig(file_name, bbox_inches='tight')

    # buckling
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    N, bins, patches = axs.hist(buckling, bins=n_bins)

    fracs = N / N.max()

    coolwarm_cmap = matplotlib.cm.get_cmap('coolwarm')
    norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)
    j = 0
    for thisfrac, thispatch in zip(fracs, patches):
        bucklingVar = patches[j].xy[0] + patches[j]._width / 2
        color = coolwarm_cmap(norm(bucklingVar))
        thispatch.set_facecolor(color)

        j += 1

    axs.set_xlabel(r'$max(-\sigma / (\gamma A), 0)$', fontsize=20)
    axs.set_ylabel(r'number', fontsize=20, rotation=0)

    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)

    height = N.max() / 2
    axs.text(1.0, height, r'Buckling', color = 'k', fontsize=20, alpha=0.3)
    axs.plot([1, 1], [0, height], '-', color='k', alpha=0.3)

    file_name = "./solution/" + case_name_list[i] + "/" + "buckling_hist.pdf"
    fig.savefig(file_name, bbox_inches='tight')