import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['figure.figsize'] = (20,25)
plt.style.use('fivethirtyeight')

starting_data = pd.read_csv('starting_hits.txt', sep=" ")
dead_data = pd.read_csv('output_hits.txt', sep=" ")

methods = ["CNN", "BETWEEN", "NEAR"]
reco_data = {}
# true_positives = {}
# false_positives = {}
tp_list = []
fp_list = []

for method in methods:
    reco_data[method] = pd.read_csv('reco_hits_'+method+'.txt', sep=" ")
    reco_data[method] = pd.merge(reco_data[method], starting_data, on=['ID', 'SIDE', 'LAYER', 'ROW'], how='left', indicator='IS_PRESENT')
    reco_data[method]['IS_PRESENT'] = np.where(reco_data[method].IS_PRESENT == 'both', True, False)
    false_positives, true_positives = reco_data[method].IS_PRESENT.value_counts().sort_index().tolist()
    #     false_positives[method], true_positives[method] = reco_data[method].IS_PRESENT.value_counts().sort_index().tolist()
    tp_list.append(true_positives)
    fp_list.append(false_positives)

df_stats = pd.DataFrame({'True Positives': tp_list, 'False Positives': fp_list}, index=methods)
ax = df_stats.plot.bar(rot=0)
fig = ax.get_figure()
fig.savefig('Reco_stats.png')

fraction_missed = dead_data.count()[0]/starting_data.count()[0]
print("With no method, we get", fraction_missed, " of original hits")

gt_per_ev = starting_data.groupby('ID').count().SIDE.to_numpy()

survived_per_ev = dead_data.groupby('ID').count().SIDE.to_numpy()

missed_per_ev = gt_per_ev - survived_per_ev

tp_per_ev_list = []
fp_per_ev_list = []
avg_true_reco = []
avg_false_reco = []
for method in methods:
    tp = reco_data[method].groupby('ID').sum().reindex(dead_data.groupby('ID').count().index).IS_PRESENT.fillna(0).to_numpy()
    fp = (reco_data[method].groupby('ID').count().reindex(dead_data.groupby('ID').count().index).IS_PRESENT.fillna(0) - tp).to_numpy()
    tp_per_ev_list.append(tp)
    fp_per_ev_list.append(fp)
    avg_true_reco.append(np.average((survived_per_ev + tp)/gt_per_ev))
    avg_false_reco.append(np.average(fp/gt_per_ev))

avg_per_ev = dict(zip(methods, np.vstack((avg_true_reco, avg_false_reco)).T))

df_stats_per_ev = pd.DataFrame(avg_per_ev, index=['Avg TPR per Event', 'Avg FPR per Event'])
ax2 = df_stats_per_ev.plot.bar(rot=0)
fig2 = ax2.get_figure()
plt.axhline(y=fraction_missed, color='r', linestyle='--')
fig2.savefig('Reco_stats_per_ev.png')

df_stats_per_ev = pd.DataFrame(avg_per_ev, index=['Avg TPR per Event', 'Avg FPR per Event'])
ax2 = df_stats_per_ev.plot.bar(rot=0)
fig2 = ax2.get_figure()
plt.axhline(y=fraction_missed, color='r', linestyle='--')
fig2.savefig('Reco_stats_per_ev.png')

df_stats_per_ev2 = pd.DataFrame({'% True recovered hits': avg_true_reco, '% False recovered hits': avg_false_reco}, index=methods)
ax3 = df_stats_per_ev2.plot.bar(rot=0)
fig3 = ax3.get_figure()
plt.axhline(y=fraction_missed, color='r', linestyle='--')
fig3.savefig('Reco_stats_per_ev2.png')
