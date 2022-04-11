import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import numpy as np
import pandas as pd
import csv


def fmtcols(mylist, cols):
    lines = (", ".join(mylist[i:i+cols]) for i in list(range(0,len(mylist),cols)))
    return '\n'.join(lines)


data = pd.read_excel('./data.xlsx', 
                    sheet_name='brut',
                    index_col=0,
                    usecols=['jangerpos', 'R2_Teratome', 'A3_Prognome', 'N1_Melanome', 'SYMBOL', 'Tier', 'Consequence'])

data.rename(columns={"R2_Teratome": "Mature Teratoma", "A3_Prognome": "Melanocytic Neuroectodermal Transformation (MNT)", "N1_Melanome": "Metastatic Anaplastic MNT"}, inplace=True)


""" Remove non-coding mutations """
data = data.loc[["intergenic" not in c for c in data.loc[:,'Consequence']],:]
data = data.loc[["intron" not in c for c in data.loc[:,'Consequence']],:]
data = data.loc[["UTR" not in c for c in data.loc[:,'Consequence']],:]
data = data.loc[["synonymous" not in c for c in data.loc[:,'Consequence']],:]
data = data.loc[["stream" not in c for c in data.loc[:,'Consequence']],:]


""" Define sets of mutations that are present for each condition """
A = set(data.loc[data.notna().loc[:,"Mature Teratoma"],:].index)
B = set(data.loc[data.notna().loc[:,"Melanocytic Neuroectodermal Transformation (MNT)"],:].index)
C = set(data.loc[data.notna().loc[:,"Metastatic Anaplastic MNT"],:].index)


""" Create venn diagram """
v = venn3([A,B,C], ("Mature Teratoma", "Melanocytic Neuroectodermal \nTransformation (MNT)", "Metastatic Anaplastic MNT"))



""" Remove initial text labels """
v.get_label_by_id('110').set_text('')
v.get_label_by_id('011').set_text('')
v.get_label_by_id('101').set_text('')
v.get_label_by_id('111').set_text('')



""" Replace labels by list of genes that are specific to each condition """
v.get_label_by_id('100').set_text(fmtcols(list(data.loc[A-B-C, 'SYMBOL'].loc[data.loc[A-B-C, 'Tier'].notna()]) + ["... (n={})".format(len(A-B-C))], 3))
v.get_label_by_id('100').set_weight("bold")
v.get_label_by_id('100').set_position((-0.33,0.1))
v.get_label_by_id('010').set_text(fmtcols(list(data.loc[B-A-C, 'SYMBOL'].loc[data.loc[B-A-C, 'Tier'].notna()]) + ["... (n={})".format(len(B-A-C))], 3))
v.get_label_by_id('010').set_weight("bold")
v.get_label_by_id('010').set_position((0.33,0.1))
v.get_label_by_id('001').set_text(fmtcols(list(data.loc[C-B-A, 'SYMBOL'].loc[data.loc[C-B-A, 'Tier'].notna()]) + ["... (n={})".format(len(C-B-A))], 2))
v.get_label_by_id('001').set_weight("bold")
v.get_label_by_id('001').set_position((0,-0.33))



abc_VAF1 = list(data.loc[A&B&C,'Mature Teratoma'])
abc_VAF2 = list(data.loc[A&B&C,'Melanocytic Neuroectodermal Transformation (MNT)'])
abc_VAF3 = list(data.loc[A&B&C,'Metastatic Anaplastic MNT'])


""" Split text labels for changing fonts, style and color for each Tier """
abc_zip = list(zip(data.loc[A&B&C,'Tier'], data.loc[A&B&C,'SYMBOL'], abc_VAF1, abc_VAF2, abc_VAF3))
abc = ['' if (t==1 or t==2) else '{} ({}%)'.format(n, round(max(v1,v2,v3),1)) for (t, n, v1, v2, v3) in abc_zip]
abc1 = ['{} ({}%)'.format(n, round(max(v1,v2,v3),1)) if t==1 else '' for (t, n, v1, v2, v3) in abc_zip]
abc2 = ['{} ({}%)'.format(n, round(max(v1,v2,v3),1)) if t==2 else '' for (t, n, v1, v2, v3) in abc_zip]

abc_xy = (0, 150)

plt.annotate('\n'.join(abc), xy=v.get_label_by_id('111').get_position() +
             np.array([0, 0]), xytext=abc_xy, ha='center',
             textcoords='offset points', 
             bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
             arrowprops=dict(arrowstyle='->',              
                             connectionstyle='arc',color='black'))
plt.annotate('\n'.join(abc2), xy=v.get_label_by_id('111').get_position() +
             np.array([0, 0]), xytext=abc_xy, ha='center', weight='bold',
             textcoords='offset points', 
             )                             
plt.annotate('\n'.join(abc1), xy=v.get_label_by_id('111').get_position() +
             np.array([0, 0]), xytext=abc_xy, ha='center', weight='bold', color='red',
             textcoords='offset points', 
             )     


ab_VAF1 = list(data.loc[A&B-C,'Mature Teratoma'])
ab_VAF2 = list(data.loc[A&B-C,'Melanocytic Neuroectodermal Transformation (MNT)'])
ab_VAF3 = list(data.loc[A&B-C,'Metastatic Anaplastic MNT'])

ab_zip = list(zip(data.loc[A&B-C,'Tier'], data.loc[A&B-C,'SYMBOL'], ab_VAF1, ab_VAF2, ab_VAF3))
ab = ['' if (t==1 or t==2) else '{} ({}%)'.format(n, round(max(v1,v2),1)) for (t, n, v1, v2, v3) in ab_zip]
ab1 = ['{} ({}%)'.format(n, round(max(v1,v2),1)) if t==1 else '' for (t, n, v1, v2, v3) in ab_zip]
ab2 = ['{} ({}%)'.format(n, round(max(v1,v2),1)) if t==2 else '' for (t, n, v1, v2, v3) in ab_zip]


ab_xy = (250, -200)
plt.annotate('\n'.join(ab), xy=v.get_label_by_id('110').get_position() +
             np.array([0, 0]), xytext=ab_xy, ha='center',
             textcoords='offset points', 
             bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.1),
             arrowprops=dict(arrowstyle='->',              
                             connectionstyle='arc',color='black'))
plt.annotate('\n'.join(ab2), xy=v.get_label_by_id('110').get_position() +
             np.array([0, 0]), xytext=ab_xy, ha='center',
             textcoords='offset points',  weight='bold',
             )
plt.annotate('\n'.join(ab1), xy=v.get_label_by_id('110').get_position() +
             np.array([0, 0]), xytext=ab_xy, ha='center',
             textcoords='offset points',  weight='bold', color='red',
             )                                                          

ac_VAF1 = list(data.loc[A&C-B,'Mature Teratoma'])
ac_VAF2 = list(data.loc[A&C-B,'Melanocytic Neuroectodermal Transformation (MNT)'])
ac_VAF3 = list(data.loc[A&C-B,'Metastatic Anaplastic MNT'])

ac_zip = list(zip(data.loc[A&C-B,'Tier'], data.loc[A&C-B,'SYMBOL'], ac_VAF1, ac_VAF2, ac_VAF3))
ac = ['' if (t==1 or t==2) else '{} ({}%)'.format(n, round(max(v1,v3),1)) for (t, n, v1, v2, v3) in ac_zip]
ac1 = ['{} ({}%)'.format(n, round(max(v1,v3),1)) if t==1 else '' for (t, n, v1, v2, v3) in ac_zip]
ac2 = ['{} ({}%)'.format(n, round(max(v1,v3),1)) if t==2 else '' for (t, n, v1, v2, v3) in ac_zip]

ac_xy = (-150, -100)
plt.annotate('\n'.join(ac), xy=v.get_label_by_id('101').get_position() +
             np.array([0, 0]), xytext=ac_xy, ha='center',
             textcoords='offset points', 
             bbox=dict(boxstyle='round,pad=0.5', fc='purple', alpha=0.1),
             arrowprops=dict(arrowstyle='->',              
                             connectionstyle='arc',color='black'))
plt.annotate('\n'.join(ac2), xy=v.get_label_by_id('101').get_position() +
             np.array([0, 0]), xytext=ac_xy, ha='center',
             textcoords='offset points', weight = 'bold',
             #bbox=dict(boxstyle='round,pad=0.5', fc='purple', alpha=0.1),
             )
plt.annotate('\n'.join(ac1), xy=v.get_label_by_id('101').get_position() +
             np.array([0, 0]), xytext=ac_xy, ha='center',
             textcoords='offset points', weight = 'bold', color = 'red',
             #bbox=dict(boxstyle='round,pad=0.5', fc='purple', alpha=0.1),
             )                                                          

bc_VAF1 = list(data.loc[B&C-A,'Mature Teratoma'])
bc_VAF2 = list(data.loc[B&C-A,'Melanocytic Neuroectodermal Transformation (MNT)'])
bc_VAF3 = list(data.loc[B&C-A,'Metastatic Anaplastic MNT'])

bc_zip = list(zip(data.loc[B&C-A,'Tier'], data.loc[B&C-A,'SYMBOL'], bc_VAF1, bc_VAF2, bc_VAF3))
bc = ['' if (t==1 or t==2) else '{} ({}%)'.format(n, round(max(v2,v3),1)) for (t, n, v1, v2, v3) in bc_zip]
bc1 = ['{} ({}%)'.format(n, round(max(v2,v3),1)) if t==1 else '' for (t, n, v1, v2, v3) in bc_zip]
bc2 = ['{} ({}%)'.format(n, round(max(v2,v3),1)) if t==2 else '' for (t, n, v1, v2, v3) in bc_zip]

bc_xy = (100, -100)
plt.annotate('\n'.join(bc), xy=v.get_label_by_id('011').get_position() +
             np.array([0, 0]), xytext=bc_xy, ha='center',
             textcoords='offset points', 
             bbox=dict(boxstyle='round,pad=0.5', fc='blue', alpha=0.1),
             arrowprops=dict(arrowstyle='->',              
                             connectionstyle='arc',color='black'))     
plt.annotate('\n'.join(bc2), xy=v.get_label_by_id('011').get_position() +
             np.array([0, 0]), xytext=bc_xy, ha='center',
             textcoords='offset points', weight = 'bold',
             #bbox=dict(boxstyle='round,pad=0.5', fc='blue', alpha=0.1),
             )     
plt.annotate('\n'.join(bc1), xy=v.get_label_by_id('011').get_position() +
             np.array([0, 0]), xytext=bc_xy, ha='center',
             textcoords='offset points', weight = 'bold', color = 'red',
             #bbox=dict(boxstyle='round,pad=0.5', fc='blue', alpha=0.1),
             )                                                 



a_VAF1 = list(data.loc[A-B-C,'Mature Teratoma'])
a_VAF2 = list(data.loc[A-B-C,'Melanocytic Neuroectodermal Transformation (MNT)'])
a_VAF3 = list(data.loc[A-B-C,'Metastatic Anaplastic MNT'])

b_VAF1 = list(data.loc[B-A-C,'Mature Teratoma'])
b_VAF2 = list(data.loc[B-A-C,'Melanocytic Neuroectodermal Transformation (MNT)'])
b_VAF3 = list(data.loc[B-A-C,'Metastatic Anaplastic MNT'])

c_VAF1 = list(data.loc[C-A-B,'Mature Teratoma'])
c_VAF2 = list(data.loc[C-A-B,'Melanocytic Neuroectodermal Transformation (MNT)'])
c_VAF3 = list(data.loc[C-A-B,'Metastatic Anaplastic MNT'])




""" Define scores """

scores = []                                      

scores.extend([(g, 1)for (_,g,x,y,z) in list(zip(data.loc[A-B-C,'Tier'], data.loc[A-B-C,'SYMBOL'], a_VAF1, a_VAF2, a_VAF3))])
scores.extend([(g, 2)for (_,g,x,y,z) in list(zip(data.loc[B-A-C,'Tier'], data.loc[B-A-C,'SYMBOL'], b_VAF1, b_VAF2, b_VAF3))])
scores.extend([(g, 3)for (_,g,x,y,z) in list(zip(data.loc[C-B-A,'Tier'], data.loc[C-B-A,'SYMBOL'], c_VAF1, c_VAF2, c_VAF3))])
scores.extend([(g, 3)for (_,g,x,y,z) in ab_zip])
scores.extend([(g, 4)for (_,g,x,y,z) in ac_zip])
scores.extend([(g, 5)for (_,g,x,y,z) in bc_zip])
scores.extend([(g, 6)for (_,g,x,y,z) in abc_zip])



with open('venn_scores.rnk','w', newline='') as out:
    csv_out=csv.writer(out, delimiter='\t')
    for row in scores:
        csv_out.writerow(row)


plt.show()