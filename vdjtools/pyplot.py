import seaborn as sns, matplotlib.pyplot as plt
import math
import numpy as np 

# clonotypes count
input_path='/SGRNJ06/randd/USER/cjj/celedev/vdj_bulk/20230131demo/plot/vdjtools_result/'#输入路径
output_path='/SGRNJ06/randd/USER/cjj/celedev/vdj_bulk/20230131demo/plot/vdjtools_result/'#输出路径
sample_info='basicstats.txt'#记录有样本分组的样本信息
########################################################################################################################
list_pos=[]
for line in open(input_path+sample_info,'r'):
    info=line[:-1].split('\t')
    if info[0]!='ID':
        if str(info[1])=='Class1':
            list_pos.append(info[0])

list_class=[]
list_sample=[]
list_richness=[]
for line in open(output_path+'basicstats.txt','r'):
    info=line[:-1].split('\t')
    if info[0]!='sample_id':
        if info[0].split('.')[0] in list_pos:
            list_class.append('Class1')
        else:
            list_class.append('Class2')
        list_sample.append(info[0].split('.')[0])
        list_richness.append(float(info[3]))

# ax1 = sns.barplot(x=list_class,y=list_richness,capsize=.2,palette=(sns.xkcd_rgb["dark red"],sns.xkcd_rgb["marine blue"]))
# ax1_1 = sns.stripplot(x=list_class,y=list_richness,size=5,palette=(sns.xkcd_rgb["dark red"],sns.xkcd_rgb["marine blue"]))
# ax1.spines['top'].set_visible(False)
# ax1.spines['right'].set_visible(False)
#plt.savefig(output_path+'Richness_barplot_Class.pdf', dpi=300, bbox_inches="tight")
#plt.show()
list_sample = [i + "_clonotypes" for i in list_sample]
ax2 = sns.barplot(x=list_sample,y=list_richness,hue=list_sample,dodge=False)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
plt.legend(bbox_to_anchor=(1.05, 1),loc=2,borderaxespad=0)
plt.title("Clonotypes Count", fontsize=15)
plt.savefig(output_path+'Richness_barplot_Sample.pdf', dpi=300, bbox_inches="tight")
plt.show()


# CDR3 Clonotypes Length
input_path='/SGRNJ06/randd/USER/cjj/celedev/vdj_bulk/20230131demo/plot/vdjtools_result/'#输入路径
output_path='/SGRNJ06/randd/USER/cjj/celedev/vdj_bulk/20230131demo/plot/vdjtools_result/'#输出路径
sample_info='basicstats.txt'#记录有样本分组的样本信息
dict_sample_LengthCount={}
list_pos=[]
for line in open(input_path+sample_info,'r'):
    info=line[:-1].split('\t')
    id=info[0]
    if id == "sample_id":
        continue
    if id!='ID':
        if str(info[1])=='Class1':
            list_pos.append(id)
        read_file=open(input_path+id+'.txt','r')
        dict_length_count={}
        for line in read_file:
            if line[0]!='c':
                length=len(line[:-1].split('\t')[3])
                if length < 30:
                    if length not in dict_length_count.keys():
                        dict_length_count[length]=1
                    else:
                        dict_length_count[length]+=1
        dict_sample_LengthCount[id]=dict_length_count
        read_file.close()
list_length=[]
list_class=[]
list_sample=[]
list_count=[]
for sample in dict_sample_LengthCount.keys():
    dict_length_count=dict_sample_LengthCount[sample]
    for length in dict_length_count.keys():
        list_length.append(length)
        list_sample.append(sample)
        if sample in list_pos:
            list_class.append('Class1')
        else:
            list_class.append('Class2')
        list_count.append(dict_length_count[length])

# ax1 = sns.barplot(x=list_length,y=list_count,hue=list_class,errwidth=0.5,capsize=0.5,palette=(sns.xkcd_rgb["dark red"],sns.xkcd_rgb["marine blue"]))
# ax1_1 = sns.stripplot(x=list_length,y=list_count,hue=list_class,size=2,dodge=True,palette=(sns.xkcd_rgb["dark red"],sns.xkcd_rgb["marine blue"]))
# ax1.set_xticklabels(ax1.get_xticklabels(), rotation=-90,fontsize=6)
# ax1.spines['top'].set_visible(False)
# ax1.spines['right'].set_visible(False)
# #plt.savefig(output_path+'Length_barplot_Class.pdf', dpi=300, bbox_inches="tight")
# plt.show()
list_sample = [i.split('.')[0] + "_clonotypes" for i in list_sample]
ax2 = sns.barplot(x=list_length,y=list_count,hue=list_sample,errwidth=0.5,capsize=0.5)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=-0,fontsize=6)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
plt.title("CDR3 Clonotypes Length", fontsize=15)
plt.legend(bbox_to_anchor=(1.05, 1),loc=2,borderaxespad=0)
plt.savefig(output_path+'Length_barplot_Sample.pdf', dpi=300, bbox_inches="tight")
plt.show()


# CDR3 diversity / evenness / clonality

input_path='/SGRNJ06/randd/USER/cjj/celedev/vdj_bulk/20230131demo/plot/vdjtools_result/'#输入路径
output_path='/SGRNJ06/randd/USER/cjj/celedev/vdj_bulk/20230131demo/plot/vdjtools_result/'#输出路径
sample_info='basicstats.txt'#记录有样本分组的样本信息
########################################################################################################################
def Shannon_entropy(list_frequency):
    sum=0
    for frequency in list_frequency:
        sum+=frequency*math.log(frequency)
    H=-sum
    return H

def Pielou_evenness(H,richness):
    E=H/math.log(richness)
    return E

def Clonality(E):
    C=1-E
    return C

list_pos=[]
list_characteristic=[]
list_value=[]
list_sample=[]
list_class=[]
for line in open(input_path+sample_info,'r'):
    info = line[:-1].split('\t')
    id = info[0]
    if id == "sample_id":
        continue
    if id != 'ID':
        if info[1] == 'Class1':
            list_pos.append(id)
        read_file = open(input_path + id + ".txt", 'r')
        richness=0
        TotalCount=0
        for line in read_file:
            if line[0] != 'c':
                count = int(line[:-1].split('\t')[0])
                richness+=1
                TotalCount+=count
        read_file.close()
        read_file = open(input_path + id + ".txt", 'r')
        list_frequency = []
        for line in read_file:
            if line[0] != 'c':
                count = int(line[:-1].split('\t')[0])
                list_frequency.append(count/TotalCount)
        diversity=Shannon_entropy(list_frequency)
        evenness=Pielou_evenness(diversity,richness)
        clonality=Clonality(evenness)
        list_characteristic+=['diversity','evenness','clonality']
        list_value+=[diversity,evenness,clonality]
        for time in range(0,3):
            list_sample.append(id)
            if id in list_pos:
                list_class.append('Class1')
            else:
                list_class.append('Class2')
        read_file.close()

# print(list_class)
# ax1 = sns.barplot(x=list_characteristic,y=list_value,hue=list_sample,errwidth=0.5,capsize=0.3,palette=(sns.xkcd_rgb["dark red"],sns.xkcd_rgb["marine blue"]))
# ax1_1 = sns.stripplot(x=list_characteristic,y=list_value,size=2,hue=list_sample,dodge=True,palette=(sns.xkcd_rgb["dark red"],sns.xkcd_rgb["marine blue"]))
# ax1.spines['top'].set_visible(False)
# ax1.spines['right'].set_visible(False)
# plt.savefig(output_path+'diversity_evenness_clonality_barplot_class.pdf')
# plt.show()

# ax2 = sns.violinplot(x=list_characteristic,y=list_value,hue=list_class,palette=('white','white'))
# ax2_1 = sns.stripplot(x=list_characteristic,y=list_value,size=4,hue=list_class,dodge=True,palette=(sns.xkcd_rgb["dark red"],sns.xkcd_rgb["marine blue"]))
# ax2.spines['top'].set_visible(False)
# ax2.spines['right'].set_visible(False)
# plt.savefig(output_path+'diversity_evenness_clonality_violinplot_class.pdf')
# plt.show()
list_sample = [i.split('.')[0] + "_clonotypes" for i in list_sample]
ax3 = sns.barplot(x=list_characteristic,y=list_value,hue=list_sample,errwidth=0.5,capsize=0.5)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
plt.title("CDR3 diversity / evenness / clonality", fontsize=15)
plt.legend(bbox_to_anchor=(1.05, 1),loc=2,borderaxespad=0)
plt.savefig(output_path+'diversity_evenness_clonality_sample.pdf', dpi=300, bbox_inches="tight")
plt.show()


# CDR3 Clonotypes Abundance Proportion
input_path='/SGRNJ06/randd/USER/cjj/celedev/vdj_bulk/20230131demo/plot/vdjtools_result/'#输入路径
output_path='/SGRNJ06/randd/USER/cjj/celedev/vdj_bulk/20230131demo/plot/vdjtools_result/'#输出路径
sample_info='basicstats.txt'#记录有样本分组的样本信息
dict_sample_ReadsCounts={}
list_pos=[]
for line in open(input_path+sample_info,'r'):
    info = line[:-1].split('\t')
    id = info[0]
    if id == "sample_id":
        continue
    if id != 'ID':
        if str(info[1]) == 'Class1':
            list_pos.append(id)
        read_file = open(input_path + id + ".txt", 'r')
        dict_read_count={'1':0,'2-3':0,'4-10':0,'11-30':0,'31-100':0,'101-MAX':0}
        number_clonotypes=0
        for line in read_file:
            if line[0] != 'c':
                number_clonotypes+=1
                info=line[:-1].split('\t')
                count=info[0]
                if int(count) == 1:
                    dict_read_count['1']+=1
                elif int(count) in [2,3]:
                    dict_read_count['2-3']+=1
                elif 4 <= int(count) <= 10:
                    dict_read_count['4-10']+=1
                elif 11 <= int(count) <= 30:
                    dict_read_count['11-30']+=1
                elif 31 <= int(count) <= 100:
                    dict_read_count['31-100']+=1
                elif int(count) >= 101:
                    dict_read_count['101-MAX']+=1
        for key in dict_read_count:
            dict_read_count[key]=dict_read_count[key]/number_clonotypes*100
        dict_sample_ReadsCounts[id]=dict_read_count
        read_file.close()

list_sample=[]
list_class=[]
list_read=[]
list_count=[]
dict_read_SampleCount={'1': {},'2-3': {},'4-10': {},'11-30': {},'31-100': {},'101-MAX': {}}
print(dict_sample_ReadsCounts)
for sample in dict_sample_ReadsCounts.keys():
    dict_read_count=dict_sample_ReadsCounts[sample]
    for read in dict_read_count.keys():
        list_sample.append(sample)
        if sample in list_pos:
            list_class.append('Class1')
        else:
            list_class.append('Class2')
        list_read.append(read)
        count=float(dict_read_count[read])
        list_count.append(count)
        dict_sample_count=dict_read_SampleCount[read]
        dict_sample_count[sample]=count
        dict_read_SampleCount[read]=dict_sample_count

# ax1 = sns.barplot(x=list_read,y=list_count,hue=list_class,errwidth=1,capsize=0.2,palette=(sns.xkcd_rgb["dark red"],sns.xkcd_rgb["marine blue"]))
# ax1_1 = sns.stripplot(x=list_read,y=list_count,hue=list_class,dodge=True,palette=(sns.xkcd_rgb["dark red"],sns.xkcd_rgb["marine blue"]))
# ax1.set_xticklabels(ax1.get_xticklabels(), rotation=-90,fontsize=6)
# ax1.spines['top'].set_visible(False)
# ax1.spines['right'].set_visible(False)
# plt.savefig(output_path+'CountProportion_lappedplot_Class.pdf',dpi=300, bbox_inches="tight")#
# plt.show()

fig, ax = plt.subplots()
list_1=[]
list_2_3=[]
list_4_10=[]
list_11_30=[]
list_31_100=[]
list_101_MAX=[]
for read in dict_read_SampleCount.keys():
    dict_sample_count=dict_read_SampleCount[read]
    list_sample=[]
    for sample in dict_sample_count.keys():
        list_sample.append(sample)
        if read=='1':
            list_1.append(dict_sample_count[sample])
        elif read=='2-3':
            list_2_3.append(dict_sample_count[sample])
        elif read=='4-10':
            list_4_10.append(dict_sample_count[sample])
        elif read=='11-30':
            list_11_30.append(dict_sample_count[sample])
        elif read=='31-100':
            list_31_100.append(dict_sample_count[sample])
        elif read=='101-MAX':
            list_101_MAX.append(dict_sample_count[sample])
list_1=np.array(list_1)
list_2_3=np.array(list_2_3)
list_4_10=np.array(list_4_10)
list_11_30=np.array(list_11_30)
list_31_100=np.array(list_31_100)
list_101_MAX=np.array(list_101_MAX)
width = 0.6
print(list_1)
ax.bar(list_sample,list_1,width, label='1')
ax.bar(list_sample,list_2_3,width,label='2_3',bottom=list_1)
ax.bar(list_sample,list_4_10,width,label='4-10',bottom=list_1+list_2_3)
ax.bar(list_sample,list_11_30,width,label='11-30',bottom=list_1+list_2_3+list_4_10)
ax.bar(list_sample,list_31_100,width,label='31-100',bottom=list_1+list_2_3+list_4_10+list_11_30)
ax.bar(list_sample,list_101_MAX,width,label='101-MAX',bottom=list_1+list_2_3+list_4_10+list_11_30+list_31_100)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend()
list_sample = [i.split('.')[0] + "_clonotypes" for i in list_sample]
ax.set_xticklabels(list_sample, rotation=45, fontsize=10)
plt.legend(bbox_to_anchor=(1.05, 1),loc=2,borderaxespad=0)
plt.title("CDR3 Clonotypes Abundance Proportion", fontsize=15)
plt.savefig(output_path+'CountProportion_barplot_sample.pdf',dpi=300, bbox_inches="tight")
plt.show()


# count gene number
path='/SGRNJ06/randd/USER/cjj/celedev/vdj_bulk/20230131demo/plot/vdjtools_result/'

sample_info=open(path+'tmp/sample_info.txt','w')
sample_info.write('ID\tlabel\tbatch\n')
list_files=[]
for line in open(path+'/metadata.txt','r'):
    print(line)
    if line.split('\t')[0] == 'file_name':
        continue
    else:
        file=line.split('\t')[1]
        sample=file.split('.')[0]
        list_files.append(file+'.txt')
        sample_info.write(sample+'\t'+'Class1'+'\t'+'1'+'\n')
sample_info.close()
list_sample=[]
list_V=[]
list_J=[]
list_dict_V=[]
list_dict_J=[]
list_files = ["sampleE_clonotypes.txt", "sampleF_clonotypes.txt", "sampleG_clonotypes.txt", "sampleH_clonotypes.txt"]
for file in list_files:
    read_file=open(path+file,'r')
    sample=file.split('.')[0]
    list_sample.append(sample)
    dict_CDR={}
    dict_V={}
    dict_J={}
    for line in read_file:
        if line[0]!='c':
            info=line[:-1].split('\t')
            name_V=info[4]
            name_J=info[6]
            if name_V not in list_V:
                list_V.append(name_V)
            if name_J not in list_J:
                list_J.append(name_J)
            if name_V not in dict_V.keys():
                dict_V[name_V]=1
            else:
                dict_V[name_V]+=1
            if name_J not in dict_J.keys():
                dict_J[name_J]=1
            else:
                dict_J[name_J]+=1
    list_dict_V.append(dict_V)
    list_dict_J.append(dict_J)
    read_file.close()

word = 'item'
for sample in list_sample:
    word+='\t'+sample
word+='\n'
V_file=open('/SGRNJ06/randd/USER/cjj/celedev/vdj_bulk/20230131demo/plot/vdjtools_result/tmp/matrix_V.txt','w')
V_file.write(word)
for V in list_V:
    check = []
    info=str(V)
    for dict_V in list_dict_V:
        if V in dict_V.keys():
            info+='\t'+str(dict_V[V])
            check.append(float(dict_V[V]))
        else:
            info+='\t'+'0'
            check.append(0)
    info+='\n'
    if max(check)!=0:
        V_file.write(info)
V_file.close()
print('V gene expression matrix have been build!')
J_file=open('/SGRNJ06/randd/USER/cjj/celedev/vdj_bulk/20230131demo/plot/vdjtools_result/tmp/matrix_J.txt','w')
J_file.write(word)
for J in list_J:
    check = []
    info=str(J)
    for dict_J in list_dict_J:
        if J in dict_J.keys():
            info+='\t'+str(dict_J[J])
            check.append(float(dict_J[J]))
        else:
            info+='\t'+'0'
            check.append(0)
    info+='\n'
    if max(check)!=0:
        J_file.write(info)
J_file.close()
print('J gene expression matrix have been build!')
print('expression matrix have been build!')