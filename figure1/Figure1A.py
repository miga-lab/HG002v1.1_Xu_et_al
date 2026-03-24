import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
import numpy as np
import argparse

# Set up argparse for output file
parser = argparse.ArgumentParser()

#this will compare the lengths of entries in two bed files 
parser.add_argument('--bed', '-b', type=str, action='store', help='censat annotations')
#parser.add_argument('--horhap', '-f', type=str, action='store', help='horhap file')
parser.add_argument('--cdr', '-c', type=str, action='store', help='cdr predictions')


args = parser.parse_args()
bed = args.bed
#horhap = args.horhap
cdr = args.cdr

# Load the BED file
censat = pd.read_csv(bed, sep="\t", header=None)
censat.columns = ["chrom", "start", "end", "name", "score", "strand", "start2", "end2", "color"]


# load the horhaps 
# horhap = pd.read_csv(horhap, sep="\t", header=None)
# horhap.columns = ["chrom", "start", "end", "name", "score", "strand", "start2", "end2", "color"]


# finally load cdr 
cdr = pd.read_csv(cdr, sep="\t", header=None)
cdr.columns = ["chrom", "start", "end", "name", "score", "strand", "start2", "end2", "color"]




active = censat[censat['name'].str.contains('active', case=False)]

normalizer = censat.groupby('chrom')['start'].min().to_dict() # find the minimum value as a start position


# Calculate a normalizing factor based on the start of the active array 

# Split the RGB color codes into lists of integers
censat['color'] = censat['color'].apply(lambda x: [int(c) / 255 for c in x.split(',')])
# horhap['color'] = horhap['color'].apply(lambda x: [int(c) / 255 for c in x.split(',')])
cdr['color'] = cdr['color'].apply(lambda x: [int(c) / 255 for c in x.split(',')])


# normalize the values 

normalizer = censat.groupby('chrom')['start'].min().to_dict() # find the minimum value as a start position

for df in [censat, cdr]:
    df[['start', 'end']] = df.apply(
    lambda row: [
        row['start'] - normalizer[row['chrom']],
        row['end'] - normalizer[row['chrom']],
    ], axis=1, result_type='expand')


# Determine the maximum normalized end for x-axis limits
max_normalized_end = censat['end'].max()
min_normalized_end = censat['start'].min()


# Set x-axis limits with 500,000 base pairs padding
x_min = min_normalized_end - 50000
x_max = max_normalized_end + 50000


# Plotting the annotations
fig, ax = plt.subplots(figsize=(20, 46))

# create a list of the chromosomes 
# set up the ordering for the plot 
chroms = ['Y', 'X'] + list(range(23, 0, -1)) 
chrom_labels = censat['chrom'].unique()  
sorted_chroms = []

for chr in chroms:
    chr = 'chr' + str(chr) + "_"
    sublist = [element for element in chrom_labels if chr in element]
    sorted_chroms.extend(sublist[::-1])



y_pos_dict = {chrom: i*2.5  for i, chrom in enumerate(sorted_chroms)}  # Create a dictionary for y positions

mat_orient = [-0.5, 0.5]
pat_orient = [0.75, -0.2]



for chrom in sorted_chroms:
    if "MATERNAL" in chrom:
        pos = mat_orient
    else:
        pos = pat_orient
    chrom_data = censat[censat['chrom'] == chrom]
    y_pos = y_pos_dict[chrom]+pos[0]  # Get the correct y position from the dictionary
    for _, row in chrom_data.iterrows():
        start = row['start']
        end = row['end']
        color = row['color']
        ax.barh(y_pos, end - start, left=start, height=1, color=color)
    
    # horhap_data = horhap[horhap['chrom'] == chrom]
    # y_pos = y_pos_dict[chrom]+pos[1]
    # for _, row in horhap_data.iterrows():
    #     start = row['start']
    #     end = row['end']
    #     color = row['color']
    #     ax.barh(y_pos, end - start, left=start, height=0.4, color=color)

    cdr_data = cdr[cdr['chrom'] == chrom]
    y_pos = y_pos_dict[chrom]+pos[1]
    for _, row in cdr_data.iterrows():
        start = row['start']
        end = row['end']
        color = row['color']
        ax.barh(y_pos, end - start, left=start, height=0.5, color=color)


# adding a small outline to the chrY because it is so hard to see 
y_start = (censat.loc[censat['chrom'] == 'chrY_PATERNAL', 'start'].min())
y_end = (censat.loc[censat['chrom'] == 'chrY_PATERNAL', 'end'].max())
y_length =  int(y_end - y_start)

ax.barh(y_pos_dict['chrY_PATERNAL']+pat_orient[0], width=y_length, left=y_start, height=1, color = 'none', edgecolor='black', linewidth=0.1)

# adding labels for chromosome 1 so folks know maternal and paternal 

# find the max end for chromosome 1 
max_chr1_pat = censat.loc[censat['chrom'] == 'chr1_PATERNAL', 'end'].max()
max_chr1_mat = censat.loc[censat['chrom'] == 'chr1_MATERNAL', 'end'].max()

ax.text(max_chr1_mat+100000,y_pos_dict['chr1_MATERNAL']+mat_orient[0]-0.25,  'Maternal', horizontalalignment = 'left', fontsize=18)
ax.text(max_chr1_pat+100000,y_pos_dict['chr1_PATERNAL']+pat_orient[0]-0.25,  'Paternal', horizontalalignment = 'left', fontsize=18)


# adding a ruler 
#calculate where to place it 
rul_xpos_left = x_max/2

# 5MB bar 
ax.barh((len(chrom_labels)*18)/10, 5000000, left=rul_xpos_left, height=1, color='black')
ax.text((rul_xpos_left+2500000),((len(chrom_labels)*18)/10)+0.9,  '5 MB', horizontalalignment = 'center', fontsize=18)

# tick marks for bar 
ax.barh((len(chrom_labels)*18)/10, 100000, left=rul_xpos_left, height=1.7, color='black')
ax.barh((len(chrom_labels)*18)/10, 100000, left=rul_xpos_left+1000000, height=1.5, color='black')
ax.barh((len(chrom_labels)*18)/10, 100000, left=rul_xpos_left+2000000, height=1.5, color='black')
ax.barh((len(chrom_labels)*18)/10, 100000, left=rul_xpos_left+3000000, height=1.5, color='black')
ax.barh((len(chrom_labels)*18)/10, 100000, left=rul_xpos_left+4000000, height=1.5, color='black')
ax.barh((len(chrom_labels)*18)/10, 100000, left=rul_xpos_left+5000000, height=1.7, color='black')

# Set the y-ticks with chromosome names and increase font size
ax.set_yticks([])

# Adjust the y-axis limits to reduce space at the top and bottom
ax.set_ylim(-1, len(chrom_labels)*2.5)


# Set x-axis limits based on the calculated min and max
ax.set_xlim(-x_min, x_max)

# Remove x-ticks and x-tick labels
ax.xaxis.set_ticks([])
ax.xaxis.set_ticklabels([])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

# Save the plot without the key
plt.savefig("Fig1A/HG002_censat_cdr.png", dpi=300, bbox_inches='tight')

# Creating the legend
color_table = {
    "Active HOR": "153,0,0",
    "Inactive HOR": "255,102,0",
    "Diverged HOR": "255,146,0",
    "monomeric alpha": "255,204,153",
    "mixedAlpha": "204,0,0",
    "Ribosomal DNA": "102,47,144",
    "Beta Satellite": "250,153,255",
    "Gamma Satellite": "172,51,199",
    "Other CenSat": "0,204,204",
    "HSat1A": "0,222,96",
    "HSat1B": "27,153,139",
    "HSat2": "0,128,250",
    "HSat3": "51,81,137",
    "Centromere Transition": "224,224,224",
    "GAP": "0,0,0",
    "CDR": "245,66,227"
}

legend_elements = [mpatches.Patch(facecolor=[int(x)/255 for x in color.split(',')], label=element)
                   for element, color in color_table.items()]

# Plot the legend
fig, ax = plt.subplots(figsize=(3, 4.2))
ax.legend(handles=legend_elements, loc='center', frameon=False)
ax.axis('off')
ax.set_title("CenSat Key", loc='center', fontsize=16)
plt.savefig("Fig1A/CenSat_legend_sorted.png", dpi=300, bbox_inches='tight')


# # Creating the HSAT subtypes 
# color_table = {
#     "2A1": "153,0,0",
#     "2A2": "204,51,51",
#     "2B": "255,102,102",
#     "3A1": "0,51,153",
#     "3A2": "51,102,204",
#     "3A3": "102,153,255",
#     "3A4": "153,204,255",
#     "3A5": "51,153,255",
#     "3A6": "0,102,204",
#     "3B1": "102,0,153",
#     "3B2": "153,51,204",
#     "3B3": "204,102,255",
#     "3B4": "255,153,255",
#     "3B5": "153,0,204",
#     "None": "128,128,128"

# }

# legend_elements = [mpatches.Patch(facecolor=[int(x)/255 for x in color.split(',')], label=element)
#                    for element, color in color_table.items()]

# # Plot the legend
# fig, ax = plt.subplots(figsize=(5, 5))
# ax.legend(handles=legend_elements, loc='center', frameon=False)
# ax.axis('off')
# plt.savefig("Fig1A/HSAT_legend_sorted.png", dpi=300, bbox_inches='tight')