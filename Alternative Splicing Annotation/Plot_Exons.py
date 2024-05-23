import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import ScalarFormatter
def annotate_cluster(cluster_file, exon_file, output_file):
    exon_info = {} 
    with open(exon_file, 'r', encoding='utf-8') as f:  
        for line in f:  
            chrom, start, end = line.strip().split('\t')  
            start, end = int(start), int(end) 
            if chrom not in exon_info.keys(): 
                exon_info[chrom] = [] 
            exon_info[chrom].append((start, end))    
    for key in exon_info:
        exon_info[key] = sorted(exon_info[key])


    cluster_positions = {}
    cluster_infosum = {}
    with open(cluster_file, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            cluster_info = parts[0]
            chrom, start, end, cluster = cluster_info.split(':')
            start, end = int(start), int(end)

            if (chrom, cluster) not in cluster_positions:
                cluster_positions[(chrom, cluster)] = []
            cluster_positions[(chrom, cluster)].append(start)
            cluster_positions[(chrom, cluster)].append(end)

            if (chrom, cluster) not in cluster_infosum:
                cluster_infosum[(chrom, cluster)] = []
            cluster_infosum[(chrom, cluster)].append((start, end))
    for key in cluster_positions:
        cluster_positions[key] = sorted(cluster_positions[key])

    for key in cluster_infosum:
        cluster_infosum[key] = sorted(cluster_infosum[key])


    with open(output_file, 'w', encoding='utf-8') as out:
        cluster_exons = {}
        for (chrom, cluster), positions in cluster_positions.items():
            min_position = min(positions) 
            max_position = max(positions)
            for chrom_exon, exons in exon_info.items():
                if chrom == chrom_exon:  
                    for i in range(len(exons) - 1):
                        min_exon_start, min_exon_end = exons[i]
                        min__1_exon_start, min__1_exon_end = exons[i + 1]  
                        if (min_exon_start <= min_position <= min_exon_end):
                            if (chrom_exon, cluster) not in cluster_exons.keys(): 
                                cluster_exons[chrom_exon, cluster] = []
                            cluster_exons[chrom_exon, cluster].append((min_exon_start, min_exon_end))
                            break                       
                        elif (min_exon_end < min_position < min__1_exon_start):
                            if (chrom_exon, cluster) not in cluster_exons.keys():
                                cluster_exons[chrom_exon, cluster] = []
                            cluster_exons[chrom_exon, cluster].append((min_exon_start, min_exon_end))
                            cluster_exons[chrom_exon, cluster].append((min__1_exon_start, min__1_exon_end))
                            break
                    for j in range(i, len(exons)):  
                        max_exon_start, max_exon_end = exons[j]
                        if j + 1 < len(exons):
                            max_1_exon_start, max_1_exon_end = exons[j + 1]
                        else:
                            max_1_exon_start, max_1_exon_end = max_exon_start, max_exon_end 
                        if (chrom_exon, cluster) not in cluster_exons.keys():
                            cluster_exons[chrom_exon, cluster] = []
                        cluster_exons[chrom_exon, cluster].append((max_exon_start, max_exon_end))
                        if (max_exon_start <= max_position <= max_exon_end):  
                            if (chrom_exon, cluster) not in cluster_exons.keys(): 
                                cluster_exons[chrom_exon, cluster] = []
                            cluster_exons[chrom_exon, cluster].append((max_exon_start, max_exon_end))
                            break
                        elif (max_1_exon_start <= max_position <= max_1_exon_end):  
                            if (chrom_exon, cluster) not in cluster_exons.keys():
                                cluster_exons[chrom_exon, cluster] = []
                            cluster_exons[chrom_exon, cluster].append((max_1_exon_start, max_1_exon_end))
                            break
                        
                        elif (max_exon_end < max_position < max_1_exon_start): 
                            if (chrom_exon, cluster) not in cluster_exons.keys():  
                                cluster_exons[chrom_exon, cluster] = []
                            cluster_exons[chrom_exon, cluster].append((max_1_exon_start, max_1_exon_end))
                            cluster_exons[chrom_exon, cluster].append((max_exon_start, max_exon_end))
                            break

        for key in cluster_infosum:
            values = sorted(set(cluster_infosum[key]))
            out.write(f"{key}:cluster_sum: {values}\n")
            if key in cluster_positions:
                values_positions = sorted(set(cluster_positions[key]))
                out.write(f"{key}:cluster_positions: {values_positions}\n")
                if key in cluster_exons:
                    values_exon = sorted(set(cluster_exons[key]))
                    out.write(f"{key}:cluster_exons: {values_exon}\n")
                    out.write(f"{key}:ASevent:\n\n")


    #plot
    pdf_file = "exon_clusters.pdf"
    with PdfPages(pdf_file) as pdf:
        for key, values_exon in cluster_exons.items():
            fig, ax = plt.subplots(figsize=(15, 5))
            for exon_start, exon_end in values_exon:
                rect = patches.Rectangle((exon_start, 0.9), exon_end - exon_start, 0.1, edgecolor='black',facecolor='blue')
                ax.add_patch(rect)
            min_position = min(set(cluster_positions[key]))
            max_position = max(set(cluster_positions[key]))
            ax.plot([min_position, min_position], [0.8, 1.2], color='red', linestyle='--', label='Min Position')
            ax.plot([max_position, max_position], [0.8, 1.2], color='green', linestyle='--', label='Max Position')
            y_base = 0.5
            y_step = 0.05
            for i, (start, end) in enumerate(cluster_infosum[key]):
                y_pos = y_base - i * y_step
                ax.plot([start, end], [y_pos, y_pos], color='orange')
                ax.plot([start, start], [y_pos - 0.5, y_pos + 0.5], color='purple', linestyle='--') 
                ax.plot([end, end], [y_pos - 0.5, y_pos + 0.5], color='purple', linestyle='--')
            ax.set_xlabel('Position')
            ax.set_ylabel('Exon')
            ax.set_title(f'{key}')
            ax.legend(loc='upper right')
            ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
            ax.ticklabel_format(style='plain', axis='x')
            plt.xticks(rotation=45)
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


exon_file = "HumanMergeExon.bed"
cluster_file = "AS.txt"
output_file = "annotated_position.txt"
annotate_cluster(cluster_file,exon_file,output_file)
