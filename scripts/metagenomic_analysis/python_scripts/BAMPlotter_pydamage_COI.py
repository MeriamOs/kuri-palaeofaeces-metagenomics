#!/usr/bin/env python3
import os
import argparse

parser = argparse.ArgumentParser(
    description="""aDNA BAM Mapping Plots For Microbial Genomes -  Meriam Guellil  -  September 2021 v2.0""",
    epilog="""Outputs coverage, edit distance and deamination plots for each BAM header"""
)
parser.add_argument('-b', metavar='BAM file', dest='bamR', required=True, type=str,
                    help='Indexed BAM file, which should ideally not be filtered based on MAPQ (required)')
parser.add_argument('-d', metavar='Misincorporation file', dest='deam', required=False, type=str,
                    help='mapDamage2 misincorporation.txt file or new CSV (optional)')
parser.add_argument('-o', metavar='Output File', dest='out', required=True, type=str,
                    help='Output file with extension for desired format (e.g. pdf, svg, png) (required)')
parser.add_argument('-i', metavar='Headers for output', dest='headlist', required=False, type=str, nargs='+',
                    help='Space seperated list of BAM headers to filter for (optional)')
parser.add_argument('-q', metavar='desired MQ threshold', dest='mqf', required=False, default="30", type=str,
                    help='MQ threshold for quality filtered coverage plot (default: 30)')
args = parser.parse_args()

import datetime
import numpy as np
import pysam
import pysamstats
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from io import StringIO
from tqdm import tqdm

date = datetime.datetime.now().strftime("%d/%m/%Y")

# Load BAM
bam = pysam.AlignmentFile(args.bamR, "rb")
head = {n["SN"]: n["LN"] for n in bam.header['SQ']}

# LOAD MISINCORPORATION
deam_df = None
deam_format = None
if args.deam:
    deam_df = pd.read_csv(args.deam, header=0, delimiter=',', encoding='utf-8')
    if 'Chr' in deam_df.columns and 'End' in deam_df.columns:
        deam_format = 'mapDamage2'
    else:
        deam_format = 'newcsv'

# filter if available
if args.headlist:
    head = {k: v for k, v in head.items() if k in args.headlist}

# Filter head for >10 mapped reads
filtered_head = {}
for ID, LEN in head.items():
    try:
        read_count = int(pysam.view("-c", "-F", "4", args.bamR, "-u", ID))
    except Exception:
        read_count = 0
    if read_count > 10:
        filtered_head[ID] = LEN

# LOAD DEPTH
depth = pd.read_csv(StringIO(pysam.depth(args.bamR, "-aa")), delimiter='\t', encoding='utf-8',
                    names=["id", "pos", "dp"])
depth30 = pd.read_csv(StringIO(pysam.depth(args.bamR, "-aa", "-Q", args.mqf)), delimiter='\t', encoding='utf-8',
                      names=["id", "pos", "dp"])

# Calculate maximum depth per sample
max_depth_per_sample = depth.groupby('id')['dp'].max().to_dict()

plt.rcParams.update({'font.size': 11, 'figure.titlesize': 16, 'figure.titleweight': "bold", 'legend.fontsize': 11})
fig = plt.figure(constrained_layout=False)
gs = gridspec.GridSpec(nrows=max(1, len(filtered_head)), ncols=4, figure=fig)
fig.set_figheight(((3.5 * max(1, len(filtered_head)))))
fig.set_figwidth(25)
fig.suptitle('BAM File: ' + os.path.basename(args.bamR))

count = 0

pbar = tqdm(filtered_head.items(), desc="Progress", unit=" Plot", total=len(filtered_head))
for ID, LEN in pbar:
    count += 1
    deamn3 = {}
    deamn5 = {}

    # get mapped read counts (ensure int)
    try:
        read_count = int(pysam.view("-c", "-F", "4", args.bamR, "-u", ID))
    except Exception:
        read_count = 0
    try:
        read_count_mqf = int(pysam.view("-c", "-F", "4", "-q", args.mqf, args.bamR, "-u", ID))
    except Exception:
        read_count_mqf = 0

    if int(read_count) > 0:
        # Coverage Stats
        depthF = depth[depth['id'].values == ID]
        per_covered = (((depthF['dp'].values != 0).astype(int).sum(axis=0)) / LEN) * 100
        perc2x = (((depthF['dp'].values >= 2).astype(int).sum(axis=0)) / LEN) * 100
        meandepth = depthF["dp"].mean() if not depthF.empty else 0.0
        ## Coverage Stats MQ threshold
        depthMQ30 = depth30[depth30['id'].values == ID]

        # DoC Stats with safe binning (preserve mean pos per bin)
        bin_size = max(1, int(LEN * 0.002))
        Cov_df = depthF[['pos', 'dp']].copy()
        Cov_df['bin'] = Cov_df['pos'] // bin_size
        Cov_Simpl = Cov_df.groupby('bin', as_index=False).mean()  # columns: bin, pos, dp

        Cov_df_MQ30 = depthMQ30[['pos', 'dp']].copy()
        Cov_df_MQ30['bin'] = Cov_df_MQ30['pos'] // bin_size
        Cov_Simpl_MQ30 = Cov_df_MQ30.groupby('bin', as_index=False).mean()

        # MQ Stats
        bamh = bam.fetch(ID)
        MQ_list = [read.mapping_quality for read in bamh]
        mean = float(np.average(MQ_list)) if len(MQ_list) > 0 else 0.0

        # MQ Zero List (positions where rms_mapq == 0 and reads_all != 0)
        mq = pysamstats.load_mapq_binned(bam, chrom=ID, window_size=bin_size)
        MQ = {i.pos: 0 for i in mq if (getattr(i, "rms_mapq", None) == 0 and getattr(i, "reads_all", 0) != 0)}

        # ED Stats MQ0 (use get_tag safely)
        bamh = bam.fetch(ID)
        NM_list = []
        for read in bamh:
            try:
                NM_list.append(read.get_tag("NM"))
            except KeyError:
                continue
        if len(NM_list) > 0:
            ED_Dict = {x: NM_list.count(x) / len(NM_list) for x in range(0, 6)}
            mean_ED = float(np.average(NM_list))
        else:
            ED_Dict = {x: 0 for x in range(0, 6)}
            mean_ED = 0.0

        ## ED Stats MQ filtered
        bamh = bam.fetch(ID)
        NM_list_MQ30 = []
        for read in bamh:
            try:
                if read.mapping_quality >= int(args.mqf):
                    NM_list_MQ30.append(read.get_tag("NM"))
            except KeyError:
                continue
        if len(NM_list_MQ30) > 0:
            ED_Dict_MQ30 = {x: NM_list_MQ30.count(x) / len(NM_list_MQ30) for x in range(0, 6)}
            mean_ED_MQ30 = float(np.average(NM_list_MQ30))
        else:
            ED_Dict_MQ30 = {x: 0 for x in range(0, 6)}
            mean_ED_MQ30 = 0.0

        # Deamin Stats
        if args.deam and deam_format is not None:
            if deam_format == 'mapDamage2':
                for i in range(1, 6):
                    df_5p = deam_df[(deam_df['Chr'].values == ID) & (deam_df['End'].values == "5p") & (deam_df['Pos'].values == i)]
                    if not df_5p.empty and df_5p['C'].sum() > 0:
                        mean5 = df_5p['C>T'].sum() / df_5p['C'].sum()
                    else:
                        mean5 = 0
                    deamn5.update({i: mean5})
                    df_3p = deam_df[(deam_df['Chr'].values == ID) & (deam_df['End'].values == "3p") & (deam_df['Pos'].values == i)]
                    if not df_3p.empty and df_3p['G'].sum() > 0:
                        mean3 = df_3p['G>A'].sum() / df_3p['G'].sum()
                    else:
                        mean3 = 0
                    deamn3.update({i: mean3})
            else:  # new CSV format
                deamn5 = {}
                deamn3 = {}
                ref_row = deam_df.loc[deam_df['reference'] == ID]
                if not ref_row.empty:
                    for i in range(1, 6):
                        col_name = f'CtoT-{i - 1}'
                        if col_name in ref_row.columns:
                            val = float(ref_row.iloc[0][col_name])
                            deamn5[i] = val

        # Plotting
        f = fig.add_subplot(gs[(count - 1), :3])

        # If binned data exists, plot mean position (pos) vs mean dp
        if not Cov_Simpl_MQ30.empty:
            plt.plot(Cov_Simpl_MQ30["pos"], Cov_Simpl_MQ30['dp'], color='#148F77', alpha=1)
            if len(Cov_Simpl_MQ30) > 1:
                f.fill_between(Cov_Simpl_MQ30["pos"], Cov_Simpl_MQ30['dp'], color='#148F77', alpha=0.4)

        if not Cov_Simpl.empty:
            plt.plot(Cov_Simpl["pos"], Cov_Simpl['dp'], color='gray', alpha=0.6)
            if len(Cov_Simpl) > 1:
                f.fill_between(Cov_Simpl["pos"], Cov_Simpl['dp'], color='gray', alpha=0.2)

        f.axes = plt.gca()
        max_dp = max_depth_per_sample.get(ID, max(1, int(np.nanmax(depthF['dp'].values)) if not depthF.empty else 10))
        f.axes.set_ylim([0, max_dp])
        f.axes.set_xlim([0, LEN])
        ax2 = f.twinx()
        ax2.set_ylim([-1, 1])
        ax2.get_yaxis().set_ticks([])

        # MQ zeros (scatter) - keys might be 0-based positions depending on pysamstats
        if MQ:
            plt.scatter(list(MQ.keys()), list(MQ.values()), alpha=0.2)

        # Text annotations (fixed unmatched-paren issue)
        f.text((LEN / 100) * 0.6, max_dp * 0.95, ID, ha='left', va='top', weight="bold")
        f.text((LEN / 100) * 0.6, max_dp * 0.9,
               'Mean Depth (MQ>=0): ' + "{0:.2f}".format(meandepth),
               ha='left', va='top')
        f.text((LEN / 100) * 0.6, max_dp * 0.85,
               'Mapped Reads (MQ>=0): ' + f"{int(read_count):,}",
               ha='left', va='top')
        f.text((LEN / 100) * 0.6, max_dp * 0.8,
               'Mapped Reads (MQ>=' + args.mqf + '): ' + f"{int(read_count_mqf):,}",
               ha='left', va='top')
        f.text((LEN / 100) * 0.6, max_dp * 0.75,
               'Cov% (MQ>=0): ' + "{0:.2f}".format(per_covered),
               ha='left', va='top')
        f.text((LEN / 100) * 0.6, max_dp * 0.7,
               '%2X (MQ>=0): ' + "{0:.2f}".format(perc2x),
               ha='left', va='top')
        f.text((LEN / 100) * 0.6, max_dp * 0.65,
               'Mean MAPQ (>=0): ' + "{0:.2f}".format(mean),
               ha='left', va='top')

        f.axes.set_ylabel('')
        f.axes.set_xlabel('')
        legend_elements1 = [
            Line2D([0], [0], marker='o', color='#148F77', label='CovMQ>=' + args.mqf,
                   markerfacecolor='#148F77', markersize=6),
            Line2D([0], [0], marker='o', color='gray', label='CovMQ>=0',
                   markerfacecolor='gray', markersize=6),
            Line2D([0], [0], marker='o', color='#DB414D', label='MQ0',
                   markerfacecolor='#DB414D', markersize=6, alpha=0.3)
        ]
        f.legend(handles=legend_elements1, bbox_to_anchor=(1, 1), loc="upper right", fancybox=False, framealpha=0.4)

        e = fig.add_subplot(gs[(count - 1), 3])
        plt.bar(list(ED_Dict_MQ30.keys()), list(ED_Dict_MQ30.values()), color="#148F77", alpha=0.5, width=0.8)
        plt.bar(list(ED_Dict.keys()), list(ED_Dict.values()), alpha=0.3, color="gray", width=0.5)

        if args.deam and deamn5:
            ax3 = e.twinx()
            plt.plot(list(deamn5.keys()), list(deamn5.values()), color="#E15261", alpha=0.7)
            # If deamn3 has values, plot them too
            if deamn3:
                plt.plot(list(deamn3.keys()), list(deamn3.values()), color="#535F9F", alpha=0.7)
            e.text(5.6, max(list(ED_Dict.values()) or [0]), 'Mean ED (MQ>=0): ' + "{0:.3f}".format(mean_ED),
                   ha='right', va='top')
            legend_elements2 = [
                Line2D([0], [0], marker='o', color="#148F77", label='EDistMQ>' + args.mqf,
                       markerfacecolor='#148F77', markersize=6),
                Line2D([0], [0], marker='o', color="gray", label='EDistMQ<' + args.mqf,
                       markerfacecolor='gray', markersize=6, alpha=0.3),
                Line2D([0], [0], marker='o', color="#E15261", label='5pCtoT',
                       markerfacecolor='#E15261', markersize=6),
                Line2D([0], [0], marker='o', color="#535F9F", label='3pGtoA',
                       markerfacecolor='#535F9F', markersize=6)
            ]
            e.legend(handles=legend_elements2, frameon=False, bbox_to_anchor=(1, 0.9), loc="upper right")
        else:
            e.text(5.6, max(list(ED_Dict.values()) or [0]), 'Mean ED (MQ>=0): ' + "{0:.3f}".format(mean_ED),
                   ha='right', va='top')
            legend_elements2 = [
                Line2D([0], [0], marker='o', color="#148F77", label='EDistMQ>' + args.mqf,
                       markerfacecolor='#148F77', markersize=6),
                Line2D([0], [0], marker='o', color="gray", label='EDistMQ<' + args.mqf,
                       markerfacecolor='gray', markersize=6, alpha=0.3)
            ]
            e.legend(handles=legend_elements2, frameon=False, bbox_to_anchor=(1, 0.9), loc="upper right")

    else:
        # No mapping case
        f = fig.add_subplot(gs[(count - 1), :3])
        f.set_ylim([0, 10])
        f.set_xlim([0, LEN])
        f.set_ylabel('')
        f.set_xlabel('')
        f.text((LEN / 100) * 1.5, 9.5, ID, ha='left', va='top', weight="bold")
        f.text((LEN / 100) * 20, 9.5, "NO MAPPING!", ha='left', va='top', weight="bold")
        e = fig.add_subplot(gs[(count - 1), 3])

# Footer and save
fig.text(0.99, 0.99, date, ha='right', va='top', fontsize=9)
fig.tight_layout()
fig.savefig(args.out, bbox_inches='tight')
