import os
import pandas as pd
import argparse
import sys
import re
import warnings
#----------------------------------------------------------------------------------------#
def transform_args(value) :
    mapping = {
        'P' : 'Pathogenic.VUS',
        'A' : 'All.Variants'
    }
    if value in mapping :
        return mapping[value]
    else : 
        raise argparse.ArgumentTypeError(f"Invalid choice: {value} (choose from {list(mapping.keys())})")
#----------------------------------------------------------------------------------------#
def mkMAF(file, args, path) :
    global MAF
    for sheet in args.sheets : 
        tmp = pd.read_excel(f'{path}/{file}',sheet_name=sheet, engine='openpyxl')
        tmp = tmp[tmp['select'].isin(args.filter)]
        tmp.rename(columns={'gene':'Hugo_Symbol', 'ref':'Reference_Allele', 'alt':'Tumor_Seq_Allele2', 'VAF.var.freq':'i_TumorVAF_WU', 'NM':'i_transcript_name'}, inplace = True)
        tmp['Chromosome'] = tmp['chrom.pos'].apply(lambda x:x.split(':')[0].replace('chr',""))
        tmp['Start_Position'] = tmp['chrom.pos'].apply(lambda x:x.split(':')[1].split('-')[0])
        tmp['Start_Position'] = pd.to_numeric(tmp['Start_Position'])
        tmp['End_Position'] = tmp['chrom.pos'].apply(lambda x:x.split(':')[1].split('-')[1])
        tmp['End_Position'] = pd.to_numeric(tmp['End_Position'])
        tmp['Tumor_Sample_Barcode'] = file.split('.results.xlsx')[0]
        tmp['Variant_Classification'] = ""
        tmp['Variant_Type'] = ""
        for i in range(len(tmp['HGVSc'])) :
            if pd.isna(tmp['HGVSc'].iloc[i]):
                pass
            
            elif '>' in tmp['HGVSc'].iloc[i] :
                tmp.loc[tmp.index[i], 'Variant_Type'] = 'SNP'
                if pd.isna(tmp['HGVSp'].iloc[i]) :
                    if '+' in tmp['HGVSc'].iloc[i] or '-' in tmp['HGVSc'].iloc[i] :
                        intron = r'[+-](\d+)'
                        match = re.findall(intron, tmp['HGVSc'].iloc[i])
                        if match :
                            l = int(match[0])
                            if l == 1 or l == 2 :
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = 'Splice_Site'
                            else :
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = 'Intron_variant'
                        else :
                            continue
                        
                elif '=' in tmp['HGVSp'].iloc[i] :
                    tmp.loc[tmp.index[i], 'Variant_Classification'] = 'Silent_Mutation'
                elif 'Ter' in tmp['HGVSp'].iloc[i] :
                    tmp.loc[tmp.index[i], 'Variant_Classification'] = 'Nonsense_Mutation'
                else :
                    tmp.loc[tmp.index[i], 'Variant_Classification'] = 'Missense_Mutation'

            elif 'delins' in tmp['HGVSc'].iloc[i] :
                tmp.loc[tmp.index[i], 'Variant_Type'] = 'DELINS'
                if pd.isna(tmp['HGVSp'].iloc[i]) :
                    if '+' in tmp['HGVSc'].iloc[i] or '-' in tmp['HGVSc'].iloc[i] :
                        intron = r'[+-](\d+)'
                        match = re.findall(intron, tmp['HGVSc'].iloc[i])
                        if match :
                            l = int(match[0])
                            if l == 1 or l == 2 :
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = 'Splice_Site'
                            else :
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = 'Intron_variant'
                        else :
                            continue
                else :
                    if 'fs' in tmp['HGVSp'].iloc[i] :
                        tmp.loc[tmp.index[i], 'Variant_Classification'] = 'Frame_Shift_Indel'
                    else :
                        tmp.loc[tmp.index[i], 'Variant_Classification'] = 'In_Frame_Indel'

            elif 'del' in tmp['HGVSc'].iloc[i] :
                tmp.loc[tmp.index[i], 'Variant_Type'] = 'DEL'
                if pd.isna(tmp['HGVSp'].iloc[i]) :
                    if '+' in tmp['HGVSc'].iloc[i] or '-' in tmp['HGVSc'].iloc[i] :
                        intron = r'[+-](\d+)'
                        match = re.findall(intron, tmp['HGVSc'].iloc[i])
                        if match :
                            l = int(match[0])
                            if l == 1 or l == 2 :
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = 'Splice_Site'
                            else :
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = 'Intron_variant'
                        else :
                            continue
                else :
                    if 'fs' in tmp['HGVSp'].iloc[i] :
                        tmp.loc[tmp.index[i], 'Variant_Classification'] = 'Frame_Shift_Del'
                    else :
                        tmp.loc[tmp.index[i], 'Variant_Classification'] = 'In_Frame_Del'

            elif 'ins' in tmp['HGVSc'].iloc[i] :
                tmp.loc[tmp.index[i], 'Variant_Type'] = 'INS'
                if pd.isna(tmp['HGVSp'].iloc[i]) :
                    if '+' in tmp['HGVSc'].iloc[i] or '-' in tmp['HGVSc'].iloc[i] :
                        intron = r'[+-](\d+)'
                        match = re.findall(intron, tmp['HGVSc'].iloc[i])
                        if match :
                            l = int(match[0])
                            if l == 1 or l == 2 :
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = 'Splice_Site'
                            else :
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = 'Intron_variant'
                        else :
                            continue
                else :
                    if 'fs' in tmp['HGVSp'].iloc[i] :
                        tmp.loc[tmp.index[i], 'Variant_Classification'] = 'Frame_Shift_Ins'
                    else :
                        tmp.loc[tmp.index[i], 'Variant_Classification'] = 'In_Frame_Ins'
            
            elif 'dup' in tmp['HGVSc'].iloc[i] :
                # Chromosomal -> Structure_variant로 정의해야하는지 -> 그렇다면 deletion 길이 보고 특정 길이를 넘어서면 그것도 Structure_variant로 따로 분류해야하는 건 아닌지 + 교수님 파이프라인이 SNV랑 INDEL만 볼 수 있는거면 inversion뿐만 아니라 dup도 관찰할 수 없어야하는 거 아닌지..? dup 결과가 나왔다면 확인은 하지 못했지만 inv도 존재할 가능성..?
                tmp.loc[tmp.index[i], 'Variant_Type'] = 'INS' 
                tmp.loc[tmp.index[i], 'Variant_Classification'] = 'Duplication'
                
        tmp['Tumor_Sample_Barcode'] = file.split('.')[0]
        tmp.rename(columns={'HGVSp':'Protein_Change'}, inplace = True)
        tmp['i_TumorVAF_WU'] = pd.to_numeric(tmp['i_TumorVAF_WU'])
        tmp = tmp.assign(i_TumorVAF_WU=(tmp['i_TumorVAF_WU'] * 100).round(4))
        tmp = tmp[['Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Variant_Type','Tumor_Sample_Barcode','Protein_Change','i_TumorVAF_WU','i_transcript_name']]
        MAF = pd.concat([MAF, tmp], axis=0, ignore_index=True)
    # MAF.to_excel(f'{path}/OnGo.xlsx',index = False)
    MAF.to_csv(f'{path}/EasyOnco.maf', index = False, sep ='\t')
#----------------------------------------------------------------------------------------#
def run() : 
    warnings.simplefilter(action='ignore', category=FutureWarning)
    parser = argparse.ArgumentParser(description='OnGoPlotter Usage')
    parser.add_argument("-i", "--input", dest = "input", action = "store", nargs='+')
    parser.add_argument("-f", "--selected_filter", dest = "filter", action = "store", nargs='+', type = str, default=['v'])
    parser.add_argument("-s","--sheet", dest = "sheets", action = "store", nargs='+', choices=['P','A'], required = True, help='Pathogenic.VUS(P),All.Variants(A)')
    args = parser.parse_args()
    path = os.getcwd()
    args.sheets = list(map(transform_args, args.sheets))
    global MAF
    MAF = pd.DataFrame(columns=['Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Variant_Type','Tumor_Sample_Barcode','Protein_Change','i_TumorVAF_WU','i_transcript_name'])
    for File in args.input :
        mkMAF(File, args, path)

    command = 'Rscript ../Oncoplotter_v1.R'
    os.system(command)
#----------------------------------------------------------------------------------------#
if __name__ == '__main__' : 
    run()
