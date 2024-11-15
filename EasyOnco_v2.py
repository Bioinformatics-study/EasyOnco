import os
import pandas as pd
import argparse
import sys
import re
import warnings
from datetime import datetime
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
        col_names = tmp.columns.tolist()
        col_names[0] = 'select'
        tmp.columns = col_names
        pattern = '|'.join([re.sub(r'\*', r'.*', f) for f in args.filter]) 
        tmp['select'] = tmp['select'].astype(str)
        tmp = tmp[tmp['select']!='nan']
        tmp = tmp[tmp['select']!='NaN']
        tmp = tmp[tmp['select'].str.contains(pattern, regex=True, na=False)]
        tmp.rename(columns={'gene':'Hugo_Symbol', 'ref':'Reference_Allele', 'alt':'Tumor_Seq_Allele2', 'VAF.var.freq':'i_TumorVAF_WU', 'NM':'i_transcript_name'}, inplace = True)
        tmp['Chromosome'] = tmp['chrom.pos'].apply(lambda x: x.split(':')[0].replace('chr', "") if isinstance(x, str) else "")
        tmp['Start_Position'] = tmp['chrom.pos'].apply(lambda x: x.split(':')[1].split('-')[0] if isinstance(x, str) and ':' in x else "")
        tmp['Start_Position'] = pd.to_numeric(tmp['Start_Position'], errors='coerce')
        tmp['End_Position'] = tmp['chrom.pos'].apply(lambda x: x.split(':')[1].split('-')[1] if isinstance(x, str) and ':' in x else "")
        tmp['End_Position'] = pd.to_numeric(tmp['End_Position'], errors='coerce')
        
        tmp['Tumor_Sample_Barcode'] = file.split('.')[0]
        tmp['Variant_Classification'] = ""
        tmp['Variant_Type'] = ""
        for i in range(len(tmp['HGVSc'])) :
            if pd.isna(tmp['HGVSc'].iloc[i]):
                pass
            
            elif '>' in tmp['HGVSc'].iloc[i] :
                tmp.loc[tmp.index[i], 'Variant_Type'] = 'SNP'
                if pd.isna(tmp['HGVSp'].iloc[i]) :
                    if '+' in tmp['HGVSc'].iloc[i] or '-' or '*' in tmp['HGVSc'].iloc[i] :
                        UTR_pattern =  r'^(?P<dash>c\.-)|(?P<star>c\.\*)|(?P<one_dash>c\.1-)'
                        match = re.match(UTR_pattern, tmp['HGVSc'].iloc[i])
                        if match :
                            if match.group('dash'):
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = "5_prime_UTR_variant"
                            elif match.group('star') :
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = "3_prime_UTR_variant"
                            elif match.group('one_dash') :
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = "Promoter_variant"
                        else :
                            intron_pattern = r'.+(\d+)[+-](\d+).*'
                            match = re.findall(intron_pattern,tmp['HGVSc'].iloc[i])
                            if match :
                                l = int(match[0][1])
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
                    if '+' in tmp['HGVSc'].iloc[i] or '-' or '*' in tmp['HGVSc'].iloc[i] :
                        UTR_pattern =  r'^(?P<dash>c\.-)|(?P<star>c\.\*)|(?P<one_dash>c\.1-)'
                        match = re.match(UTR_pattern, tmp['HGVSc'].iloc[i])
                        if match :
                            if match.group('dash'):
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = "5_prime_UTR_variant"
                            elif match.group('star') :
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = "3_prime_UTR_variant"
                            elif match.group('one_dash') :
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = "Promoter_variant"
                        else :
                            intron_pattern = r'(\d+)[+-](\d+).*'
                            match = re.findall(intron_pattern,tmp['HGVSc'].iloc[i])
                            if match :
                                l = int(match[0][1])
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
                    if '+' in tmp['HGVSc'].iloc[i] or '-' or '*' in tmp['HGVSc'].iloc[i] :
                        UTR_pattern =  r'^(?P<dash>c\.-)|(?P<star>c\.\*)|(?P<one_dash>c\.1-)'
                        match = re.match(UTR_pattern, tmp['HGVSc'].iloc[i])
                        if match :
                            if match.group('dash'):
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = "5_prime_UTR_variant"
                            elif match.group('star') :
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = "3_prime_UTR_variant"
                            elif match.group('one_dash') :
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = "Promoter_variant"
                        else :
                            intron_pattern = r'(\d+)[+-](\d+).*'
                            match = re.findall(intron_pattern,tmp['HGVSc'].iloc[i])
                            if match :
                                l = int(match[0][1])
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
                    if '+' in tmp['HGVSc'].iloc[i] or '-' or '*' in tmp['HGVSc'].iloc[i] :
                        UTR_pattern =  r'^(?P<dash>c\.-)|(?P<star>c\.\*)|(?P<one_dash>c\.1-)'
                        match = re.match(UTR_pattern, tmp['HGVSc'].iloc[i])
                        if match :
                            if match.group('dash'):
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = "5_prime_UTR_variant"
                            elif match.group('star') :
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = "3_prime_UTR_variant"
                            elif match.group('one_dash') :
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = "Promoter_variant"
                        else :
                            intron_pattern = r'(\d+)[+-](\d+).*'
                            match = re.findall(intron_pattern,tmp['HGVSc'].iloc[i])
                            if match :
                                l = int(match[0][1])
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
                tmp.loc[tmp.index[i], 'Variant_Type'] = 'INS'
                if pd.isna(tmp['HGVSp'].iloc[i]) :
                    if '+' in tmp['HGVSc'].iloc[i] or '-' or '*' in tmp['HGVSc'].iloc[i] :
                        UTR_pattern =  r'^(?P<dash>c\.-)|(?P<star>c\.\*)|(?P<one_dash>c\.1-)'
                        match = re.match(UTR_pattern, tmp['HGVSc'].iloc[i])
                        if match :
                            if match.group('dash'):
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = "5_prime_UTR_variant"
                            elif match.group('star') :
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = "3_prime_UTR_variant"
                            elif match.group('one_dash') :
                                tmp.loc[tmp.index[i], 'Variant_Classification'] = "Promoter_variant"
                        else :
                            intron_pattern = r'(\d+)[+-](\d+).*'
                            match = re.findall(intron_pattern,tmp['HGVSc'].iloc[i])
                            if match :
                                l = int(match[0][1])
                                if l == 1 or l == 2 :
                                    tmp.loc[tmp.index[i], 'Variant_Classification'] = 'Splice_Site'
                                else :
                                    tmp.loc[tmp.index[i], 'Variant_Classification'] = 'Intron_variant'
                            else :
                                continue
                else :
                    tmp.loc[tmp.index[i], 'Variant_Classification'] = 'Duplication'
                
        # tmp['Tumor_Sample_Barcode'] = file.split('.')[0]
        tmp.rename(columns={'HGVSp':'Protein_Change'}, inplace = True)
        tmp['i_TumorVAF_WU'] = pd.to_numeric(tmp['i_TumorVAF_WU'])
        tmp = tmp.assign(i_TumorVAF_WU=(tmp['i_TumorVAF_WU'] * 100).round(4))
        tmp = tmp[['Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Variant_Type','Tumor_Sample_Barcode','Protein_Change','i_TumorVAF_WU','i_transcript_name']]
        MAF = pd.concat([MAF, tmp], axis=0, ignore_index=True)

#----------------------------------------------------------------------------------------#
def run() : 
    path = os.getcwd()
    args.sheets = list(map(transform_args, args.sheets))
    global MAF
    MAF = pd.DataFrame(columns=['Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Variant_Type','Tumor_Sample_Barcode','Protein_Change','i_TumorVAF_WU','i_transcript_name'])
    for File in args.input :
        if 'oncoplot_options' in File :
            pass
        else :
            with open(args.log_file, 'a') as log_file:
                now = datetime.now()
                log_file.write(f'Processing file: {File}\t({now.strftime("%Y-%m-%d %H:%M:%S")})\n')
                try:
                    mkMAF(File, args, path)
                except Exception as e:
                    log_file.write(f'Error file: {File}\t({now.strftime("%Y-%m-%d %H:%M:%S")})\n')
                    raise
    Sample_Names = [x.split('.')[0] for x in args.input]
    VariantO = list(MAF['Tumor_Sample_Barcode'].unique())
    VariantX = [x for x in Sample_Names if x not in VariantO]
    for variantx in VariantX :
        new_row = {'Tumor_Sample_Barcode': variantx}
        MAF = pd.concat([MAF, pd.DataFrame([new_row])], ignore_index=True)
    MAF.to_csv(f'{path}/{args.output}', index = False, sep ='\t')

    if args.mafonly == 'n' :
        EasyOnco_path = str(__file__).split('/')[:-1]
        EasyOnco_path = '/'.join(EasyOnco_path)
        command = f'Rscript {EasyOnco_path}/Oncoplotter_v3.R {args.output}'
        os.system(command)
    else :
        pass
#----------------------------------------------------------------------------------------#
if __name__ == '__main__' : 
    path = os.getcwd()
    warnings.simplefilter(action='ignore', category=FutureWarning)
    parser = argparse.ArgumentParser(description='OnGoPlotter Usage')
    parser.add_argument("-i", "--input", dest = "input", action = "store", nargs='+') #required=True
    parser.add_argument("-f", "--selected_filter", dest = "filter", action = "store", nargs='+',type = str) #required=True 
    parser.add_argument("-s","--sheet", dest = "sheets", action = "store", nargs='+', choices=['P','A'], help='Pathogenic.VUS(P),All.Variants(A)') #required=True
    parser.add_argument("-o","--output", dest = "output", action = "store", default = 'EasyOnco.maf')
    parser.add_argument("-m", "--maf-only", dest="mafonly", action="store", choices=['y', 'n'], default='n')
    parser.add_argument("-e", "--easyonco", dest="easyonco", action="store")
    parser.add_argument('--log_file', type=str, default=f'{path}/EasyOnco.log', help='Path to the log file.')
    args = parser.parse_args()

    if args.easyonco :
        EasyOnco_path = str(__file__).split('/')[:-1]
        EasyOnco_path = '/'.join(EasyOnco_path)
        command = f'Rscript {EasyOnco_path}/Oncoplotter_v3.R {args.easyonco}'
        os.system(command)
    else :
        run()