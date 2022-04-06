
#!/usr/bin/env python3
# -*- coding: utf-8 -*-



description ="""
'this script is for establish sequence matrix from .gz meta_json of TCGA data STAR-counts, 
you may assign specific value on [biotype] to extract mRNA or lncRNA ,eg biotype = protein_coding, that would make mRNAmatrix'
#注意事项 1.压缩文件应该是.gz 结尾的文件 2.确保当前环境安装了 pandas和tqdm
@author: ApriPro zhuj ex hangzhou 
"""
biotype_types ="""
'IG_V_pseudogene', 'IG_J_pseudogene', 
'transcribed_unitary_pseudogene', 'scRNA', 'IG_J_gene', 
'TR_J_pseudogene', 'vault_RNA', 'TEC', 'snRNA', 'ribozyme', 
'scaRNA', 'translated_unprocessed_pseudogene', 'TR_V_pseudogene',
 'IG_V_gene', 'IG_pseudogene', 'Mt_tRNA', 'TR_D_gene', 'sRNA', 
 'IG_C_gene', 'misc_RNA', 'Mt_rRNA', 'snoRNA',
 'transcribed_processed_pseudogene', 
'translated_processed_pseudogene', 'unprocessed_pseudogene', 'IG_C_pseudogene', 
'pseudogene', 'IG_D_gene', 'rRNA', 'processed_pseudogene', 'TR_V_gene',
 'miRNA', 'polymorphic_pseudogene', 'transcribed_unprocessed_pseudogene', 
 'unitary_pseudogene', 'protein_coding', 'lncRNA', 'TR_J_gene', 'rRNA_pseudogene', 'TR_C_gene
 '"""



from tqdm import tqdm as tq
import os, time
import pandas as pd
import argparse as agp


parser = agp.ArgumentParser(description = description)
parser.add_argument('-j', '--json_path', type = str, help = 'meta.json file path, rudimentary.')
parser.add_argument('-t', '--tar_path', type = str, help = 'path of counts file, suffix .gz, rudimentary')
parser.add_argument('-b', '--biotype', type = str, help = f'sequence biotype, protein_coding(mRNA) as default, you may choose one from the following list(your input data could contains more types that are not indicated in the list):  {biotype_types}', default='protein_coding')

cmd_argus = parser.parse_args()


def sep(f):
    name = f.__code__.co_name
    start = time.time()
    print(f'start:|______{name}______')
    def wrapper(*argus,**kwargus):
        
        res = f(*argus, **kwargus)
        print(f'______{name}______|end|{time.time()-start}')
        
        return res
    return wrapper

@sep
def tar_(file = None, target_dir= 'unzipped'):
    if not os.path.exists(file):
        raise FileNotFoundError('could not find file when trying tar it')
    
    else:
        try:
            os.mkdir(f'{target_dir}')
            
        except FileExistsError:
            raise Warning(f'Directory exxists, over-writing the directory {target_dir}')
            
        finally:
            os.system(f'tar -zxf {file} -C {target_dir}')

    return os.path.abspath(target_dir)


@sep
def get_the_tsv(directory):
    
    data = []
    
    for root, directory, files in os.walk(directory):
        for file in files:
            if file.endswith('.tsv'):
                data.append(os.path.join(root,file))
    return data


def parser_tsv(tsv):
    
    df = pd.read_table(tsv, skiprows=5)
    tem = list(df.columns)
    tem[0:4]=['ens','gene','biotype','raw_data']
    
    df.columns = tem
    
    
    return df.iloc[:,1:4]


def choose_biotype_make_dict(df, biotype = 'protein_coding'):

    biotype = biotype
    df = df[df['biotype'] == biotype]
   
    assert len(df.index) != 0

    genes, raw_data = df.gene, df.raw_data
    return dict(zip(genes, raw_data))


@sep
def parser_json_for_sample_id(json_path = None):
    
    json_dict = eval(open(json_path).read())
    
    id_df = pd.json_normalize(data= json_dict,record_path='associated_entities')

    file_name_df = pd.json_normalize(data= json_dict,)
    file_name, entity_submitter_id = file_name_df.file_name, id_df.entity_submitter_id

    
    return dict(zip(file_name,entity_submitter_id))


@sep
def normal_sample_detector(df=None, columns=None, stat=True):
    from collections import Counter
    columns = df.columns

    columns = list(columns)

    test_v = list()
    for i in columns:
        try:

            if i.split('-')[3][0] == '1':

                test_v.append(False)
            else:
                test_v.append(True)

        except (Exception) as e:
            print(e, '\n', 'notation:', '\n',
                  '1.Encountering errors in columns-matching;', '\n',
                  f'2.For the sample named \'{i}\' in which TCGA sample_naming_method may be altered')
    if all(test_v):
        print('ALL NORMAL WARNINGS')
    if stat:
        print( Counter(test_v))
        #要存的话 加下面这句
        pd.DataFrame([Counter(test_v)]).to_csv('report.csv')
        
    return test_v
    print(Counter(test_v))
    
    
@sep
def rearrange_columns_of_matrix(matrix):
        temp_col = normal_sample_detector(df=matrix)
        matrix.loc['temp_v'] = temp_col
        new_col_order = matrix.T.sort_values(by='temp_v').index
        matrix = matrix.loc[:, new_col_order].drop('temp_v')
        print('matrix columns rearranged, tumor samples follows the normal ones.')
        return matrix


if __name__ =='__main__':

    biotype = cmd_argus.biotype

    json_path = cmd_argus.json_path

    tar_path = cmd_argus.tar_path
    

    dir = tar_(tar_path)
    print(f'unzppied in to {dir}')
   
    res = get_the_tsv(dir)
    def tsvpools():
        for i in tq(res):
            v = parser_tsv(i)
            d = choose_biotype_make_dict(v,biotype=biotype)
            yield d

    res2 = [os.path.basename(i) for i in res]
    
    data = dict(zip(res2,tsvpools()))
    meta_dict = parser_json_for_sample_id(json_path=json_path)

    data = pd.DataFrame(data)
    data.columns = data.columns.map(meta_dict)
    data=rearrange_columns_of_matrix(data)
    data.index.name='id'

    data.to_csv(f'{biotype}_mRNAmatrix.txt')
