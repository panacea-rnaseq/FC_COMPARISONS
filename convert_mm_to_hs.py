import sys
import pandas as pd
from collections import defaultdict
from indra.databases.uniprot_client import um
from indra.databases.hgnc_client import get_hgnc_from_mouse, get_hgnc_name

path_to_mgi = sys.argv[1]
#path_to_mgi = "/Users/sbunga/PycharmProjects/INDRA/ligandReceptorInteractome/input/mgi_symbols.txt"
mgi = []
with open(path_to_mgi, 'r') as infile:
    for gene in infile:
        mgi.append(gene.strip())

def mgi_to_hgnc_name(genes):
    filtered_mgi = defaultdict(set)
    for gene in genes:
        if gene in mouse_gene_name_to_mgi:
            filtered_mgi[(gene)].add(mouse_gene_name_to_mgi[gene])
        else:
            filtered_mgi[(gene)].add('None')
    
    hgnc_gene_dict = dict()
    seen_genes = set()
    for key, value in filtered_mgi.items():
        mgi_id = next(iter(value))
        hgnc_id = get_hgnc_from_mouse(mgi_id)
        
        if hgnc_id == None or hgnc_id == 'None':
            hgnc_gene_dict[key] = 'None'
        else:
            hgnc_symbol = get_hgnc_name(hgnc_id)
            if hgnc_symbol not in seen_genes:
                hgnc_gene_dict[key] = hgnc_symbol
            seen_genes.add(hgnc_symbol)
    return hgnc_gene_dict

mouse_gene_name_to_mgi = {v: um.uniprot_mgi.get(k)
                          for k, v in um.uniprot_gene_name.items()
                          if k in um.uniprot_mgi}


n = mgi_to_hgnc_name(mgi)
mm_genes = list(n.keys())
hs_genes = list(n.values())

n = {
        'Mm_genes' : mm_genes,
        'Hs_genes': hs_genes
    }
n_df = pd.DataFrame(n)

n_df.to_csv("mm_hs_genes.csv", sep=",", header=True, index=False)
