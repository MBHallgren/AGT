import os
import sys

def fix_names(item):
    if '\'' in item:
        item = item.replace('\'', '')
    if '(' in item:
        item = '\'' + item + '\''
    if ')' in item:
        item = item.replace(')', '')
    return item

def realign_consensus(output, prefix, database, keep):
    non_perfect_hits = []
    alignment_dict = {}
    headers = ''

    headers, alignment_dict = load_kma_res_file('{}/{}.res'.format(output, prefix))


    for item in alignment_dict:
        if alignment_dict[item][3] != 100.00 or alignment_dict[item][4] != 100.00 or alignment_dict[item][5] != 100.00:
            item = fix_names(item)
            non_perfect_hits.append('>' + item)


    with open('{}/{}.fsa'.format(output, prefix), 'r') as f:
        flag = False
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                header = '>' + fix_names(line[1:])
                if line in non_perfect_hits:
                    flag = True
                else:
                    flag = False
            if flag:
                with open('{}/{}.fsa'.format(output, output_name), 'a') as f:
                    print (line, file=f)

    print (non_perfect_hits)
    sys.exit()

    for item in non_perfect_hits:
        cmd = 'kma -i {}/{}.fsa -o {}/{} -t_db {} -1t1 -proxi -0.95'.format(output, item[1:], output, item[1:], database)
        print (cmd)
        os.system(cmd)

    eval_realignments(output, prefix, headers, alignment_dict, non_perfect_hits)

    if not keep:
        for item in non_perfect_hits:
            if '\'' in item:
                item = item.replace('\'', '')
            if '(' in item:
                item = '\'' + item + '\''
            os.system('rm {}/{}*'.format(output, item[1:]))
        os.system('rm {}/old_*'.format(output, prefix))

def load_kma_res_file(file):
    kma_dict = dict()
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header = line.rstrip()
            else:
                line = line.rstrip()
                line = line.split('\t')
                kma_dict[line[0].strip()] = []
                for item in line[1:-1]:
                    kma_dict[line[0].strip()].append(float(item.rstrip()))
                kma_dict[line[0].strip()].append(line[-1].strip())
    return header, kma_dict
def eval_realignments(output, prefix, headers, alignment_dict, non_perfect_hits):
    realignment_dict = {}

    for item in alignment_dict:
        if alignment_dict[item][3] == 100.00 and alignment_dict[item][4] == 100.00 and alignment_dict[item][5] == 100.00: #Perfect alignment for Template_Identity	Template_Coverage	Query_Identity
            realignment_dict[item] = alignment_dict[item]

    for item in non_perfect_hits:
        headers, currect_gene_dict = load_kma_res_file('{}/{}.res'.format(output, item[1:]))
        original_gene = item[1:]
        for gene in currect_gene_dict:
            if gene not in realignment_dict:
                realignment_dict[gene] = alignment_dict[original_gene]
                realignment_dict[gene][3] = float(currect_gene_dict[gene][3]) #Replace template identity
                realignment_dict[gene][4] = float(currect_gene_dict[gene][4]) #Replace template coverage
                realignment_dict[gene][5] = float(currect_gene_dict[gene][5])  # Replace query identity
                realignment_dict[gene][6] = float(currect_gene_dict[gene][6])  # Replace Query_Coverage
            else:
                realignment_dict[gene][3] = float(currect_gene_dict[gene][3]) #Replace template identity
                realignment_dict[gene][4] = float(currect_gene_dict[gene][4]) #Replace template coverage
                realignment_dict[gene][5] = float(currect_gene_dict[gene][5])  # Replace query identity
                realignment_dict[gene][6] = float(currect_gene_dict[gene][6])  # Replace Query_Coverage
                realignment_dict[gene][7] = max(float(alignment_dict[original_gene][7]), float(currect_gene_dict[gene][7])) # Select max depth
                realignment_dict[gene][8] = max(float(alignment_dict[original_gene][8]), float(currect_gene_dict[gene][8])) # Select max q_value

    keys = list(realignment_dict.keys())
    keys.sort()

    with open('{}/final_{}.res'.format(output, prefix), 'w') as f:
        print (headers, file=f)
        for item in keys:
            print_string = '{}\t{}\t{}\t{}\t{}\t' \
                           '{}\t{}\t{}\t{}\t{}\t{}'\
            .format(item, int(realignment_dict[item][0]), int(realignment_dict[item][1]), int(realignment_dict[item][2]),
                    realignment_dict[item][3], realignment_dict[item][4], realignment_dict[item][5],
                    realignment_dict[item][6], realignment_dict[item][7], realignment_dict[item][8],
                    realignment_dict[item][9])

            print (print_string, file=f)

    os.system('mv {}/{}.res {}/old_{}.res'.format(output, prefix, output, prefix))
    os.system('mv {}/final_{}.res {}/{}.res'.format(output, prefix, output, prefix))


