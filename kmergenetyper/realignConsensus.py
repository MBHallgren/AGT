import os
import sys

def realign_consensus(output, prefix, database, keep):
    non_perfect_hits = []
    alignment_dict = {}
    headers = ''
    with open('{}/{}.res'.format(output, prefix), 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('#'):
                line = line.split('\t')
                for item in line:
                    item = item.strip()
                if float(line[4]) < 100.00:
                    non_perfect_hits.append('>' + line[0].strip())
                alignment_dict[line[0]] = line[1:]
            else:
                headers = line

    alignment_dict = reformat_dict(alignment_dict)

    with open('{}/{}.fsa'.format(output, prefix), 'r') as f:
        flag = False
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                header = line
                if line in non_perfect_hits:
                    flag = True
                else:
                    flag = False
            if flag:
                with open('{}/{}.fsa'.format(output, header[1:]), 'a') as f:
                    print (line, file=f)

    for item in non_perfect_hits:
        os.system('kma -i {}/{}.fsa -o {}/{} -t_db {} -1t1 -proxi -0.95'.format(output, item[1:], output, item[1:], database))

    eval_realignments(output, prefix, headers, alignment_dict, non_perfect_hits)

    if not keep:
        for item in non_perfect_hits:
            os.system('rm {}/{}*'.format(output, item[1:]))
        os.system('rm {}/old_*'.format(output, prefix))

def reformat_dict(input_dict):
    output_dict = {}
    for item in input_dict:
        output_dict[item.strip()] = []
        for value in input_dict[item]:
            output_dict[item.strip()].append(value.strip())
    return output_dict
def eval_realignments(output, prefix, headers, alignment_dict, non_perfect_hits):
    realignment_dict = {}

    for item in alignment_dict:
        if float(alignment_dict[item][3]) == 100.00:
            realignment_dict[item] = alignment_dict[item]

    for item in non_perfect_hits:
        with open('{}/{}.res'.format(output, item[1:]), 'r') as f:
            original_gene = item[1:]
            for line in f:
                line = line.rstrip()
                if not line.startswith('#'):
                    gene = line.split('\t')[0]
                    if gene not in realignment_dict:
                        realignment_dict[gene] = alignment_dict[original_gene]
                        t_id_1 = realignment_dict[gene][3]
                        t_id_2 = line.split('\t')[4]
                        if float(t_id_1) < float(t_id_2):
                            realignment_dict[gene][3] = t_id_2 # Replace template identity
                            realignment_dict[gene][4] = line.split('\t')[5]  # Replace template coverage
                            realignment_dict[gene][5] = line.split('\t')[6]  # Replace template coverage
                    else:
                        t_id_1 = realignment_dict[gene][3]
                        t_id_2 = line.split('\t')[4]
                        if float(t_id_1) < float(t_id_2):
                            realignment_dict[gene][3] = t_id_2 #Replace template identity
                            realignment_dict[gene][4] = line.split('\t')[5] #Replace template coverage
                            realignment_dict[gene][5] = line.split('\t')[6]  # Replace template coverage
    realignment_dict = reformat_dict(realignment_dict)

    keys = list(realignment_dict.keys())
    keys.sort()

    with open('{}/final_{}.res'.format(output, prefix), 'w') as f:
        print (headers, file=f)
        for item in keys:
            print_list = [item] + realignment_dict[item]
            print_list = "\t".join(print_list)
            print (print_list, file=f)

    os.system('mv {}/{}.res {}/old_{}.res'.format(output, prefix, output, prefix))
    os.system('mv {}/final_{}.res {}/{}.res'.format(output, prefix, output, prefix))


