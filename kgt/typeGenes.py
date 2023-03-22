import os
import sys

def type_genes(args):
    set_up_output_directory(args.output)
    if args.illumina != []:
        args.illumina.sort()
        for i in range(0, len(args.illumina), 2):
            prefix = derive_prefix(args.illumina[i])
            input_seqs = ' '.join(args.illumina[i:i+2])
            os.system('kma -ipe {} -o {} -t_db {} -ill -md {} -ID 95'.format(input_seqs, args.output, args.t_db, args.md))
            #realign_consensus(args.output, prefix, args.t_db, input_seqs)
    if args.nanopore != []:
        for item in args.nanopore:
            prefix = derive_prefix(item)
            os.system('kma -i {} -o {} -t_db {} -ont -md {} -ID 95'.format(item, args.output, args.t_db, args.md))
            #realign_consensus(args.output, prefix, args.t_db, item)

def derive_prefix(file):
    return os.path.basename(file).split('.')[0]

def set_up_output_directory(output):
    if not os.path.exists(output):
        os.makedirs(output)