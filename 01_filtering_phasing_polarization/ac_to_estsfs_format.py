import re
import random

orders = ['A','C','G','T']

def allele_sum(x):
    if x in ['./.', '', None]:
        return None
    return sum(int(i) for i in re.split(r'[|/]', x))


def ac_to_est_format(path_s, path_o, f_out, f_major):

    with open(path_s) as fs, open(path_o) as fo:

        next(fs)  # skip header
        next(fo)

        for s_line, o_line in zip(fs, fo):
            s = s_line.strip().split('\t')
            o = o_line.strip().split('\t')
            chrom, pos, ref, alt = s[:4]
            geno_sum = 0
            non_missing = 0
            for g in s[4:]:
                v = allele_sum(g)
                if v is not None:
                    geno_sum += v
                    non_missing += 1
            # major allele
            if geno_sum > non_missing:
                major = 'ALT'
            elif geno_sum < non_missing:
                major = 'REF'
            else:
                major = 'nomajor'
            f_major.write(
                f"{chrom}\t{pos}\t{ref}\t{alt}\t{geno_sum}\t{non_missing}\t{major}\n"
            )
            # focal
            s_out = [0,0,0,0]
            s_out[orders.index(alt)] = geno_sum
            s_out[orders.index(ref)] = non_missing*2 - geno_sum
            # outgroups
            o1_out=[0,0,0,0]
            o2_out=[0,0,0,0]
            for g, out in zip(o[4:6], [o1_out,o2_out]):
                v = allele_sum(g)
                if v == 1:
                    allele = random.choice([alt,ref])
                    out[orders.index(allele)] = 1
                elif v == 2:
                    out[orders.index(alt)] = 1
                else:
                    out[orders.index(ref)] = 1
            f_out.write(
                ",".join(map(str,s_out))+"\t"+
                ",".join(map(str,o1_out))+"\t"+
                ",".join(map(str,o2_out))+"\n"
            )
import argparse
parser = argparse.ArgumentParser(
    prog="ac_to_estsfs_format.py",
    description="transform allele count file to est-sfs input"
)

parser.add_argument(
    '-com','--outgroup_combination',
    help='prefix of files'
)

args = parser.parse_args()
f_out=open(f"{args.outgroup_combination}_data_file.txt",'w')
f_major=open(f"{args.outgroup_combination}_major.txt",'w')
f_major.write("CHROM\tPOS\tREF\tALT\tgeno_sum_s\tnonmissing\tmajor\n")
path_s=f"{args.outgroup_combination}.samples.count"
path_o=f"{args.outgroup_combination}.out.count"
ac_to_est_format(path_s,path_o,f_out,f_major)
f_out.close()
f_major.close()
