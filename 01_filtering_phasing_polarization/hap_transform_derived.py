import os
import shutil

def hap_tranform_basedon_derived(hap_path, map_path, ref_derived_id_set):
    # Read hap file line by line
    with open(hap_path) as f:
        hap_lines = f.readlines()

    # Read map file line by line
    with open(map_path) as f:
        map_lines = f.readlines()

    # Create ids list: chrom:pos
    id_index = [line.strip().split()[0] + ":" + line.strip().split()[1] for line in map_lines]

    # Transform haplotypes
    transformed_hap = []
    for idx, line in enumerate(hap_lines):
        alleles = list(map(int, line.strip().split()))
        if id_index[idx] in ref_derived_id_set:
            alleles = [1 - a for a in alleles]  # flip allele
        transformed_hap.append(" ".join(map(str, alleles)))

    return transformed_hap

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        prog="hap_transform_derived.py",
        description="Recode phased haplotypes so that 0 = ancestral and 1 = derived"
    )
    parser.add_argument('-sp', "--species", action='store', help='species')
    parser.add_argument('-chr', "--chromosome", action='store', help='chromosome')
    args = parser.parse_args()

    # Load ancestral list into set for faster lookup
    with open('./alt_ancestral.txt') as f:
        ref_derived_id_set = set(line.strip() for line in f)

    sp = args.species
    chr = args.chromosome
    print(f"Processing chromosome: {chr}")

    hap_path = f"./hapmap/{sp}/{chr}.phase.hap"
    map_path = f"./hapmap/{sp}/{chr}.phase.map"

    transformed_hap = hap_tranform_basedon_derived(hap_path, map_path, ref_derived_id_set)

    # Save results
    out_dir = f"./new_hapmap/{sp}/"
    os.makedirs(out_dir, exist_ok=True)
    with open(f"{out_dir}/{chr}.phase.hap", "w") as f:
        f.write("\n".join(transformed_hap))

    # Copy other files
    for ext in ['map', 'samples', 'legend']:
        shutil.copy(f"./hapmap/{sp}/{chr}.phase.{ext}", f"{out_dir}/{chr}.phase.{ext}")
