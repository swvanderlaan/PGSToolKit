#!/hpc/local/Rocky8/dhl_ec/software/mambaforge3/envs/pgstoolkit/bin/python

### DEPRECATED
### #!/hpc/local/CentOS7/dhl_ec/software/python37/bin/python3.7

import gzip
import argparse

# Construct argument parser
ap = argparse.ArgumentParser()

# Add the arguments to the parser
ap.add_argument("-b", "--base", required=True, help="Path to the base file containing the GWAS summary statistics")
ap.add_argument("-i", "--bf_idcol", required=True, help="Name of the column containing the variant IDs in the base file")
ap.add_argument("-s", "--stats", required=True, help="Path to the stats file containing the MAF and info scores")
ap.add_argument("-a", "--stats_idcol", required=True, help="Name of the column containing the variant IDs in the stats file")
ap.add_argument("-c", "--stats_mafcol", required=True, help="Name of the column containing the minor allele frequencies in the stats file")
ap.add_argument("-d", "--stats_infocol", required=True, help="Name of the column containing the info scores in the stats file")
ap.add_argument("-m", "--mafthres", required=True, help="Minimum minor allele frequency used for filtering variants")
ap.add_argument("-n", "--infothres", required=True, help="Minimum info score used for filtering variants")
ap.add_argument("-o", "--out", required=True, help="Path to the file to which the filtered variants are written")
ap.add_argument("-r", "--rsum", required=True, help="Path to which a summary of the results will be written")

args = vars(ap.parse_args())

base_variant_ids = set()
base_variants_count, stats_occurrence_count, removed_count = 0, 0, 0

# Read the base file to get a set of variant IDs
with gzip.open(args['base'], 'rb') as f:
    header = f.readline().decode("utf-8").split()
    base_id_col_index = header.index(args['bf_idcol'])

    while True:
        try:
            # for UKB CAD file
            # base_variant_ids.add(f.readline().decode("utf-8").split()[0].split('_')[0])

            # Add the IDs of the variants in the base file to a set
            base_variant_ids.add(f.readline().decode("utf-8").split()[base_id_col_index])
            
            # Count the amount of variants in the base file
            base_variants_count += 1
        except:
            break

base_unique_IDs_count = len(base_variant_ids)

# Read the stats file and remove IDs in the base_variant_ids set if they do not meet the provided thresholds
with gzip.open(args['stats'], 'rb') as f:

    # Get the indices of the required columns
    header = f.readline().decode("utf-8").split()
    stats_id_col_index = header.index(args['stats_idcol'])
    stats_maf_col_index = header.index(args['stats_mafcol'])
    stats_info_col_index = header.index(args['stats_infocol'])

    while True:
        try:
            line = f.readline().decode("utf-8").split()

            # Check if the stats variant occurs in the base file
            if line[stats_id_col_index] in base_variant_ids:

                # Count the amount of base file variants that occur within the stats file
                stats_occurrence_count += 1
                
                # Check if the variant meets the info and maf requirements
                if not (float(line[stats_maf_col_index]) >= float(args['mafthres']) and float(line[stats_info_col_index]) >= float(args['infothres'])):
                    base_variant_ids.remove(line[stats_id_col_index])

                    # Count the amount of times a variant was removed from the list
                    removed_count += 1
        except:
            break

# Read the base file again, this time to write the variants that occur in the base_variant_ids set to a file
with gzip.open(args['base'], 'rb') as f:

    with gzip.open(args['out'], 'wb') as o:

        # Write the header to the new file
        header = f.readline().decode("utf-8")
        o.write(header.encode())

        # Determine the index of the ID column again
        base_id_col_index = header.split().index(args['bf_idcol'])
    
        while True:
            try:
                line = f.readline().decode("utf-8")
                
                # for UKB CAD file
                # if line.split()[0].split('_')[0] in base_variant_ids:
                #     o.write(line.encode())
                
                if line.split()[base_id_col_index] in base_variant_ids:
                    o.write(line.encode())
            except:
                break

# Write a summary of the results to a file
with open(args['rsum'], 'w') as f:
    f.write(str(base_variants_count) + ' variants found in basefile\n')
    f.write('\t- unique variants: ' + str(base_unique_IDs_count) + '\n')
    f.write('\t- duplicate variants: ' + str(base_variants_count - base_unique_IDs_count) + '\n\n')
    f.write(str(stats_occurrence_count) + ' / ' + str(base_unique_IDs_count) + ' unique variants appeared in the stats file\n')
    f.write(str(removed_count) + ' / ' + str(stats_occurrence_count) + ' unique variants that appeared in the stats file did not meet the thresholds\n')
    f.write(str(stats_occurrence_count - removed_count) + ' / ' + str(stats_occurrence_count) + ' unique variants that appeared in the stats file met the thresholds\n\n')
