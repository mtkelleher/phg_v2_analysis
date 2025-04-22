import gzip
import os
import glob
import argparse

def parse_alt_allele_field(line):
    # Remove the alt‑allele prefix, trailing '\n', and the closing '>'
    line = line.removeprefix('##ALT=<').rstrip('\n').removesuffix('>')

    # Store the expected info keys in alt-allele field
    expected_keys = [
        'ID', 'Description', 'Source', 'SampleName',
        'Regions', 'Checksum', 'RefChecksum', 'RefRange'
    ]

    # Start alt-allele field dict
    alt_allele_field_dict = {}

    # Split by '='
    first_split = line.split('=')

    # Raise error if an unexpected extra '=' or missing key
    if len(first_split) != 9:
        raise ValueError(
            'This parser assumes one "=" per key‑value pair in the ALT header. '
            'An unexpected "=" was found in a value '
            '(probably in Description, Source, SampleName, Regions or RefRange).'
        )

    # Set the current key to the first thing in split
    current_key = first_split[0]

    # Iterate through info from second split to second to last
    for info in first_split[1:-1]:
        # Split by ','
        second_split = info.split(',')

        # Next key is the last value from this split
        next_key = second_split[-1]

        # Store info up to the length of the next key and comma
        alt_allele_field_dict[current_key] = info[:-len(next_key) - 1]

        # Set current key to next key
        current_key = next_key

    # Add last one
    alt_allele_field_dict[current_key] = first_split[-1]

    # Check that expected_keys were used
    for expected_key in expected_keys:
        if expected_key not in alt_allele_field_dict:
            print('Missing expected key', expected_key, 'in alt-allele', line)

    return alt_allele_field_dict

def parse_data_lines(line):
    # Split line by \t
    split = line.strip().split('\t')

    # Use split list to populate data_lines_dict
    data_lines_dict = {
        'CHROM':split[0],
        'POS':split[1],
        'ID':split[2],
        'REF':split[3],
        'ALT':split[4],
        'QUAL':split[5],
        'FILTER':split[6],
        'INFO':split[7],
        'FORMAT':split[8],
        'ALLELE_VALUE':split[9]
    }
    
    return data_lines_dict

def process_reference_ranges_file(ref_ranges_file):
    # Start the reference ranges list
    reference_ranges_list = []
    
    with open(ref_ranges_file, 'r') as reference_ranges:
        # Iterate through each line of the reference ranges file formatted as Sequence name, Start position (bp), End position (bp), Range ID, Score (always 0), Strand information
        for line in reference_ranges:
            # Split the line formatted 
            split = line.strip().split('\t')
            # Get reference range contig
            reference_range_contig = str(split[0])
            # Get reference range start position and convert from 0-based bed format to vcf 1-based format
            reference_range_start = int(split[1]) + 1
            # Get reference range end position
            reference_range_end = int(split[2])
    
            # Store reference ranges in a list as a tuple of contig, start, end positions
            reference_ranges_list.append((reference_range_contig, reference_range_start, reference_range_end))

    # Print statement
    print('Number of reference ranges in ' + ref_ranges_file + ': ' + str(len(reference_ranges_list)))
    
    return reference_ranges_list

def process_merged_parents_vcf(merged_parents_vcf_path, reference_ranges_list):
    # Open vcf file, handling gzip compression
    open_func = gzip.open if merged_parents_vcf_path.endswith('.gz') else open
    with open_func(merged_parents_vcf_path, 'rt') as merged_parents_vcf:  # 'rt' ensures text mode for gzip
        # Print statement
        print('Processing:', merged_parents_vcf_path)
        
        # Set data_lines to False
        data_lines = False

        # Set num_variants to 0
        num_variants = 0
    
        # Start empty header matrix
        header_matrix = []
    
        # Start empty reference range data lines info dictionary and set current_reference_range_range to None
        reference_range_data_lines = {}
        current_reference_range = (None, None, None)
        reference_ranges_in_vcf_ordered_list = []
    
        # Iterate through each line of the merged_parents_vcf
        for line in merged_parents_vcf:
            # If not a data line
            if not data_lines:
                # Store header lines not starting with '#CHROM'
                if not line.startswith('#CHROM'):
                    header_matrix.append(line.strip())
                
                # Is line that starts with #CHROM
                else:
                    # Set data_lines to True
                    data_lines = True
    
                    # Split the tab separated line
                    split = line.strip().split('\t')
    
                    # Set header_line as empty string and matrix_column_to_parent_key to an empty dict
                    matrix_column_to_parent_key = {}
                    header_line = ''
    
                    # Iterate through fields in split
                    for i in range(len(split)):
                        # If fixed field 0
                        if i == 0:
                            header_line = split[0]
                        # If a fixed field 1 - 8
                        elif 1 <= i <= 8:
                            header_line += '\t' + split[i]
                        # Else it is a genotype field
                        else:
                            # Get parent in column and set matrix_column_to_parent_key[i] equal to parent
                            parent = split[i]
                            matrix_column_to_parent_key[i] = parent
    
            # Else it is a data line
            else:
                # Add 1 var to num_variants
                num_variants += 1

                # Split the tab separated line
                split = line.strip().split('\t')
    
                contig = str(split[0])
                pos = int(split[1])
    
                # If the contig is not the same, or the pos is not greater than or equal to current reference range start and less than or equal to current reference range end
                if contig != current_reference_range[0] or not current_reference_range[1] <= pos <= current_reference_range[2]:
                    # Find reference range that holds the position and set as current reference range
                    for reference_range in reference_ranges_list:
                        # If the contig is the same, and the pos is greater than or equal to current reference range start and less than or equal to current reference range end
                        if contig == reference_range[0] and reference_range[1] <= pos <= reference_range[2]:
                            # Set current_reference_range
                            current_reference_range = reference_range
        
                            # Append the current reference range to reference_ranges_in_vcf_ordered_list
                            reference_ranges_in_vcf_ordered_list.append(current_reference_range)
    
                            # Create a dictionary for the reference range in the reference_range_data_lines dictionary
                            reference_range_data_lines[current_reference_range] = {}
    
                            # Create a list for fixed_fields in the reference range dictionary within the reference_range_data_lines dictionary
                            reference_range_data_lines[current_reference_range]['fixed_fields'] = []
    
                            # For each genotype field column create a list in the reference range dictionary within the reference_range_data_lines dictionary
                            for column in matrix_column_to_parent_key.keys():
                                reference_range_data_lines[current_reference_range][column] = []
    
                            # Stop search for current_reference_range
                            break
    
                # Iterate through fields in split
                for i in range(len(split)):
                    # If fixed field 0
                    if i == 0:
                        # Start fixed fields string
                        fixed_fields = str(split[i])
                    # If a fixed field 1 - 8
                    elif 1 <= i <= 8:
                        # Add to fixed fields string
                        fixed_fields += '\t' + str(split[i])
                    # Else it is a genotype field
                    else:
                        # Get genotype for genotype filed and store append to the list for the genotype
                        genotype = str(split[i])
                        reference_range_data_lines[current_reference_range][i].append(genotype)
                # Append the full fixed fields string to the list for that reference range
                reference_range_data_lines[current_reference_range]['fixed_fields'].append(fixed_fields)

    # Print statements
    print('Number of reference ranges in ' + merged_parents_vcf_path + ': ' + str(len(reference_ranges_in_vcf_ordered_list)))
    print('Number of variants in ' + merged_parents_vcf_path + ': ' + str(num_variants))

    return reference_ranges_in_vcf_ordered_list, matrix_column_to_parent_key, reference_range_data_lines, header_matrix, header_line

def get_imputed_haplotype_information(hvcf_directory):
    # Imputed sample alt dictionary to store the sampleName and RefRange
    imputed_samples_alt_dicts = {}

    # imputed sample list
    imputed_sample_list = []
    
    # Find all h.vcf and h.vcf.gz files in the directory
    hvcf_files = glob.glob(os.path.join(hvcf_directory, "*.h.vcf")) + glob.glob(os.path.join(hvcf_directory, "*.h.vcf.gz"))

    # Specify prefix ALT fields for identifying alt lines
    alt_id_prefix = '##ALT=<'

    for hvcf in hvcf_files:
        # Start a temporary imputed_sample_alt_dict and imputed_sample_header_dict
        imputed_sample_alt_dict = {}
        imputed_sample_header_dict = {}

        # Open vcf file, handling gzip compression
        open_func = gzip.open if hvcf.endswith('.gz') else open
        with open_func(hvcf, 'rt') as file:  # 'rt' ensures text mode for gzip
            next(file) # Skip header line which is ##fileformat=VCFv4.2

            # Print which hvcf file is being processed
            print('Processing hvcf:', hvcf)
            
            # Data lines start after info lines and after the first line starting with #CHROM
            data_lines = False
        
            for line in file:
                if data_lines == True: # If a data line
                    # If the id and ref range match then keep it, else don't
                    data_lines_dict = parse_data_lines(line)

                    ALT = data_lines_dict['ALT'].removeprefix('<').removesuffix('>') # Remove prefix and suffix around ALT
                    CHROM = data_lines_dict['CHROM']
                    POS = int(data_lines_dict['POS'])

                    # Check if there is an ALT at this reference range in the data lines to look up
                    if ALT != '.':
                        # imputed_sample_header_dict[id] = ((ref_range_contig, ref_range_start, ref_range_end), parent_sample_name)
                        current_haplotype_in_header = imputed_sample_header_dict[ALT]

                        # Only store info if current haplotype's CHROM and POS are equal to the reference range in the alt header
                        if current_haplotype_in_header[0][0] == CHROM and current_haplotype_in_header[0][1] == POS:
                            # dict needs the key to be dict[(contig, start, end)] = parent_sample_name
                            imputed_sample_alt_dict[current_haplotype_in_header[0]] = current_haplotype_in_header[1]
        
                else: # If not a data line, these are the info lines which are before the data lines
                    if line.startswith(alt_id_prefix): # If an ALT line
                        # Extract info from the ALT line
                        alt_allele_field_dict = parse_alt_allele_field(line)
                        # RefRange
                        ref_range = alt_allele_field_dict['RefRange']
                        split_contig_ref_range = ref_range.split(':')
                        ref_range_contig = split_contig_ref_range[0]
                        split_pos_ref_range = split_contig_ref_range[1].split('-')
                        ref_range_start = int(split_pos_ref_range[0])
                        ref_range_end = int(split_pos_ref_range[1])
                        # SampleName
                        parent_sample_name = alt_allele_field_dict['SampleName']
                        # ID
                        id = alt_allele_field_dict['ID']

                        # In header dict, for ID get reference range contig, start, end, and parent sample name
                        imputed_sample_header_dict[id] = ((ref_range_contig, ref_range_start, ref_range_end), parent_sample_name)

                    else:
                        # Other info lines in the file are not needed
                        if line.startswith('#CHROM'):
                            data_lines = True # The next line is the start of data lines

                            # Get the sampleName
                            split = line.strip().split('\t')
                            imputed_sample_name = split[9]

                            # Append the sampleName to imputed_sample_list
                            imputed_sample_list.append(imputed_sample_name)

            # Add temporary dict to imputed_samples_alt_dicts
            imputed_samples_alt_dicts[imputed_sample_name] = imputed_sample_alt_dict.copy()
                            
    return imputed_samples_alt_dicts, imputed_sample_list

def main():
    parser = argparse.ArgumentParser(
        description="""
        This script processes a merged parents VCF (from `phg merge-gvcf`), a reference_ranges.bed file,
        and a directory of imputed hVCFs to produce a merged imputed VCF.

        For each variant in the merged parent VCF, genotype calls are pulled from the parents and grouped
        by reference range. The script then matches these with corresponding entries in the imputed hVCFs
        by examining their ALT headers and data lines. Imputed haplotypes are only kept if the reference
        range assigned matches between the entry in the ALT headers and the data lines because we want only
        one call for a reference range per sample. If the haplotype call for a reference range is missing,
        the genotype is marked as '.'; if a sample is imputed to the reference, its genotype is marked as
        '0'.

        Required file formats:

        - Merged Parents VCF (Required)
        - Ref ranges BED (Required):
            Example: chr1	0	34116	intergenic_chr1:0-34116	0	+
        - Imputed hVCF (Required):
            Example header and data lines: 
            ##fileformat=VCFv4.2
            ##ALT=<ID=...,Description="haplotype data for line: HP301",Source="phg_v2/vcf_dbs/assemblies.agc",SampleName=HP301,Regions=chr4:1-100001,Checksum=...,RefRange=chr4:1-100000,RefChecksum=...>
            ...
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ZeaSynDH-2474
            chr1	1	.	T	<9ba10f19f09db4e1ce3d4b238d075444>	.	.	END=108053	GT	1

        Usage:
            python3 imputed_merged_vcf.py --ref_ranges_file phg_v2/output/ref_ranges.bed \\
                --merged_parents_vcf_path phg_v2/merged_parents/merged_parents.vcf \\
                --imputed_hvcf_directory phg_v2/output/read_mappings/vcf_files \\
                --reference_sample_name B73 \\
                --merged_imputed_vcf_path phg_v2/output/merged_imputed.vcf
        """,
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "--ref_ranges_file",
        required=True,
        help="Path to the ref_ranges.bed file. Required."
    )
    parser.add_argument(
        "--merged_parents_vcf_path",
        required=True,
        help="Path to the merged parents VCF. Required."
    )
    parser.add_argument(
        "--imputed_hvcf_directory",
        required=True,
        help="Path to the directory containing imputed hVCF files. Required."
    )
    parser.add_argument(
        "--reference_sample_name",
        required=False,
        help="Sample name for the reference as it appears in imputed hVCF headers. Optional."
    )
    parser.add_argument(
        "--merged_imputed_vcf_path",
        required=False,
        help="Path to output merged imputed VCF. If not provided, defaults to <merged_parents_vcf_path>/merged_imputed.vcf. Optional."
    )

    args = parser.parse_args()

    # Validate paths
    if not os.path.exists(args.ref_ranges_file):
        parser.error(f"File not found: {args.ref_ranges_file}")
    if not os.path.exists(args.merged_parents_vcf_path):
        parser.error(f"File not found: {args.merged_parents_vcf_path}")
    if not os.path.isdir(args.imputed_hvcf_directory):
        parser.error(f"Directory not found: {args.imputed_hvcf_directory}")

    # Set default output path if none provided
    merged_imputed_vcf_path = args.merged_imputed_vcf_path
    if not merged_imputed_vcf_path:
        merged_imputed_vcf_path = os.path.join(
            os.path.dirname(args.merged_parents_vcf_path), "merged_imputed.vcf"
        )

    # Function calls
    reference_ranges_list = process_reference_ranges_file(args.ref_ranges_file)
    reference_ranges_in_vcf_ordered_list, matrix_column_to_parent_key, reference_range_data_lines, header_matrix, header_line = process_merged_parents_vcf(
        args.merged_parents_vcf_path, reference_ranges_list
    )
    imputed_samples_alt_dicts, imputed_sample_list = get_imputed_haplotype_information(args.imputed_hvcf_directory)
    
    # Print statement
    print('Column to parent key:', matrix_column_to_parent_key)

    # Iterate through imputed samples
    for imputed_sample in imputed_sample_list:
        # Add the imputed sample to the header line
        header_line += '\t' + imputed_sample
    
        # Print statement
        print('Getting gentoype calls for:', imputed_sample)
    
        # Iterate through each reference range as they appear in merged_parents_vcf
        for reference_range in reference_ranges_in_vcf_ordered_list:
            # Temporary pointer to reference_range in the reference_range_data_lines dict
            reference_range_data = reference_range_data_lines[reference_range]
            # Temporary pointer to fixed fields list for reference range
            fixed_fields = reference_range_data['fixed_fields'] 
            # Get imputed parent for imputed sample at this reference range, if no haplotype for this imputed parent set the imputed_parent to None
            imputed_parent = imputed_samples_alt_dicts[imputed_sample].get(reference_range, None)
    
    
            # If the reference_sample_name is set, check if the imputed parent matches the reference
            if args.reference_sample_name is not None and imputed_parent == args.reference_sample_name:
                # Iterate through fixed fields for reference range
                for i, val in enumerate(fixed_fields):
                    # If the imputed parent is the reference, add 0 genotype calls
                    fixed_fields[i] = val + '\t0'
    
            elif imputed_parent == None:
                # Iterate through fixed fields for reference range
                for i, val in enumerate(fixed_fields):
                    # If the imputed parent is None, add . genotype calls
                    fixed_fields[i] = val + '\t.'
            
            # The imputed parent is not the reference
            else:
                # Search for the column in matrix_column_to_parent_key for this parent
                for column, parent in matrix_column_to_parent_key.items():
                    if imputed_parent == parent:
                        current_column = column
                        # Temporarily store the genotype fields for this imputed parent within this reference range
                        genotype_fields = reference_range_data[current_column]
            
                        # Stop searching for the column number
                        break
                    
                # Else parent is not found
                else:
                    raise ValueError(
                        f"The imputed parent '{imputed_parent}' from the hVCF was not found in the merged parents VCF. "
                        f"No matching column exists in matrix_column_to_parent_key."
                    )
    
                # Iterate through fixed fields for reference range
                for i, val in enumerate(fixed_fields):
                    # Add the genotype fields from the imputed parent to fixed fields for the reference range
                    fixed_fields[i] = val + '\t' + genotype_fields[i]

        # Clear memory from that imputed sample in the imputed_samples_alt_dicts
        del imputed_samples_alt_dicts[imputed_sample]
    
    # Write the data to output
    with open(merged_imputed_vcf_path, 'w') as merged_imputed_vcf:
        # Write header matrix
        for line in header_matrix:
            merged_imputed_vcf.writelines(line + '\n')
    
        # Write header line
        merged_imputed_vcf.writelines(header_line + '\n')
        
        # Iterate through each reference range as they appear in merged_parents_vcf
        for reference_range in reference_ranges_in_vcf_ordered_list:
            for line in reference_range_data_lines[reference_range]['fixed_fields']:
                merged_imputed_vcf.write(line + '\n')
    
        # Print statement
        print('Merged imputed vcf written to:', merged_imputed_vcf_path)

    return

if __name__ == "__main__":
    main()
