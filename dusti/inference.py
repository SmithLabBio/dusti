#import sys # for taking arguments
import itertools # for getting all combinations of quartets from species list
import random # for sampling quartets
import os # for getting lists of files
import dendropy # for reading alignments
#import copy # for making subalignments
import numpy as np # for handling arrays
from collections import Counter # for site pattern count dictionaries
import argparse # to parse user arguments
from scipy.linalg import svd
import math
import time # to track time
#from operator import itemgetter
import logging # for logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

def check_input(input_folder):

    logging.info(f"Conducting quartet inference for alignments in folder: {input_folder}")

    phylips = [x for x in os.listdir(input_folder) if x.endswith('.phy')]
    fastas = [x for x in os.listdir(input_folder) if x.endswith('.fa') or x.endswith('.fasta')]

    if len(phylips) == 0 and len(fastas) == 0:
        logging.error("ERROR: There are no phylip (ends with .phy) or fasta (ends with .fa or .fasta) files detected in the input directory.")
        exit(1)
    
    elif len(phylips) > 0 and len(fastas) > 0:
        logging.error("ERROR: Both phylip and fasta files are detected in the input directory. All alignments must be in the same format.")
        exit(1)
    
    elif len(phylips) > 0:
        logging.info(f"Using {len(phylips)} alignments for inference.")
        return "phylip"
    
    else:
        logging.info(f"Using {len(fastas)} alignments for inference.")
        return "fasta"

def get_species_map(map_file):

    sp_dict = {}
    match_type = []

    with open(map_file, 'r') as f:
        for line in f.readlines():
            if len(line.split()) != 2:
                logging.error("ERROR: The mapping file must contain two columns separated by spaces or tabs.")
                exit(1)
            sp_name = line.split()[0]
            gene_name = line.split()[1]
            if gene_name.endswith('*'):
                match_type.append('prefix')
                if gene_name.startswith('*'):
                    logging.error("ERROR: dusti cannot match both prefixes and suffixes in the mapping file.")
            elif gene_name.startswith('*'):
                match_type.append('suffix')
            else:
                match_type.append('exact')

            if sp_name not in sp_dict:
                sp_dict[sp_name] = gene_name.strip('*')
            else:
                sp_dict[sp_name].append(gene_name.strip('*'))


    if len(set(match_type)) > 1:
        logging.error("ERROR: Mixed matching used in mapping file. Please choose either exact matching, prefix matching, or suffix matching.")


    return(sp_dict, match_type[0])

def sample_quartets(map_dict, max_quartets, seed):

    logging.info(f"Sampling up to {max_quartets} quartets.")

    # get all combinations
    combos = list(itertools.combinations(map_dict.keys(), 4))

    # if combos is longer than max_quartets, sample max_quartets combos from combos
    if len(combos) > max_quartets:
        if seed is not None:
            random.seed(seed)
            logging.info(f"Random seed set to: {seed}")
        combos = random.sample(combos, max_quartets)
    return(combos)

def to_matrix(alignment):
    sequences = [str(alignment[taxon]) for taxon in alignment]
    taxon_list = [taxon for taxon in alignment]
    alignment_matrix = np.array([list(sequence) for sequence in sequences])
    return(alignment_matrix, taxon_list)

def count_sites_alignment(alignment):

    # get site strings
    site_strings = [tuple(column) for column in alignment.T]
    
    # build the dictionary
    site_dict = Counter(site_strings)

    return(site_dict)

def subset_dictionary(quartet, taxon_names, species_matches, alignment_dictionary):

    # find all quartet combos
    species_list = [[],[],[],[]]

    for species_name in range(len(quartet)):
        for idname in range(len(taxon_names)):
            if species_matches[idname] == quartet[species_name]:
                species_list[species_name].append(taxon_names[idname])

    # generate combos of the species list
    combinations = list(itertools.product(*species_list))

    # for each combination, get the site patterns
    list_of_combo_sitedicts = []
    for combo in combinations:
        
        # get indices
        index_0 = taxon_names.index(combo[0])
        index_1 = taxon_names.index(combo[1])
        index_2 = taxon_names.index(combo[2])
        index_3 = taxon_names.index(combo[3])

        # Step 1: Create a new dictionary with modified keys
        modified_dict = {}

        for old_key, value in alignment_dictionary.items():
            new_key = (old_key[index_0], old_key[index_1], old_key[index_2], old_key[index_3])
            modified_dict[new_key] = modified_dict.get(new_key, 0) + value

        list_of_combo_sitedicts.append(modified_dict)
    
    # average across combinations for this quartet
    sums = Counter()
    for itemset in list_of_combo_sitedicts:
        sums.update(itemset)

    # average site counts over the number of combos for these taxa
    average_site_patterns = {x: float(sums[x])/len(list_of_combo_sitedicts) for x in sums.keys()}

    return(average_site_patterns)

def calculate_site_pattern_counts(sampled_quartets, input_folder, map_dict, match_type, file_type):

    logging.info("Calculating site pattern counts.")

    # inverse dictionary
    ivd = {v: k for k, v in map_dict.items()}


    # get list of alignments in folder
    alignments = os.listdir(input_folder)
    if file_type=="phylip":
        alignments = [x for x in alignments if x.endswith('.phy')]
    else:
        alignments = [x for x in alignments if x.endswith('.fasta') or x.endswith('.fa')]

    dendropy_matrices = []
    dendropy_taxa = []

    # create quartet dictionary
    quartet_site_dictionary = {quartet: list() for quartet in sampled_quartets}

    for alignment in alignments:
        
        # read alignment, and translate to matrix
        dendropy_alignment = dendropy.DnaCharacterMatrix.get(path=os.path.join(input_folder,alignment), schema=file_type)
        dendropy_matrix, dendropy_taxon = to_matrix(dendropy_alignment)

        # change taxon names to match species names
        taxon_names = [str(x).strip("'") for x in dendropy_taxon]
        if match_type == "prefix":
            species_matches = [ivd[x.split('_')[0]] for x in taxon_names]
        elif match_type == "suffix":
            species_matches = [ivd[x.split('_')[-1]] for x in taxon_names]
        else:
            species_matches = [ivd[x] for x in taxon_names]

        #dendropy_matrices.append(dendropy_matrix)
        #dendropy_taxa.append(dendropy_taxon)

        # build dictionary of site patterns for alignment
        alignment_dictionary = count_sites_alignment(dendropy_matrix)

        # iterate over quartets and subset dictionary
        for quartet in sampled_quartets:

            temp_quartet_dictionary = subset_dictionary(quartet, taxon_names, species_matches, alignment_dictionary)
            quartet_site_dictionary[quartet].append(temp_quartet_dictionary)
            del temp_quartet_dictionary

    # sum across alignments for a quartet:
    summed_quartet_site_dictionary = dict.fromkeys(sampled_quartets, [])
    for key, value in quartet_site_dictionary.items():
        sums = Counter()
        for itemset in value:
            sums.update(itemset)
        summed_quartet_site_dictionary[key] = sums
        del sums

    return(summed_quartet_site_dictionary)

def find_svdq(quartet_site_dictionary):

    logging.info("Finding quartets with lowest SVD score.")

    quartet_trees = []

    for i in quartet_site_dictionary:

        # set up relationships (list of tuples)
        pairs = list(itertools.combinations(i, 2))[0:3]
        
        flattening_matrices = []

        for item in pairs:

            # get the complement (these are 3,4)
            complement = tuple([x for x in i if not x in item])
            
            # get indices
            index_1 = i.index(item[0])
            index_2 = i.index(item[1])
            index_3 = i.index(complement[0])
            index_4 = i.index(complement[1])
            
            # get combos for iteration
            row_combos = list(itertools.product('ACGT', repeat=2))
            col_combos = list(itertools.product('ACGT', repeat=2))

            # create empty matrix to populate
            this_matrix = np.empty((16, 16), dtype=float)  # You can change the dtype if needed

            row_count = -1
            for j in row_combos:
                row_count+=1
                col_count = -1
                for k in col_combos:
                    col_count+=1
                    list_for_search = [None,None,None,None]
                    list_for_search[index_1] = j[0]
                    list_for_search[index_2] = j[1]
                    list_for_search[index_3] = k[0]
                    list_for_search[index_4] = k[1]
                    this_cell = quartet_site_dictionary[i][tuple(list_for_search)]
                    this_matrix[row_count, col_count] = this_cell

            flattening_matrices.append(this_matrix)

        svdscores = []
        for l in flattening_matrices:
            U, s, V = svd(l)
            svds = s[10:]
            svdscore_current = math.sqrt(sum([num ** 2 for num in svds]))
            svdscores.append(svdscore_current)


        # select matrix relationship with lowest rank as quartet relationship
        index_of_lowest = svdscores.index(min(svdscores))
        relationship_of_lowest = pairs[index_of_lowest]
        complement_of_lowest = tuple([x for x in i if not x in relationship_of_lowest])
        tree = '((%s,%s),(%s,%s));' % (relationship_of_lowest[0], relationship_of_lowest[1], complement_of_lowest[0], complement_of_lowest[1])
        quartet_trees.append(tree)

    return(quartet_trees)

def find_parsimony(quartet_site_dictionary):

    logging.info("Finding most parsimonious quartets.")

    quartet_trees = []

    for i in quartet_site_dictionary:

        # set up relationships (list of tuples)
        pairs = list(itertools.combinations(i, 2))[0:3]

        parsimonies = {}

        for item in pairs:

            # get the complement (these are 3,4)
            complement = tuple([x for x in i if not x in item])
            
            # get indices
            index_1 = i.index(item[0])
            index_2 = i.index(item[1])
            index_3 = i.index(complement[0])
            index_4 = i.index(complement[1])

            this_total = 0
            for key in quartet_site_dictionary[i].keys():
                if (key[index_1] == key[index_2]) and (key[index_3] == key[index_4]) and (key[index_1] != key[index_3]):
                    this_total += quartet_site_dictionary[i][key]
            
            parsimonies[item] = this_total


        ## select highest count as most parsimonious
        max_key = max(parsimonies, key=lambda k: parsimonies[k])
        complement_of_max = tuple([x for x in i if not x in max_key])
        tree = '((%s,%s),(%s,%s));' % (max_key[0], max_key[1], complement_of_max[0], complement_of_max[1])
        quartet_trees.append(tree)

    return(quartet_trees)


def quartet_puzzling(svdq_quartet_trees, output_directory, qfm, prefix):

    # make output file
    outfile = open(os.path.join(output_directory, f"{prefix}quartets.trees"), 'w')
    for item in svdq_quartet_trees:
        outfile.write(item)
        outfile.write (' 1')
        outfile.write('\n')
    outfile.close()
    
    command = 'java -jar %s -i %s -o %s' % (qfm, os.path.join(output_directory, f"{prefix}quartets.trees"), os.path.join(output_directory, f"{prefix}.tre"))
    os.system(command)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Quartet MaxCut Python Script')

    # Define command-line arguments
    parser.add_argument('-q', '--max_quartets', type=int, default=np.inf, help='Maximum number of quartets')
    parser.add_argument('-a', '--map', type=str, help='Path to the species to gene mapping file.')
    parser.add_argument('-o', '--output_directory', type=str, default='dusti_results', help='Output folder directory')
    parser.add_argument('--qfm', type=str, help='Path to wQFM jar file.')
    parser.add_argument('-i', '--input_folder', type=str, help='Path to the input folder containing gene alignments in phylip format (will analyze all files in folder ending in .phy)')
    parser.add_argument('-s', '--seed', type=int, help='Random seed for reproducibility')
    parser.add_argument('--force', action='store_true', help='Allow overwriting of the output folder (default: False)')
    parser.add_argument('--svd', action='store_true', help="Perform inference based on SVD.")
    parser.add_argument('--parsimony', action='store_true', help="Perform inference based on parsimony.")

    # Parse the arguments
    args = parser.parse_args()

    # create output directory
    if os.path.exists(args.output_directory):
        if args.force:
            logging.warning(f"Output folder '{args.output_directory}' exists and will be overwritten.")
            os.system(f"rm -r {args.output_directory}")
            os.mkdir(args.output_directory)
        else:
            logging.error(f"Output folder '{args.output_directory}' already exists. Use --force to overwrite.")
            exit(1)
    else:
        os.mkdir(args.output_directory)

    return args

def main():

    start_time = time.time()

    # set parameters
    args = parse_arguments()

    # check input folder contents
    file_type = check_input(args.input_folder)

    # species map
    map_dict, match_type = get_species_map(args.map)

    # sample quartets
    sampled_quartets = sample_quartets(map_dict = map_dict, max_quartets = args.max_quartets, seed = args.seed)

    ## count site patterns
    quartet_site_dictionary = calculate_site_pattern_counts(sampled_quartets = sampled_quartets, input_folder = args.input_folder, match_type = match_type, map_dict = map_dict, file_type=file_type)

    if args.svd:
        # find best (svdq)
        svdq_quartet_trees = find_svdq(quartet_site_dictionary)

        # do quartet puzzling
        quartet_puzzling(svdq_quartet_trees, args.output_directory, args.qfm, 'svd')

    if args.parsimony:
        # find best (svdq)
        parsimony_quartet_trees = find_parsimony(quartet_site_dictionary)

        # do quartet puzzling
        quartet_puzzling(parsimony_quartet_trees, args.output_directory, args.qfm, 'parsimony')

    # time
    end_time = time.time()
    elapsed_time = end_time - start_time
    logging.info(f"Elapsed Time: {elapsed_time} seconds")

if __name__ == "__main__":
    main()