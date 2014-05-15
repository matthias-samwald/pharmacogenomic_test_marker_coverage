'''
Created on 06.03.2014
@author: Matthias Samwald, Medical University of Vienna
'''

import itertools
import vcf
import pprint
from collections import defaultdict

AutoVivification = lambda: defaultdict(AutoVivification) # Autovivifying dictionary, makes it easy to create nested dictionaries without running into problems of non-existing keys.

def convert_nested_dd(dd):
    '''Converts a nested defaultdict back into a native dictionary.'''
    return {k:convert_nested_dd(v) for k,v in dd.items()} if isinstance(dd, defaultdict) else dd
    
gene_definitions = AutoVivification() # captures definitions of alleles/haplotypes for several genes
assay_rsid_coverage = dict() # mapping assays to the rsids covered by each assay
assay_allele_coverage = AutoVivification() # mapping assays to the alleles covered by each assay
sample_rsid_coverage = set() # the rsids observed in the genetic sample data (1000genomes)
mapping_of_1000genomes_record_to_population = dict()

exemplary_1000genomes_record = 'HG00096' # a random pick among the 1000genome samples -- used to check which rsids are covered in the data and which are not

assay_rsid_coverage['hypothetical_assay_covering_all_rsids_in_pharmgkb'] = set() # this will be filled later on based on data from PharmGBK
assay_rsid_coverage['dmet_plus'] = set(line.strip() for line in open('data_about_assays/dmet_plus_rsids'))
assay_rsid_coverage['florida_stanford_chip'] = set(line.strip() for line in open('data_about_assays/florida_stanford_chip_rsids'))
assay_rsid_coverage['taqman'] = set(line.strip() for line in open('data_about_assays/taqman_rsids'))
assay_rsid_coverage['veracode_adme_corepanel'] = set(line.strip() for line in open('data_about_assays/veracode_adme_corepanel_rsids'))

list_of_assays = ['hypothetical_assay_covering_all_rsids_in_pharmgkb', 'florida_stanford_chip', 'dmet_plus', 'taqman', 'veracode_adme_corepanel']

haplotype_tables = dict()         
for gene in ['CYP2C9', 'CYP2C19', 'CYP2D6', 'CYP3A5', 'DPYD', 'SLCO1B1', 'TPMT', 'UGT1A1', 'VKORC1']:
    with open('gene_definitions/' + gene) as f:
        haplotype_tables[gene] = f.read()

# read haplotype table data into nested dictionary
for gene, table in haplotype_tables.items():     
    table = table.replace(" [tag]", "") # remove [tag] information for now...
    first_line = True
    for line in table.splitlines():
        if (first_line == True):
            header_array = line.split("\t")
            first_line = False
        else:
            line_array = line.split("\t")
            snps = {}
    
            rsid_start_at_column = 3    # from this column on, rsid data can be found (to the left, there is the name of the haplotype etc)
            for index, rsid in enumerate(header_array[rsid_start_at_column:]):
                if (rsid[0:2] == "rs"):    # only add if this is an rs number and not some other identifier
                    assay_rsid_coverage['hypothetical_assay_covering_all_rsids_in_pharmgkb'].add(rsid)
                    if line_array[index + rsid_start_at_column] != "":   # only add rsid:variant pair if there actually is an entry in the table
                        snps[rsid] = line_array[index + rsid_start_at_column].strip()
                
            gene_definitions[line_array[0].strip()][line_array[2].strip()] = snps # For example: {'UGT1A1': {'*60c': {'rs4148323': 'G', 'rs1057911': 'A'}}}           

gene_definitions = convert_nested_dd(gene_definitions) # convert to native dictionary

# Dump parsed gene definitions to file
gene_definitions_dump_file = open('output/gene_definitions.txt', 'w')
pprint.PrettyPrinter(indent=4, stream = gene_definitions_dump_file).pprint(gene_definitions)
gene_definitions_dump_file.close()

# read allele coverage data of various assays into nested dictionary
for gene in ['CYP2C9', 'CYP2C19', 'CYP2D6', 'CYP3A5', 'DPYD', 'SLCO1B1', 'TPMT', 'UGT1A1', 'VKORC1']:
    first_line = True
    for line in open('assay_allele_coverage/' + gene):
        if (first_line == True):
            header_array = line.split("\t")
            first_line = False
        else:
            line_array = line.split("\t")
            for i, cell in enumerate(line_array):
                if (i > 0):
                    #assay_allele_coverage.append(list([gene.strip(), line_array[0].strip(), header_array[i].strip()])) # e.g., [CYP2C9, "*2", "dmet_plus"]
                    assay_allele_coverage[header_array[i].strip()][gene.strip()][line_array[0].strip()] = cell.strip() # e.g., {'veracode_adme_corepanel': {'CYP2D6': {'*1': 'x'}}} 
 
assay_allele_coverage = convert_nested_dd(assay_allele_coverage) # convert to native dictionary

# read mapping between 1000genome records and populations
for line in open('1000genomes/20130606_sample_info.txt'):
    line_array = line.split("\t")
    mapping_of_1000genomes_record_to_population[line_array[0]] = line_array[1]    # mapping from record ID to population
    

         
# TODO: for currently_processed_gene in ['CYP2C9', 'CYP2C19', 'CYP2D6', 'CYP3A5', 'DPYD', 'SLCO1B1', 'TPMT', 'UGT1A1', 'VKORC1']:
for currently_processed_gene in ['CYP2C9', 'CYP2C19', 'CYP2D6', 'CYP3A5', 'DPYD', 'SLCO1B1', 'TPMT', 'UGT1A1', 'VKORC1']:
    
    print("Started processing data on " + currently_processed_gene)
    
    thousand_genome_samples = AutoVivification() # data about samples from 1000 genomes dataset, each sample describes one maternal and one paternal set of haplotypes
                
    f = open('output/' + currently_processed_gene + '.txt', 'w')
    output = f
                    
    vcf_reader = vcf.Reader(open('1000genomes/' + currently_processed_gene, 'r')) 
    
    for record in vcf_reader:
        sample_rsid_coverage.add(record.ID)
        if record.ID in assay_rsid_coverage['hypothetical_assay_covering_all_rsids_in_pharmgkb']:
            for call in record.samples:
                thousand_genome_samples[str(call.sample)][currently_processed_gene]['maternal_haplotype']['snps'][record.ID] = call.gt_bases.split("|")[0] # first variant is arbitrarily named 'maternal'
                thousand_genome_samples[str(call.sample)][currently_processed_gene]['paternal_haplotype']['snps'][record.ID] = call.gt_bases.split("|")[1] # second variant is arbitrarily named 'paternal'
        
                # create empty structure to hold allele information inferred later on
                for assay in list_of_assays:
                    for maternal_or_paternal_haplotype in ['maternal_haplotype', 'paternal_haplotype']:
                        thousand_genome_samples[str(call.sample)][currently_processed_gene][maternal_or_paternal_haplotype]['allele'][assay] = set()
    
    thousand_genome_samples = convert_nested_dd(thousand_genome_samples) # convert to native dictionary
    
    rsids_in_1000genomes = set(thousand_genome_samples[exemplary_1000genomes_record][currently_processed_gene]['maternal_haplotype']['snps'].keys())
    
    print(''' 

Task 1:
For different lists of SNPs covered by different assays, group alleles into indistinguishable allels.    
TODO: Compare that to alleles claimed to be discoverable by manufacturer.

''', file = output)
        
    for assay, covered_rsids in assay_rsid_coverage.items():
        print("\n\n** ASSAY:", assay, "**", file = output)
        haplotypes = gene_definitions[currently_processed_gene]
        
        #print("GENE:", gene)
        count_overlapping = 0
        overlapping_haplotype_sets = list()
        # generate 'pruned' versions of the haplotypes that only contain rsids/SNPs covered by the specific assay (i.e., the intersection of rsids in the PharmGKB table and rsids covered by the assay)
        pruned_haplotypes = dict()
        for haplotype_id, haplotype in haplotypes.items():
            pruned_haplotypes[haplotype_id] = {rsid: haplotype[rsid] for rsid in set(covered_rsids & haplotype.keys())}
            
        # check all combinations of these pruned haplotype definitions for overlaps
        for haplotype_id_1, haplotype_id_2 in itertools.combinations(pruned_haplotypes, 2):   
            haplotype_1_snps = set(pruned_haplotypes[haplotype_id_1].items())
            haplotype_2_snps = set(pruned_haplotypes[haplotype_id_2].items()) 
            if haplotype_1_snps.issubset(haplotype_2_snps):
                #count_overlapping += 1
                found_existing_set_containing_at_least_one_of_the_overlapping_snps = False
                
                # iterate through list of overlapping haplotypes already found, if one haplotype_id in the overlapping pair is already part of a set of overlapping haplotypes, add the pair to this set
                for i, val in enumerate(overlapping_haplotype_sets):
                    if (haplotype_id_1 in val) or (haplotype_id_2 in val):
                        overlapping_haplotype_sets[i].add(haplotype_id_1)
                        overlapping_haplotype_sets[i].add(haplotype_id_2)
                        found_existing_set_containing_at_least_one_of_the_overlapping_snps = True
                if not found_existing_set_containing_at_least_one_of_the_overlapping_snps:
                    overlapping_haplotype_sets.append({haplotype_id_1, haplotype_id_2})
                    
        print(len(overlapping_haplotype_sets), "set(s) of mutually indistinguishable haplotypes: ", overlapping_haplotype_sets, file = output)
        print(len(pruned_haplotypes[haplotype_id_1].keys()), "rsids were used to define haplotypes for this gene", file = output)
        
        
        # disabled the code below -- until now, we did not constraint haplotypes by presence in 1000genomes as well -- but now we do.
        '''
        rsids_covered_by_assay_and_found_in_samples = set(pruned_haplotypes[haplotype_id_1].keys()).intersection(rsids_in_1000genomes)
        print(len(rsids_covered_by_assay_and_found_in_samples), "out of these rsids are present in the 1000genomes samples:", rsids_covered_by_assay_and_found_in_samples, file = output)
        rsids_covered_by_assay_but_not_found_in_samples = set(pruned_haplotypes[haplotype_id_1].keys()).difference(rsids_in_1000genomes)   
        print(len(rsids_covered_by_assay_but_not_found_in_samples), "out of these rsids are not present in the 1000genomes samples:", rsids_covered_by_assay_but_not_found_in_samples, "\n", file = output)
        '''
            
    print(''' 

Task 2:
For different lists of SNPs covered by different assays, calculate how many samples in 1000genomes would be assigned an allele 
that actually does not match the 'correct' allele ('correct' allele = the allele inferred based on 'hypothetical_assay_covering_all_rsids_in_pharmgkb'). 
TODO: Also generate statistics on inferred alleles, ethnicities etc.

''', file = output)
    
    # print(gene_definitions['VKORC1']['*1'])
    # print(thousand_genome_samples['HG00096']['VKORC1']['paternal_haplotype']['snps'])
    
    # for each type of assay
    for assay, covered_rsids in assay_rsid_coverage.items():
        # iterate through all 1000genome samples
        for sample, sample_data in thousand_genome_samples.items():
            # iterate through genes of each sample
            for gene, gene_data in sample_data.items():
                # generate 'pruned' versions of the PharmGKB haplotypes based on the limited set of rsids that can be obseverd by the assay
                # TODO: the current implemntation is a bit inefficient, since we are redundantly generating pruned haplotype definitions for each new sample...
                pruned_haplotypes = dict()
                for haplotype_id, haplotype in gene_definitions[gene].items():
                    pruned_haplotypes[haplotype_id] = {rsid: haplotype[rsid] for rsid in set(covered_rsids & haplotype.keys() & gene_data['maternal_haplotype']['snps'].keys())}
                # skip and print message if assay has no rsids of this gene at all (which probably means that the assay was not designed to test for this gene at all)
                if len(pruned_haplotypes) == 0:
                    print("It seems like this assay covers no markers of this gene at all; skipping...", file = output)
                    continue
                # print(pruned_haplotypes)
                for maternal_or_paternal_haplotype in ['maternal_haplotype', 'paternal_haplotype']:
                    # iterate through all the pruned haplotypes and see if they are subsets of the haplotype found in the sample
                    for pruned_haplotype_id in pruned_haplotypes:
                        if set(pruned_haplotypes[pruned_haplotype_id].items()).issubset(set(gene_data[maternal_or_paternal_haplotype]['snps'].items())):
                            assay_allele_coverage_status = assay_allele_coverage.get(assay, {}).get(gene, {}).get(pruned_haplotype_id, "")
                            if (assay_allele_coverage_status != "") or assay == 'hypothetical_assay_covering_all_rsids_in_pharmgkb':   # only call allele if it is among the callable alleles for this assay, or if we are examining the hypothetical assay
                                #print(maternal_or_paternal_haplotype, "match:", pruned_haplotype_id, file = output)
                                thousand_genome_samples[sample][gene][maternal_or_paternal_haplotype]['allele'][assay].add(pruned_haplotype_id)
    
                                
                                # print('Information on coverage of this allele by this assay:', assay_allele_coverage_status, file = output)
                                # iterate through the variants expected in the full (unpruned) definition of the matching haplotype, and flag any mismatches with the data in the sample
                                '''
                                for rsid, variant in gene_definitions[gene][pruned_haplotype_id].items():
                                    if gene_data[maternal_or_paternal_haplotype].get(rsid, False) and gene_data[maternal_or_paternal_haplotype][rsid] != variant:
                                        print("Mismatch found for", rsid, "- based on haplotype definition expected", variant, "but found", gene_data[maternal_or_paternal_haplotype][rsid], file = output)                  
                              
                        else:
                            if pruned_haplotype_id[0:3] == "*1":
                                print("Did not match wild-type, here are the allele definition for this assay and the data we see in the sample:")
                                print(pruned_haplotypes[pruned_haplotype_id])
                                print(gene_data[maternal_or_paternal_haplotype]['snps'])
                                  '''
    print('''
Task 3:

Generate statistics

''', file = output)
                                
    snp_frequency_data = AutoVivification()  #    for example: {  'rs1057910': {   'count_in_gene_definitions': {'A': 32, 'C': 2},
                                            #                                     'count_in_population_sample': {'A': 2091, 'C': 93}}}
    for haplotype_id in gene_definitions[currently_processed_gene]:
        for rsid in gene_definitions[currently_processed_gene][haplotype_id]:
            # this is quite unelegant because of the way the autovivification dictionary works (inexistant values do not default to 0)
            if isinstance(snp_frequency_data[rsid]['count_in_gene_definitions'][gene_definitions[currently_processed_gene][haplotype_id][rsid]], int):
                snp_frequency_data[rsid]['count_in_gene_definitions'][gene_definitions[currently_processed_gene][haplotype_id][rsid]] += 1
            else:
                snp_frequency_data[rsid]['count_in_gene_definitions'][gene_definitions[currently_processed_gene][haplotype_id][rsid]] = 1
             
    for sample in thousand_genome_samples:
        for gene in thousand_genome_samples[sample]:
            for maternal_or_paternal_haplotype in ['maternal_haplotype', 'paternal_haplotype']:
                for rsid in thousand_genome_samples[sample][gene][maternal_or_paternal_haplotype]['snps']:
                    if isinstance(snp_frequency_data[rsid]['count_in_population_sample'][thousand_genome_samples[sample][gene][maternal_or_paternal_haplotype]['snps'][rsid]], int):
                        snp_frequency_data[rsid]['count_in_population_sample'][thousand_genome_samples[sample][gene][maternal_or_paternal_haplotype]['snps'][rsid]] += 1
                    else:
                        snp_frequency_data[rsid]['count_in_population_sample'][thousand_genome_samples[sample][gene][maternal_or_paternal_haplotype]['snps'][rsid]] = 1
    
    snp_frequency_data = convert_nested_dd(snp_frequency_data)                                                                                                        
    
    
    '''
    print("The following SNP variants are used in a gene definition, but were never observed in a sample, even though their rs number was tested for at least one sample (this can also be a hint for strand orientation mismatch):", file = output)
    for rsid in assay_rsid_coverage['hypothetical_assay_covering_all_rsids_in_pharmgkb']:
        for snp_variant in snp_frequency_data[rsid]['count_in_gene_definitions']:
            if snp_frequency_data[rsid].get('count_in_population_sample'):
                if snp_frequency_data[rsid].get('count_in_population_sample').get(snp_variant, 0) == 0:
                    print(rsid, snp_variant, file = output)
    print("---------")
    '''
    
    
    # show alleles called 
    print('''
    
    SNP variant counts in gene definitions and in sample data.
    
    ''', file = output)
    
    # display header row
    output_line = "population" + "\t" + "sample" + "\t" + "gene" + "\t" + "maternal_or_paternal_haplotype"
    for assay in list_of_assays:
        output_line = output_line + "\t" + assay + "\t Bogus calls from " + assay 
    print(output_line, file = output)    
    
    # display data     
    for sample in thousand_genome_samples:
            for maternal_or_paternal_haplotype in ['maternal_haplotype', 'paternal_haplotype']:
                output_line = mapping_of_1000genomes_record_to_population[sample] + "\t" + sample + "\t" + currently_processed_gene + "\t" + maternal_or_paternal_haplotype
                for assay in list_of_assays:
                    output_line = output_line + "\t" + str(thousand_genome_samples[sample][currently_processed_gene][maternal_or_paternal_haplotype]['allele'][assay])
                    output_line = output_line + "\t Bogus calls: " + str(thousand_genome_samples[sample][currently_processed_gene][maternal_or_paternal_haplotype]['allele'][assay].difference(thousand_genome_samples[sample][currently_processed_gene][maternal_or_paternal_haplotype]['allele']['hypothetical_assay_covering_all_rsids_in_pharmgkb']))
                print(output_line, file = output)
                #for assay in thousand_genome_samples[sample][gene][maternal_or_paternal_haplotype]['allele']:
                # print(sample, gene, maternal_or_paternal_haplotype, assay, thousand_genome_samples[sample][gene][maternal_or_paternal_haplotype]['allele'][assay])
    
    pp = pprint.PrettyPrinter(indent=4, stream = output)
    pp.pprint(snp_frequency_data)
    # pp.pprint(thousand_genome_samples)

    f.close()                                


with open('output/overview.txt', 'w') as f:
    print("Started generating overview document")
    print(len(assay_rsid_coverage['hypothetical_assay_covering_all_rsids_in_pharmgkb']), 
          "rsids in hypothetical_assay_covering_all_rsids_in_pharmgkb:", file = f)
    print(assay_rsid_coverage['hypothetical_assay_covering_all_rsids_in_pharmgkb'], file = f)
    
    rsids_in_pharmgkb_but_not_in_1000genomes = set(assay_rsid_coverage['hypothetical_assay_covering_all_rsids_in_pharmgkb']).difference(sample_rsid_coverage)
    print("\n", len(rsids_in_pharmgkb_but_not_in_1000genomes), 
          "rsids in hypothetical_assay_covering_all_rsids_in_pharmgkb that are not found in 1000genomes data:", file = f)
    print(rsids_in_pharmgkb_but_not_in_1000genomes, file = f)
    
   
    
'''
rsids_covered_by_assay_but_not_found_in_samples = covered_rsids - set(thousand_genome_samples[exemplary_1000genomes_record][currently_processed_gene]['maternal_haplotype']['snps'].keys())    
print("\nThe following", len(rsids_covered_by_assay_but_not_found_in_samples), "rsids are covered by this assay, but are not present in the 1000genomes samples:", rsids_covered_by_assay_but_not_found_in_samples, file = output)
'''        