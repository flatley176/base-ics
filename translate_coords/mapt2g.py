import argparse
import csv
import pandas
import re

def split_cigar_into_tuples(cigar_string):
    """
    splits a cigar string into a list of easily parseable tuples
   
    input param cigar_string: a valid cigar string
    
    output: list of tuples; each tuple is a number representing the segment length
            (will need format conversion) and either of M, I or D.
    """
    return re.findall('([0-9]+)([MID])', cigar_string)


def update_coords(last_known_genomic_pos, 
                  last_known_transcript_pos, 
                  segment_length, 
                  condition):
    """
    updates 0-based genomic/transcript coordinates for condition

    input param last_known_genomic_pos: 0-based genomic start position for segment
    input param last_known_transcript_pos: 0-based transcript start position for segment
    input param segment_length: integer representing segment length
    input param condition: character representing one of M, I or D.
    
    output: a dict with  updated genomic/transcript 0-based coords for the 
            entirety of the segment. Gaps will be represented as dashes.
    """
    genome_transcript_mapping = {}
    if condition == 'M':
        segment_final_genomic_pos = last_known_genomic_pos + segment_length
        transcript_pos = last_known_transcript_pos + 1
        for genomic_pos in range(last_known_genomic_pos + 1, segment_final_genomic_pos + 1):
           genome_transcript_mapping[transcript_pos] = genomic_pos
           transcript_pos = transcript_pos + 1
    elif condition == 'D':
        segment_final_genomic_pos = last_known_genomic_pos + segment_length - 1      
        for genomic_pos in range(last_known_genomic_pos, segment_final_genomic_pos + 1):
            genome_transcript_mapping[last_known_transcript_pos] = genomic_pos + 1
    elif condition == 'I':
        segment_final_genomic_pos = last_known_genomic_pos  
        segment_final_transcript_pos = last_known_transcript_pos + segment_length
        for transcript_pos in range(last_known_transcript_pos + 1, segment_final_transcript_pos + 1):
            genome_transcript_mapping[transcript_pos] = str(segment_final_genomic_pos) + 'I'
    return genome_transcript_mapping            


def clean_genomic_coord(genomic_coord):
    """
    cleans genomic coordinates when there's an insertion

    input param genomic_coord: 0-based genomic position, may contain an 'I'
    
    output: returns genomic_coord if integer, removes the 'I' and returns the integer part of
            the string otherwise
    """
    if isinstance(genomic_coord, str):
        genomic_coord = int(genomic_coord.replace('I', ''))
    return genomic_coord
        

def build_transcript_coordinate_maps(transcript_info):
    """
    creates a genome:transcript coordinate hashmap given a four-item list corresponding to each
    row in input file 1

    input param transcript_info: comprises of the transcript name, chromosome name, 0-based 
                                      chromosome start position and cigar string
    
    output: returns a dictionary of transcript to chromosome coordinate maps
    """
    cigar_tuples = split_cigar_into_tuples(transcript_info[3])
    last_known_transcript_pos = -1
    last_known_genomic_pos = transcript_info[2] - 1
    transcript_maps = {}
    for cigar_tuple in cigar_tuples:
        updated_coords = update_coords(last_known_genomic_pos, 
                                       last_known_transcript_pos,
                                       int(cigar_tuple[0]),
                                       cigar_tuple[1])
        if cigar_tuple[1] != 'D':
            transcript_maps.update(updated_coords)
        last_known_transcript_pos = list(updated_coords)[-1]
        last_known_genomic_pos = clean_genomic_coord(updated_coords[last_known_transcript_pos])
    return [transcript_info[1], transcript_maps]     

def build_hashmap(data_maps):
    """
    builds a hashable map for the transcripts 

    input param data_maps: pandas data frame; one transcript per row, corresponds to input_file_1

    output: a dictionary keyed by each transcript. contains chr and coordinate maps (another dict)
    """
    hashmaps = {}
    for index, transcript_row in data_maps.iterrows():
        transcript = transcript_row['transcript']
        chr = transcript_row['chr']
        chr_start_pos = transcript_row['chr_start_pos']
        cigar = transcript_row['cigar']
        hashmaps[transcript] = build_transcript_coordinate_maps([transcript, 
                                                                 chr, 
                                                                 chr_start_pos, 
                                                                 cigar])
    return hashmaps    

def query_transcript(data_query, hashmaps):
    """
    converts transcript coords to genomic coords

    input param data_query: transcript coords to be queried
    input param hashmaps: ground truth established using input_file_1

    output: returns a list of genomic map results by queried transcript
    """
    results = []
    for index, query_row in data_query.iterrows():
        transcript = query_row['transcript']
        transcript_coord = query_row['queried_coordinate']
        transcript_chr = maps[transcript][0]
        transcript_genomic_pos = hashmaps[transcript][1][transcript_coord]
        result = [transcript, transcript_coord, transcript_chr, transcript_genomic_pos]
        results.append(result)    
    return results

def write_output(result_list, output_file_name):
    """
    generalized output function 

    input param result_list: a list of lists to be written to file
    input param output_file_name: string denoting name of tab-delimited output file

    output: returns nothing, just writes data to file
    """
    with open(output_file_name, 'w+') as output_file:
        writer = csv.writer(output_file, delimiter='\t')
        for item in result_list:
            writer.writerow(item)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="translate transcript coords to genome coords")
    parser.add_argument('-i1', 
                        '--input-file-1', 
                        dest="input_file_1", 
                        action='store', 
                        default=None,
                        help='input file 1', 
                        required=True)
    parser.add_argument('-i2', 
                        '--input-file-2', 
                        dest="input_file_2", 
                        action='store', 
                        default=None,
                        help='input file 2', 
                        required=True)
    parser.add_argument('-o', 
                        '--output-file', 
                        dest="output_file", 
                        action='store', 
                        default=None,
                        help='output file', 
                        required=True)

    args = parser.parse_args()
    file_1 = args.input_file_1
    file_2 = args.input_file_2
    output_file = args.output_file

    data_maps = pandas.read_csv(file_1, sep="\t", header=None)
    data_maps.columns = ['transcript', 'chr', 'chr_start_pos', 'cigar']
    maps = build_hashmap(data_maps)
    
    data_query = pandas.read_csv(file_2, sep="\t", header=None)
    data_query.columns = ['transcript', 'queried_coordinate']
    transcript_coord_maps = query_transcript(data_query, maps)
    
    write_output(transcript_coord_maps, output_file)


"""
Quick tests
update_coords(2, -1, 8, 'M')
update_coords(10, 7, 7, 'D')
update_coords(17, 7, 6, 'M')
update_coords(23, 13, 2, 'I')    
update_coords(23, 15, 2, 'M')
update_coords(25, 17, 11, 'D')
update_coords(36, 17, 7, 'M')    
#transcript = "TR1"
#chromosome = "CHR1"
#genomic_start_pos = 3
#cigar_string = "8M7D6M2I2M11D7M"

re.findall('([0-9]+)([MID])', cigar_string)
ret = build_transcript_coordinate_maps(["TR1", "CHR1", 3, "8M7D6M2I2M11D7M"])
print(ret)
"""        
           
