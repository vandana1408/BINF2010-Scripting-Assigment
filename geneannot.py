#!/usr/bin/python3

"""
BINF2010 ASSIGNMENT 
"""

import sys 
import re 
import subprocess 
import urllib.request 
import requests 
import urllib.parse
import json
import csv 
import pandas as pd 

def check_file(filename):
    """
    Makes sure that the file is in fasta format using regex. 

    Arguments: 
        filename (str) - the input filename that the user wants to check. 

    Errors: 
        TypeError      - Occurs when the filename does not end in 'fasta'
    
    Returns: 
        Returns the filename 
    """
    regex = r'^.*\.fasta$'
    if not re.fullmatch(regex, filename): 
        raise TypeError("File must be in fasta format!")
    
    return filename

def get_orfs(filename): 
    """
    Searches the file for orfs and returns the three longest orfs that are greater than 150bp.

    Arguments: 
        filename (str) - the input filename that the user wants to check. 
    
    Returns: 
        Returns a list of dictionaries containing the three longest orfs, each dictionary has keys "orf" and "length". 
        The dictionaries are sorted by the length key. 
    
    """
    subprocess.run(["getorf", filename], stderr=subprocess.DEVNULL, text=True, input="getorf_output.orf")

    ORFFILE = open("getorf_output.orf", "r")

    orfs = ORFFILE.read()

    orfs_list = orfs.split(">")

    orfs_list.remove(orfs_list[0])

    isolated_orfs_list = []

    for orf in orfs_list: 
        str_orf = str(orf)
        orf_range = re.match(r'.*\[([0-9]+) - ([0-9]+)\].*', str_orf)
        length_of_orfs = int(orf_range.group(2)) - int(orf_range.group(1))
        if length_of_orfs >= 150 or length_of_orfs <= -150: 
            orf_dict = dict(orf = str_orf, length = abs(length_of_orfs))
            isolated_orfs_list.append(orf_dict)

    longest_orfs = sorted(isolated_orfs_list, key=lambda d: d['length'], reverse=True)[:3]

    return longest_orfs


def search_pdb_url(sequence): 
    """
    Creates a url to post to a server and retrieve the highest hit. 

    Arguments: 
        sequence (str) - the input sequence that becomes a query for the pdb database. 
    
    Returns: 
        Returns a url. 
    
    """
    base_url = "https://search.rcsb.org/rcsbsearch/v1/query?json="

    search_request = { 
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 1,
                "identity_cutoff": 0.4,
                "target": "pdb_protein_sequence",
                "value": sequence
            }
        }, 
        "request_options": {
            "scoring_strategy": "sequence"
        },
        "return_type": "polymer_entity" 
    }

    real_search_req = urllib.parse.quote_plus(json.dumps(search_request))

    return base_url + real_search_req 

def get_highest_hits(longest_orfs): 
    """
    Peforms a search on the longest orfs against the pdb database to find the highest hit, stores values necessary for the csv table into separate lists, 
    where the indices represent the order. 

    Arguments: 
        longest_orfs (list) - list of orfs to be searched in the pdb database
    
    Returns: 
        Returns 5 separate lists containing the results from the search and other necessary information for the csv table. 
    
    """
    pdbs = []
    e_values = []
    strands = []
    starts = []
    ends = []

    for orfs in longest_orfs:
        remove_new_lines = str(orfs['orf'].replace("\n", ""))
        orf_raw_sequence = re.match(r'.*\[([0-9]+) - ([0-9]+)\].*sequence([A-Z]+).*', remove_new_lines)
        start = orf_raw_sequence.group(1)
        end = orf_raw_sequence.group(2)
        strand = 'FORWARD' if int(end) > int(start) else 'REVERSE'

        url_for_search = search_pdb_url(orf_raw_sequence.group(3))
        response = requests.get(url_for_search)

        try: 
            results = json.loads(response.text)
            pdb_id = results['result_set'][0]['identifier']
            e_value = results['result_set'][0]['services'][0]['nodes'][0]['match_context'][0]['evalue']
        except: 
            pdb_id = '-'
            e_value = '-'

        pdbs.append(pdb_id)
        e_values.append(e_value)
        strands.append(strand)
        starts.append(start)
        ends.append(end)

    return pdbs, e_values, strands, starts, ends

def organise_results_into_csv(pdbs, e_values, strands, starts, ends): 
    """
    Creates a csv table and sorts the rows based on the start position of each orf. 

    Arguments: 
        pdbs (list) - list of pdb_ids 
        e_values (list) - list of e_values
        strands (list) - list of strand orientations (FORWARD or REVERSE)
        starts (list) - list of start positions 
        ends (list) - list of end positions 

    Returns: 
        Returns a csv table that contains data sorted by the start positions of the longest orfs. 
    
    """
    csv_header = ['Start', 'End', 'Strand', 'PDB_ID', 'E_value']

    with open('orf_results.csv', 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(csv_header)
        for i in range(len(pdbs)):
            data = [starts[i], ends[i], strands[i], pdbs[i], e_values[i]]
            writer.writerow(data)

        f.close()

    orf_results = pd.read_csv("orf_results.csv")
    sorted_df = orf_results.sort_values(by=["Start"], ascending=True)
    final_df = sorted_df.to_csv('orf_results.csv', index=False)
    return final_df

if __name__ == "__main__": 
    #filename = sys.argv[1]
    sys.tracebacklimit = 0

    filename = check_file(sys.argv[1])

    longest_orfs = get_orfs(filename)

    pdbs, e_values, strands, starts, ends = get_highest_hits(longest_orfs)

    organise_results_into_csv(pdbs, e_values, strands, starts, ends)

    subprocess.run(["cat", "orf_results.csv"], stderr=subprocess.DEVNULL, text=True)
    subprocess.run(["rm", "orf_results.csv"], stderr=subprocess.DEVNULL, text=True)
    subprocess.run(["rm", "getorf_output.orf"], stderr=subprocess.DEVNULL, text=True)
    

