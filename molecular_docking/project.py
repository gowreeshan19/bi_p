import requests #Python to talk to the internet
import json #JSON is the language databases use to talk to each other

# 1. Define the PDB API Endpoint
url = "https://search.rcsb.org/rcsbsearch/v2/query"# PDB Search API endpointThis is the specific "address" of the RCSB PDB's search engine.

# 2. Define your parameters
EC_number = "3.4.21.4" #The ID for Trypsin
min_weight = 300 # Minimum molecular weight of ligand
max_weight = 800 # Maximum molecular weight of ligand

# 3. Build the query manually (Standard JSON format)
query_payload = { # The full query payload in JSON format
    "query": { # The main query structure
        "type": "group", # Grouping multiple conditions
        "logical_operator": "and", # All conditions must be met 
        "nodes": [ #
            {
                "type": "terminal",# Terminal node for EC number search
                "service": "text", # Service for text search
                "parameters": {
                    "attribute": "rcsb_polymer_entity.rcsb_ec_lineage.id",# Attribute to search for EC number
                    "operator": "exact_match",# Exact match operator
                    "value": EC_number# Value to match
                }
            },
            {
                "type": "terminal",# Terminal node for minimum molecular weight
                "service": "text",# Service for text search
                "parameters": {
                    "attribute": "chem_comp.formula_weight",# Attribute for molecular weight
                    "operator": "greater_or_equal",# Operator for minimum weight
                    "value": min_weight# Minimum weight value
                }
            },
            {
                "type": "terminal",# Terminal node for maximum molecular weight
                "service": "text",# Service for text search
                "parameters": {
                    "attribute": "chem_comp.formula_weight",# Attribute for molecular weight
                    "operator": "less_or_equal",# Operator for maximum weight
                    "value": max_weight# Maximum weight value
                }
            }
        ]
    },
    "return_type": "entry",# Return type of the search
    "request_options": {# Request options
        "return_all_hits": True# Return all matching hits
    }
}

# 4. Send the request
print("Searching PDB...")
try:
    response = requests.post(url, json=query_payload)# Send POST request with JSON payload
    
    if response.status_code == 200:# Check for successful response
        data = response.json()# Parse JSON response
        result_list = data.get('result_set', [])# Get the list of results
        
        # Get just the IDs
        ids = [entry['identifier'] for entry in result_list]# Extract identifiers from results
        
        print(f"Success! Found {len(ids)} structures.")# Print number of structures found
        print("First 10 results:", ids[:10])# Print first 10 IDs
    else:
        print("Error from PDB:", response.status_code)# Print error status code
        print(response.text)# Print error message

except Exception as e:# Handle connection errors
    print(f"Connection failed: {e}")# Print connection error message
