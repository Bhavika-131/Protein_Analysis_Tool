import requests
import json
import time
from urllib.parse import urljoin


##So this script first asks you what you want to study, will fetch its proper name, accession ID.
##and get its details from uniprot, then Interpro. 
##I have also included its STRING interactions 
## and also uses its gene annotationa

##1) The script prompts the user for a protein name and then fetches comprehensive
#  information about that protein from multiple bioinformatics databases and services.

##2) retrieves data from UniProt (for basic protein information), InterPro (for protein family and domain information),
#  STRING (for protein-protein interactions and functional enrichment), and QuickGO (for Gene Ontology terms and annotations).

##3) All this diverse information is compiled into a single JSON file, providing a comprehensive overview
#  of the protein's sequence, structure, function, interactions, and biological roles, which can be used for
#  further analysis or interpretation from Gemini.

def fetch_protein_info(protein_name):
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"({protein_name})",
        "format": "json",
        "fields": "accession,protein_name",
        "size": 1
    }
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        data = response.json()
        if data['results']:
            result = data['results'][0]
            return result['primaryAccession'], result['proteinDescription']['recommendedName']['fullName']['value']
        else:
            return None, None
    except requests.exceptions.RequestException as e:
        print(f"Error occurred while fetching data: {e}")
        return None, None

def fetch_uniprot_info(accession, email):
    uniprot_base_url = "https://rest.uniprot.org/uniprotkb/"
    headers = {
        "Accept": "application/json",
        "User-Agent": f"Python script (mailto:{email})"
    }
    try:
        response = requests.get(f"{uniprot_base_url}{accession}", headers=headers)
        response.raise_for_status()
        uniprot_data = response.json()
        protein_info = {
            "accession": accession,
            "entry_type": uniprot_data.get('entryType'),
            "entry_name": uniprot_data.get('uniProtkbId'),
            "protein_name": uniprot_data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value'),
            "gene_name": next((gene.get('geneName', {}).get('value') for gene in uniprot_data.get('genes', [])), None),
            "organism": uniprot_data.get('organism', {}).get('scientificName'),
            "sequence": uniprot_data.get('sequence', {}).get('value'),
            "sequence_length": uniprot_data.get('sequence', {}).get('length'),
            "function": next((comment.get('texts', [{}])[0].get('value') for comment in uniprot_data.get('comments', []) if comment.get('commentType') == 'FUNCTION'), None),
            "subcellular_locations": [loc.get('location', {}).get('value') for comment in uniprot_data.get('comments', []) for loc in comment.get('subcellularLocations', []) if comment.get('commentType') == 'SUBCELLULAR LOCATION'],
            "ec_numbers": uniprot_data.get('proteinDescription', {}).get('ecNumbers', []),
            "keywords": [kw.get('name') for kw in uniprot_data.get('keywords', [])],
            "features": [{'type': f.get('type'), 'description': f.get('description')} for f in uniprot_data.get('features', [])]
        }
        return protein_info
    except requests.exceptions.RequestException as e:
        print(f"Error fetching UniProt data: {e}")
        return None

def fetch_comprehensive_interpro_info(accession, email):
    base_url = "https://www.ebi.ac.uk/interpro/api/"
    protein_url = urljoin(base_url, f"protein/UniProt/{accession}")
    headers = {
        "Accept": "application/json",
        "User-Agent": f"Python script (mailto:{email})"
    }
    try:
        response = requests.get(protein_url, headers=headers)
        response.raise_for_status()
        interpro_data = response.json()
        comprehensive_info = {
            "accession": accession,
            "metadata": interpro_data.get('metadata', {}),
            "entries": [],
            "structures": interpro_data.get('structure', []),
            "site_matches": interpro_data.get('site_matches', []),
            "taxonomy": interpro_data.get('taxonomy', {}),
            "proteomes": interpro_data.get('proteomes', []),
            "set_info": interpro_data.get('set_info', []),
            "extra_fields": interpro_data.get('extra_fields', {}),
            "member_databases": {}
        }
        for entry_type in ['entry_subset', 'unintegrated']:
            for entry in interpro_data.get(entry_type, []):
                entry_info = {
                    "accession": entry.get('accession'),
                    "name": entry.get('name'),
                    "type": entry.get('type'),
                    "source_database": entry.get('source_database'),
                    "integrated": entry_type == 'entry_subset',
                    "go_terms": entry.get('go_terms', []),
                    "locations": entry.get('locations', []),
                    "children": entry.get('children', []),
                    "counters": entry.get('counters', {}),
                    "signatures": entry.get('signatures', []),
                    "cross_references": entry.get('cross_references', [])
                }
                comprehensive_info["entries"].append(entry_info)
        return comprehensive_info
    except requests.exceptions.RequestException as e:
        print(f"Error fetching InterPro data: {e}")
        return None

def fetch_string_info(identifiers, species, email):
    base_url = "https://string-db.org/api"
    headers = {
        "Content-Type": "application/x-www-form-urlencoded",
        "Accept": "application/json",
        "User-Agent": f"Python script (mailto:{email})"
    }
    def api_call(method, params, output_format="json"):
        url = f"{base_url}/{output_format}/{method}"
        try:
            response = requests.post(url, data=params, headers=headers)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            print(f"Error calling STRING {method}: {str(e)}")
            return None
        finally:
            time.sleep(1)
    identifiers_str = "%0d".join([identifiers]) if isinstance(identifiers, str) else "%0d".join(identifiers)
    string_info = {
        "identifiers": api_call("get_string_ids", {
            "identifiers": identifiers_str,
            "species": species,
            "limit": 1,
            "echo_query": 1
        }),
        "network": api_call("network", {
            "identifiers": identifiers_str,
            "species": species,
            "required_score": 400
        }),
        "interaction_partners": api_call("interaction_partners", {
            "identifiers": identifiers_str,
            "species": species,
            "limit": 10
        }),
        "functional_annotation": api_call("functional_annotation", {
            "identifiers": identifiers_str,
            "species": species
        }),
        "enrichment": api_call("enrichment", {
            "identifiers": identifiers_str,
            "species": species
        }),
        "ppi_enrichment": api_call("ppi_enrichment", {
            "identifiers": identifiers_str,
            "species": species
        })
    }
    return string_info

def fetch_protein_go_terms(uniprot_id, email):
    base_url = "https://www.ebi.ac.uk/QuickGO/services"
    headers = {
        "Accept": "application/json",
        "User-Agent": f"Python script (mailto:{email})"
    }
    params = {
        "geneProductId": uniprot_id,
        "aspect": "biological_process",
        "limit": 10
    }
    try:
        response = requests.get(f"{base_url}/annotation/search", params=params, headers=headers)
        response.raise_for_status()
        data = response.json()
        go_terms = list(set([annotation['goId'] for annotation in data['results']]))
        return go_terms
    except requests.exceptions.RequestException as e:
        print(f"Error fetching GO terms for {uniprot_id}: {str(e)}")
        return []

def fetch_go_info(go_terms, email):
    base_url = "https://www.ebi.ac.uk/QuickGO/services"
    headers = {
        "Accept": "application/json",
        "User-Agent": f"Python script (mailto:{email})"
    }
    def api_call(endpoint, params=None):
        url = f"{base_url}{endpoint}"
        try:
            response = requests.get(url, params=params, headers=headers)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            print(f"Error calling QuickGO {endpoint}: {str(e)}")
            return None
        finally:
            time.sleep(1)
    go_terms_str = ",".join(go_terms) if isinstance(go_terms, list) else go_terms
    go_info = {
        "term_info": api_call(f"/ontology/go/terms/{go_terms_str}"),
        "ancestors": api_call(f"/ontology/go/terms/{go_terms_str}/ancestors"),
        "children": api_call(f"/ontology/go/terms/{go_terms_str}/children"),
        "descendants": api_call(f"/ontology/go/terms/{go_terms_str}/descendants"),
        "history": api_call(f"/ontology/go/terms/{go_terms_str}/history"),
        "xrefs": api_call(f"/ontology/go/terms/{go_terms_str}/xrefs"),
        "annotations": api_call("/annotation/search", {"goId": go_terms_str, "limit": 10})
    }
    return go_info

def main():
    email = "bhavikaberwal131@gmail.com"
    protein_name = input("What protein would you like to know about? ")
    
    print(f"\nFetching information for: {protein_name}")
    accession, full_name = fetch_protein_info(protein_name)
    
    if not accession:
        print(f"No results found for '{protein_name}'")
        return
    
    print(f"Protein: {full_name}")
    print(f"Accession: {accession}")
    
    all_data = {
        "uniprot": fetch_uniprot_info(accession, email),
        "interpro": fetch_comprehensive_interpro_info(accession, email),
        "string": fetch_string_info(accession, 9606, email),  # Assuming human (9606)
        "quickgo": {}
    }
    
    go_terms = fetch_protein_go_terms(accession, email)
    all_data["quickgo"]["go_terms"] = go_terms
    for go_term in go_terms:
        all_data["quickgo"][go_term] = fetch_go_info(go_term, email)
    
    # Save the data to a JSON file
    filename = f"{accession}_comprehensive_info.json"
    with open(filename, 'w') as f:
        json.dump(all_data, f, indent=2)
    
    print(f"\nComprehensive information has been saved to {filename}")

if __name__ == "__main__":
    main()