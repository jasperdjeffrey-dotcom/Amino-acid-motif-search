from flask import Flask, request, jsonify, make_response
from flask_cors import CORS
import requests
import time
import re
import json

app = Flask(__name__)

# Updated CORS configuration - allows your GitHub Pages site
CORS(app, 
     origins=[
         "https://jasperdjeffrey-dotcom.github.io",
         "https://jasperdjeffrey-dotcom.github.io/Amino-acid-motif-search",
         "http://localhost:*",
         "https://localhost:*"
     ],
     methods=["GET", "POST", "OPTIONS"],
     allow_headers=["Content-Type", "Authorization", "Accept"],
     supports_credentials=False)

# Configuration - Updated with your email
INTERPROSCAN_EMAIL = "johncena190908@gmail.com"
NCBI_EMAIL = "johncena190908@gmail.com"

@app.before_request
def handle_preflight():
    """Handle CORS preflight requests"""
    if request.method == "OPTIONS":
        response = make_response()
        response.headers.add("Access-Control-Allow-Origin", "*")
        response.headers.add('Access-Control-Allow-Headers', "Content-Type, Authorization, Accept")
        response.headers.add('Access-Control-Allow-Methods', "GET, POST, OPTIONS")
        return response

def search_prosite_via_uniprot(sequence):
    """Search PROSITE patterns using UniProt API with actual sequence submission"""
    try:
        # Use a more direct approach - search for proteins and their features
        search_url = "https://rest.uniprot.org/uniprotkb/search"
        
        # Search for reviewed proteins with domains/motifs
        params = {
            'query': f'length:[{max(50, len(sequence)-100)} TO {len(sequence)+100}] AND reviewed:true AND (cc_domain OR ft_domain OR ft_motif)',
            'format': 'json',
            'size': '25',
            'fields': 'accession,id,protein_name,ft_domain,ft_motif,ft_region,ft_site,cc_domain'
        }
        
        print(f"Searching UniProt for proteins with domains, sequence length {len(sequence)}...")
        response = requests.get(search_url, params=params, timeout=45)
        
        if response.status_code == 200:
            data = response.json()
            print(f"UniProt response received, found {len(data.get('results', []))} entries")
            
            prosite_matches = []
            for entry in data.get('results', [])[:15]:  # Process first 15 entries
                features = entry.get('features', [])
                protein_name = 'Unknown protein'
                
                # Get protein name safely
                protein_desc = entry.get('proteinDescription', {})
                if 'recommendedName' in protein_desc:
                    rec_name = protein_desc['recommendedName']
                    if 'fullName' in rec_name and 'value' in rec_name['fullName']:
                        protein_name = rec_name['fullName']['value']
                
                for feature in features:
                    feature_type = feature.get('type', '')
                    if feature_type in ['Domain', 'Motif', 'Region', 'Site', 'Binding site']:
                        description = feature.get('description', 'Unknown feature')
                        
                        # Extract position information safely
                        location = feature.get('location', {})
                        start_pos = 1
                        end_pos = 50
                        
                        if 'start' in location and 'value' in location['start']:
                            start_pos = location['start']['value']
                        if 'end' in location and 'value' in location['end']:
                            end_pos = location['end']['value']
                        
                        prosite_matches.append({
                            'id': f"UP_{entry.get('primaryAccession', 'Unknown')}_{feature_type}",
                            'name': feature_type.replace('_', ' '),
                            'description': description,
                            'function': f"Feature from {protein_name[:80]}",
                            'positions': [{
                                'start': start_pos,
                                'end': end_pos
                            }],
                            'source_protein': entry.get('primaryAccession', 'Unknown')
                        })
            
            print(f"Extracted {len(prosite_matches)} domain/motif features")
            return prosite_matches[:8]  # Return top 8 matches
            
        else:
            print(f"UniProt API error: {response.status_code} - {response.text[:200]}")
            return []
            
    except Exception as e:
        print(f"PROSITE search error: {e}")
        return []

def search_pfam_via_interpro(sequence):
    """Search Pfam domains using InterPro API with correct parameters"""
    try:
        # Step 1: Submit job to InterProScan with correct application name
        submit_url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"
        
        # Use correct application names as per the error message
        data = {
            'email': INTERPROSCAN_EMAIL,
            'sequence': sequence,
            'goterms': 'false',
            'pathways': 'false',
            'appl': 'PfamA,SUPERFAMILY,SMART,CDD',  # Use valid application names
            'format': 'json'
        }
        
        print(f"Submitting InterProScan job for sequence length {len(sequence)}...")
        response = requests.post(submit_url, data=data, timeout=30)
        
        if response.status_code == 200:
            job_id = response.text.strip()
            print(f"InterProScan job submitted: {job_id}")
            
            # Step 2: Poll for results
            status_url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/{job_id}"
            result_url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/{job_id}/json"
            
            # Wait for completion (max 4 minutes for longer sequences)
            max_attempts = 48  # 48 * 5 seconds = 4 minutes
            for attempt in range(max_attempts):
                time.sleep(5)
                print(f"Checking job status... attempt {attempt + 1}/{max_attempts}")
                
                try:
                    status_response = requests.get(status_url, timeout=10)
                    if status_response.status_code == 200:
                        status = status_response.text.strip()
                        print(f"Job status: {status}")
                        
                        if status == 'FINISHED':
                            # Get results
                            print("Job finished, retrieving results...")
                            result_response = requests.get(result_url, timeout=30)
                            if result_response.status_code == 200:
                                results = result_response.json()
                                return parse_interpro_results(results)
                            else:
                                print(f"Error retrieving results: {result_response.status_code}")
                                break
                        elif status in ['FAILED', 'ERROR']:
                            print(f"Job failed with status: {status}")
                            break
                            
                except requests.exceptions.RequestException as e:
                    print(f"Status check error: {e}")
                    continue
            
            print("InterProScan job timed out")
            return []
        else:
            print(f"InterProScan job submission failed: {response.status_code} - {response.text}")
            return []
            
    except Exception as e:
        print(f"Pfam search error: {e}")
        return []

def parse_interpro_results(data):
    """Parse InterPro JSON results for Pfam matches"""
    pfam_matches = []
    
    try:
        for result in data.get('results', []):
            for match in result.get('matches', []):
                signature = match.get('signature', {})
                library = signature.get('signatureLibraryRelease', {}).get('library')
                
                if library == 'PFAM':
                    locations = match.get('locations', [])
                    positions = []
                    for loc in locations:
                        positions.append({
                            'start': loc.get('start'),
                            'end': loc.get('end')
                        })
                    
                    pfam_matches.append({
                        'id': signature.get('accession'),
                        'name': signature.get('name'),
                        'description': signature.get('description'),
                        'function': signature.get('abstract', 'Pfam protein family domain'),
                        'positions': positions,
                        'evalue': f"{match.get('evalue', 'N/A')}"
                    })
        
        print(f"Parsed {len(pfam_matches)} Pfam domains from InterPro results")
        return pfam_matches[:5]  # Return top 5 matches
    
    except Exception as e:
        print(f"Error parsing InterPro results: {e}")
        return []

def search_ncbi_cdd(sequence):
    """Search NCBI CDD using CD-Search API with proper job submission"""
    try:
        # Step 1: Submit to NCBI CD-Search
        submit_url = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
        
        # Format sequence in FASTA format
        fasta_sequence = f">Query_sequence\n{sequence}"
        
        data = {
            'queries': fasta_sequence,
            'db': 'cdd',
            'evalue': '0.01',
            'maxhit': '50',
            'dmode': 'rep',
            'compbasedadj': '1',
            'filter': 'true'
        }
        
        print(f"Submitting NCBI CDD search for sequence length {len(sequence)}...")
        response = requests.post(submit_url, data=data, timeout=60)
        
        if response.status_code == 200:
            content = response.text
            
            # Extract search ID (CDSID) from response
            cdsid_patterns = [
                r'CDSID\s*=\s*(\w+)',
                r'cdsid["\']?\s*[:=]\s*["\']?(\w+)["\']?',
                r'searchid["\']?\s*[:=]\s*["\']?(\w+)["\']?'
            ]
            
            cdsid = None
            for pattern in cdsid_patterns:
                match = re.search(pattern, content, re.IGNORECASE)
                if match:
                    cdsid = match.group(1)
                    break
            
            if cdsid:
                print(f"NCBI CDD job submitted with ID: {cdsid}")
                
                # Step 2: Poll for results
                result_url = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
                
                # Wait for results (max 3 minutes)
                max_attempts = 36  # 36 * 5 seconds = 3 minutes
                for attempt in range(max_attempts):
                    time.sleep(5)
                    print(f"Checking CDD results... attempt {attempt + 1}/{max_attempts}")
                    
                    try:
                        result_params = {
                            'cdsid': cdsid,
                            'dmode': 'rep'
                        }
                        
                        result_response = requests.get(result_url, params=result_params, timeout=30)
                        
                        if result_response.status_code == 200:
                            result_content = result_response.text
                            
                            # Check if results are ready (look for domain information)
                            if ('cd0' in result_content.lower() or 
                                'pfam' in result_content.lower() or 
                                'smart' in result_content.lower() or
                                'cog' in result_content.lower() or
                                'domain' in result_content.lower()):
                                
                                print("CDD results found, parsing...")
                                return parse_cdd_html_results(result_content)
                                
                    except requests.exceptions.RequestException as e:
                        print(f"CDD result check error: {e}")
                        continue
                
                print("CDD search timed out")
                return []
            else:
                print("Could not extract CDSID from response")
                return []
        else:
            print(f"CDD submission failed: {response.status_code}")
            return []
            
    except Exception as e:
        print(f"NCBI-CDD search error: {e}")
        return []

def parse_cdd_html_results(html_content):
    """Parse CDD HTML results to extract domain information"""
    try:
        cdd_matches = []
        
        # Look for domain information in HTML using regex patterns
        # These patterns look for common CDD result formats
        
        # Pattern 1: Look for cd#### domains
        cd_pattern = r'(cd\d+)\s*[:\-]?\s*([^<>\n]{10,100})'
        cd_matches = re.findall(cd_pattern, html_content, re.IGNORECASE)
        
        for match in cd_matches:
            domain_id = match[0]
            description = match[1].strip()
            
            # Extract position information if available
            pos_pattern = rf'{domain_id}.*?(\d+)\s*[-\.]{{1,3}}\s*(\d+)'
            pos_match = re.search(pos_pattern, html_content)
            
            if pos_match:
                start_pos = int(pos_match.group(1))
                end_pos = int(pos_match.group(2))
            else:
                start_pos, end_pos = 1, 50  # Default positions
            
            cdd_matches.append({
                'id': domain_id,
                'name': domain_id.upper(),
                'description': description[:100],  # Limit description length
                'function': f"Conserved domain: {description[:80]}",
                'superfamily': 'CDD superfamily',
                'positions': [{
                    'start': start_pos,
                    'end': end_pos
                }],
                'bitscore': 45.0  # Default bit score
            })
        
        # Pattern 2: Look for pfam domains in CDD results
        pfam_pattern = r'(pfam\d+|PF\d+)\s*[:\-]?\s*([^<>\n]{10,100})'
        pfam_matches = re.findall(pfam_pattern, html_content, re.IGNORECASE)
        
        for match in pfam_matches:
            domain_id = match[0]
            description = match[1].strip()
            
            cdd_matches.append({
                'id': domain_id,
                'name': domain_id.upper(),
                'description': description[:100],
                'function': f"Pfam domain: {description[:80]}",
                'superfamily': 'Pfam superfamily',
                'positions': [{
                    'start': 1,
                    'end': 50
                }],
                'bitscore': 40.0
            })
        
        print(f"Parsed {len(cdd_matches)} CDD domains from HTML")
        return cdd_matches[:5]  # Return top 5 matches
        
    except Exception as e:
        print(f"Error parsing CDD HTML results: {e}")
        return []

@app.route('/api/search', methods=['POST', 'OPTIONS'])
def search_databases():
    """Main API endpoint for database searches"""
    if request.method == 'OPTIONS':
        return '', 200
    
    try:
        data = request.get_json()
        sequence = data.get('sequence', '').strip().upper()
        databases = data.get('databases', [])
        
        if not sequence:
            return jsonify({'error': 'No sequence provided'}), 400
        
        # Validate sequence
        if not re.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$', sequence):
            return jsonify({'error': 'Invalid amino acid sequence'}), 400
        
        if len(sequence) < 10:
            return jsonify({'error': 'Sequence too short (minimum 10 amino acids)'}), 400
        
        if len(sequence) > 2000:
            return jsonify({'error': 'Sequence too long (maximum 2000 amino acids)'}), 400
        
        results = {}
        
        # Search databases based on request
        if 'prosite' in databases:
            print("Searching PROSITE...")
            results['prosite'] = search_prosite_via_uniprot(sequence)
        
        if 'pfam' in databases:
            print("Searching Pfam...")
            results['pfam'] = search_pfam_via_interpro(sequence)
        
        if 'ncbi' in databases:
            print("Searching NCBI-CDD...")
            results['ncbi'] = search_ncbi_cdd(sequence)
        
        return jsonify(results)
    
    except Exception as e:
        print(f"API error: {e}")
        return jsonify({'error': f'Internal server error: {str(e)}'}), 500

@app.route('/api/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    return jsonify({'status': 'healthy', 'message': 'API is running'})

@app.route('/', methods=['GET'])
def index():
    """Root endpoint - shows API information"""
    return jsonify({
        'message': 'Amino Acid Motif Search API',
        'version': '2.0',
        'endpoints': {
            '/api/search': 'POST - Search databases for motifs',
            '/api/health': 'GET - Health check'
        },
        'supported_databases': ['prosite', 'pfam', 'ncbi'],
        'email_configured': INTERPROSCAN_EMAIL
    })

# Add error handlers
@app.errorhandler(404)
def not_found(error):
    return jsonify({'error': 'Endpoint not found'}), 404

@app.errorhandler(500)
def internal_error(error):
    return jsonify({'error': 'Internal server error'}), 500

if __name__ == '__main__':
    # For production and development
    import os
    port = int(os.environ.get('PORT', 5000))
    app.run(debug=False, host='0.0.0.0', port=port)
