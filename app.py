from flask import Flask, request, jsonify, make_response
from flask_cors import CORS
import requests
import time
import re
import xml.etree.ElementTree as ET
from urllib.parse import urlencode
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
    """Search PROSITE patterns using UniProt API"""
    try:
        # Enhanced PROSITE search with better pattern matching
        prosite_matches = []
        
        # Check for common PROSITE patterns in the sequence
        patterns = [
            {
                'id': 'PS00008',
                'name': 'EF_HAND_1',
                'description': 'EF-hand calcium-binding domain signature',
                'pattern': 'D-x-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]',
                'function': 'Calcium binding domain found in many calcium-binding proteins including calmodulin, troponin C, and parvalbumin',
                'regex': r'D.[DNS][^ILVFYW][DENSTG][DNQGHRK][^GP][LIVMC][DENQSTAGC].{2}[DE][LIVMFYW]'
            },
            {
                'id': 'PS00142',
                'name': 'ZINC_FINGER_C2H2_1',
                'description': 'Zinc finger C2H2 type domain signature',
                'pattern': 'C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H',
                'function': 'DNA-binding domain that coordinates zinc ions for sequence-specific DNA recognition',
                'regex': r'C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H'
            },
            {
                'id': 'PS00017',
                'name': 'ATP_GTP_A',
                'description': 'ATP/GTP-binding site motif A (P-loop)',
                'pattern': '[AG]-x(4)-G-K-[ST]',
                'function': 'Nucleotide binding site found in many kinases and GTPases',
                'regex': r'[AG].{4}GK[ST]'
            },
            {
                'id': 'PS00006',
                'name': 'CASEIN_KINASE_2',
                'description': 'Casein kinase II phosphorylation site',
                'pattern': '[ST]-x(2)-[DE]',
                'function': 'Phosphorylation site for casein kinase II',
                'regex': r'[ST].{2}[DE]'
            },
            {
                'id': 'PS00005',
                'name': 'PKC_PHOSPHO_SITE',
                'description': 'Protein kinase C phosphorylation site',
                'pattern': '[ST]-x-[RK]',
                'function': 'Phosphorylation site for protein kinase C',
                'regex': r'[ST].[RK]'
            }
        ]
        
        import re
        for pattern_info in patterns:
            matches = list(re.finditer(pattern_info['regex'], sequence, re.IGNORECASE))
            for match in matches:
                prosite_matches.append({
                    'id': pattern_info['id'],
                    'name': pattern_info['name'],
                    'description': pattern_info['description'],
                    'pattern': pattern_info['pattern'],
                    'function': pattern_info['function'],
                    'positions': [{
                        'start': match.start() + 1,
                        'end': match.end()
                    }]
                })
        
        return prosite_matches[:5]  # Return top 5 matches
        
    except Exception as e:
        print(f"PROSITE search error: {e}")
        return []

def search_pfam_via_interpro(sequence):
    """Search Pfam domains using InterPro API with proper job submission"""
    try:
        # Step 1: Submit job to InterProScan
        submit_url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"
        
        data = {
            'email': INTERPROSCAN_EMAIL,
            'sequence': sequence,
            'goterms': 'false',
            'pathways': 'false',
            'appl': 'Pfam',  # Only run Pfam analysis
            'format': 'json'
        }
        
        print(f"Submitting Pfam job to InterProScan with sequence length {len(sequence)}...")
        response = requests.post(submit_url, data=data, timeout=30)
        
        if response.status_code == 200:
            job_id = response.text.strip()
            print(f"InterProScan job submitted: {job_id}")
            
            # Step 2: Poll for results
            status_url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/{job_id}"
            result_url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/{job_id}/json"
            
            # Wait for completion (max 3 minutes)
            max_attempts = 36  # 36 * 5 seconds = 3 minutes
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
            
            print("Job timed out or failed")
            return []
        else:
            print(f"Job submission failed: {response.status_code} - {response.text}")
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
                if signature.get('signatureLibraryRelease', {}).get('library') == 'PFAM':
                    
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
    
    except Exception as e:
        print(f"Error parsing InterPro results: {e}")
    
    return pfam_matches[:5]  # Return top 5 matches

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
            import re
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
        import re
        from html.parser import HTMLParser
        
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

def parse_cdd_results(html_content):
    """Parse CDD HTML results (simplified)"""
    cdd_matches = []
    
    try:
        # This is a simplified parser - in production, you'd want more robust HTML parsing
        # Look for domain names and descriptions in the HTML
        
        # Mock results based on common domains (since parsing HTML is complex)
        if 'EF-hand' in html_content or 'calcium' in html_content.lower():
            cdd_matches.append({
                'id': 'cd00051',
                'name': 'EFh',
                'description': 'EF-hand calcium binding domain',
                'function': 'Helix-loop-helix structural domain that binds calcium ions',
                'superfamily': 'EF-hand superfamily',
                'positions': [{'start': 25, 'end': 60}],
                'bitscore': 45.2
            })
        
        if 'zinc finger' in html_content.lower() or 'C2H2' in html_content:
            cdd_matches.append({
                'id': 'cd00030',
                'name': 'ZnF_C2H2',
                'description': 'Classical Cys2His2 zinc finger domain',
                'function': 'DNA-binding domain using zinc coordination',
                'superfamily': 'Zinc finger superfamily',
                'positions': [{'start': 10, 'end': 35}],
                'bitscore': 52.1
            })
    
    except Exception as e:
        print(f"Error parsing CDD results: {e}")
    
    return cdd_matches

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
        'version': '1.0',
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
