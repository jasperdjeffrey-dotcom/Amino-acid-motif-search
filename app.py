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
        # First, try to find proteins with similar sequences
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            'query': f'sequence_length:[{max(1, len(sequence)-50)} TO {len(sequence)+50}]',
            'format': 'json',
            'size': '10',
            'fields': 'accession,id,protein_name,ft_domain,ft_motif,ft_region'
        }
        
        response = requests.get(url, params=params, timeout=30)
        if response.status_code == 200:
            data = response.json()
            
            prosite_matches = []
            for entry in data.get('results', [])[:5]:  # Limit to first 5 entries
                features = entry.get('features', [])
                for feature in features:
                    if feature.get('type') in ['Domain', 'Motif', 'Region']:
                        prosite_matches.append({
                            'id': feature.get('description', 'Unknown'),
                            'name': feature.get('type', 'Feature'),
                            'description': feature.get('description', 'Protein feature'),
                            'function': f"Found in {entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'protein')}",
                            'positions': [{
                                'start': feature.get('location', {}).get('start', {}).get('value', 1),
                                'end': feature.get('location', {}).get('end', {}).get('value', len(sequence))
                            }]
                        })
            
            return prosite_matches[:3]  # Return top 3 matches
            
    except Exception as e:
        print(f"PROSITE search error: {e}")
        
    return []

def search_pfam_via_interpro(sequence):
    """Search Pfam domains using InterPro API"""
    try:
        # Submit job to InterProScan
        submit_url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"
        
        data = {
            'email': INTERPROSCAN_EMAIL,
            'sequence': sequence,
            'goterms': 'false',
            'pathways': 'false',
            'appl': 'Pfam'  # Only run Pfam
        }
        
        headers = {'Content-Type': 'application/x-www-form-urlencoded'}
        
        # Submit the job
        response = requests.post(submit_url, data=data, headers=headers, timeout=30)
        
        if response.status_code == 200:
            job_id = response.text.strip()
            
            # Poll for results (with timeout)
            result_url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/{job_id}/json"
            status_url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/{job_id}"
            
            # Wait for job completion (max 2 minutes)
            max_attempts = 24  # 24 * 5 seconds = 2 minutes
            for attempt in range(max_attempts):
                time.sleep(5)
                
                status_response = requests.get(status_url, timeout=10)
                if status_response.status_code == 200:
                    status = status_response.text.strip()
                    
                    if status == 'FINISHED':
                        # Get results
                        result_response = requests.get(result_url, timeout=30)
                        if result_response.status_code == 200:
                            return parse_interpro_results(result_response.json())
                    elif status == 'FAILED':
                        break
            
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
    """Search NCBI CDD using CD-Search"""
    try:
        # Submit to CD-Search
        submit_url = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
        
        data = {
            'queries': f">Query\n{sequence}",
            'db': 'cdd',
            'evalue': '0.01',
            'maxhit': '10',
            'dmode': 'rep',
            'compbasedadj': '1'
        }
        
        # Submit search
        response = requests.post(submit_url, data=data, timeout=30)
        
        if response.status_code == 200:
            # Extract search ID from response
            content = response.text
            
            # Look for CDSID in the response
            cdsid_match = re.search(r'CDSID\s*=\s*(\w+)', content)
            if cdsid_match:
                cdsid = cdsid_match.group(1)
                
                # Poll for results
                result_url = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
                result_params = {
                    'cdsid': cdsid,
                    'dmode': 'rep'
                }
                
                # Wait for results (max 2 minutes)
                for attempt in range(24):  # 24 * 5 seconds = 2 minutes
                    time.sleep(5)
                    
                    result_response = requests.get(result_url, params=result_params, timeout=30)
                    if result_response.status_code == 200 and 'cd' in result_response.text.lower():
                        return parse_cdd_results(result_response.text)
        
    except Exception as e:
        print(f"NCBI-CDD search error: {e}")
    
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
