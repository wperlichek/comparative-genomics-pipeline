# How to Generate a Phylogenetic Tree Using the EBI Clustal Omega API

This guide explains how to use the EBI Clustal Omega API (Job Dispatcher) to generate a phylogenetic tree from protein sequences.

## 1. Fair-Use Policy
- Submit jobs in batches of **no more than 30 at a time**.
- Wait for results before submitting more jobs.
- Always provide a valid email address.

## 2. Submitting a Job (REST API)
You can submit a job to the Clustal Omega service using a POST request:

```
POST https://www.ebi.ac.uk/Tools/services/rest/clustalo/run/
```

**Parameters:**
- `sequence`: Your input protein sequences in FASTA format.
- `email`: Your email address (required).

**Example (Python):**
```python
import httpx

url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run/"
data = {
    "sequence": ">sp|P69905|HBA_HUMAN\nVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR\n",
    "email": "your@email.com"
}
response = httpx.post(url, data=data)
job_id = response.text.strip()
```

**Note:**
Clustal Omega alignments may not be fully deterministic between runs, especially if the input sequence order changes or the service is updated. For reproducible downstream results (like conservation scores), always use the same input order and avoid regenerating alignments unless your sequences change.

## 3. Checking Job Status
Check the status of your job using:
```
GET https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}
```

**Example:**
```python
status = httpx.get(f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}").text.strip()
```

## 4. Retrieving the Phylogenetic Tree
Once the job status is `FINISHED`, retrieve the tree:
```
GET https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/phyloxml
```
- `phyloxml` returns the tree in PhyloXML format.
- You can also use `ph` for Newick format.

**Example:**
```python
tree = httpx.get(f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/phyloxml").text
```

## 5. Example Workflow
1. Submit your protein sequences in FASTA format.
2. Poll the status endpoint until the job is `FINISHED`.
3. Download the phylogenetic tree result.

## 6. References
- [EBI Job Dispatcher Web Services Docs](https://www.ebi.ac.uk/jdispatcher/docs//webservices/)
- [Clustal Omega API](https://www.ebi.ac.uk/Tools/services/rest/clustalo)

---
For more advanced usage, see the [official documentation](https://www.ebi.ac.uk/jdispatcher/docs//webservices/).
