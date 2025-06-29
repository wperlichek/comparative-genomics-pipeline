"""Test data fixtures for comparative genomics pipeline tests."""

# Sample UniProt FASTA responses
SAMPLE_UNIPROT_FASTA = {
    "P35498": """>sp|P35498|SCN1A_HUMAN Sodium channel protein type 1 subunit alpha OS=Homo sapiens OX=9606 GN=SCN1A PE=1 SV=1
MAASDSEYRTRSEAETLSITDMEAGTDVQKADGDFVQGQHQEVSKVQGTGTDSGAFQHGPQATP
TPQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQR""",
    
    "A2APX8": """>sp|A2APX8|SCN1A_MOUSE Sodium channel protein type 1 subunit alpha OS=Mus musculus OX=10090 GN=Scn1a PE=1 SV=1
MAASDSEYRTRSEAETLSITDMEAGTDVQKADGDFVQGQHQEVSKVQGTGTDSGAFQHGPQATP
TPQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQR"""
}

# Sample UniProt variant data
SAMPLE_UNIPROT_VARIANTS = {
    "P35498": {
        "features": [
            {
                "type": "Natural variant",
                "location": {"start": 177, "end": 177},
                "wildType": "R",
                "alternativeSequence": "W", 
                "description": "in GEFS+2; dbSNP:rs121909902"
            },
            {
                "type": "Natural variant", 
                "location": {"start": 1648, "end": 1648},
                "wildType": "R",
                "alternativeSequence": "H",
                "description": "in SMEI; dbSNP:rs121909901"
            },
            {
                "type": "Transmembrane region",
                "location": {"start": 121, "end": 143},
                "description": "Helical"
            }
        ]
    }
}

# Sample EBI Clustal Omega responses
SAMPLE_EBI_ALIGNMENT = """>sp|P35498|SCN1A_HUMAN
MAASDSEYRTRSEAETLSITDMEAGTDVQKADGDFVQGQHQEVSKVQGTGTDSGAFQHGPQATP
TPQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQR
>sp|A2APX8|SCN1A_MOUSE
MAASDSEYRTRSEAETLSITDMEAGTDVQKADGDFVQGQHQEVSKVQGTGTDSGAFQHGPQATP
TPQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQRQR"""

SAMPLE_EBI_TREE = "(P35498:0.00000,A2APX8:0.00000);"

# Sample gene configuration
SAMPLE_GENE_CONFIG = {
    "SCN1A": [
        {
            "species": "Homo sapiens",
            "uniprot_id": "P35498",
            "entrez_protein_id": "ENSP05155037736.1"
        },
        {
            "species": "Mus musculus",
            "uniprot_id": "A2APX8", 
            "entrez_protein_id": "ENSMUSP00215039336.1"
        },
        {
            "species": "Macaca mulatta",
            "uniprot_id": "A0A1D5QX02",
            "entrez_protein_id": "XP_001101023.1"
        }
    ],
    "KCNQ2": [
        {
            "species": "Homo sapiens",
            "uniprot_id": "O43526",
            "entrez_protein_id": "ENSP00000355839.4"
        },
        {
            "species": "Mus musculus", 
            "uniprot_id": "Q9Z0N7",
            "entrez_protein_id": "ENSMUSP00000063773.11"
        }
    ]
}

# Sample conservation data
SAMPLE_CONSERVATION_DATA = [
    {"position": 1, "entropy_gaps": 0.0, "entropy_nogaps": 0.0, "aa_counts_gaps": "M:3", "aa_counts_nogaps": "M:3"},
    {"position": 2, "entropy_gaps": 0.0, "entropy_nogaps": 0.0, "aa_counts_gaps": "A:3", "aa_counts_nogaps": "A:3"},
    {"position": 3, "entropy_gaps": 0.636, "entropy_nogaps": 0.636, "aa_counts_gaps": "A:2,T:1", "aa_counts_nogaps": "A:2,T:1"},
    {"position": 4, "entropy_gaps": 0.0, "entropy_nogaps": 0.0, "aa_counts_gaps": "S:3", "aa_counts_nogaps": "S:3"},
    {"position": 5, "entropy_gaps": 0.0, "entropy_nogaps": 0.0, "aa_counts_gaps": "D:3", "aa_counts_nogaps": "D:3"}
]

# Sample alignment with gaps
SAMPLE_ALIGNMENT_WITH_GAPS = """>seq1|HUMAN
MA-SDSEYR
>seq2|MOUSE
MAASDSEYR
>seq3|CHICKEN
MA-ADSEYR"""

# Sample PDB structure data
SAMPLE_PDB_DATA = {
    "P35498": {
        "structures": [
            {
                "pdb_id": "6AGF",
                "method": "Cryo-EM",
                "resolution": "3.20",
                "chains": ["A"],
                "title": "Structure of human Nav1.1 voltage-gated sodium channel"
            }
        ]
    }
}

# Mock HTTP responses for testing
MOCK_HTTP_RESPONSES = {
    "uniprot_fasta_success": {
        "status_code": 200,
        "text": SAMPLE_UNIPROT_FASTA["P35498"]
    },
    "uniprot_fasta_not_found": {
        "status_code": 404,
        "text": "Entry not found"
    },
    "uniprot_variants_success": {
        "status_code": 200,
        "json": SAMPLE_UNIPROT_VARIANTS["P35498"]
    },
    "ebi_job_submit_success": {
        "status_code": 200,
        "text": "clustalo-test-job-123"
    },
    "ebi_status_finished": {
        "status_code": 200,
        "text": "FINISHED"
    },
    "ebi_alignment_result": {
        "status_code": 200,
        "text": SAMPLE_EBI_ALIGNMENT
    },
    "ebi_tree_result": {
        "status_code": 200,
        "text": SAMPLE_EBI_TREE
    }
}