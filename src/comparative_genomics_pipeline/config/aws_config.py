import os
from dataclasses import dataclass
from typing import Optional


@dataclass
class S3Config:
    """S3 configuration for sequence and variant caching."""
    
    bucket_name: str
    region: str = "us-east-1"
    prefix: str = "genomics-cache"
    
    # Key patterns
    sequence_key_pattern: str = "sequences/{accession_id}.fasta"
    variant_key_pattern: str = "variants/{accession_id}.json"
    
    @classmethod
    def from_env(cls) -> "S3Config":
        """Load S3 config from environment variables."""
        bucket_name = os.getenv("GENOMICS_S3_BUCKET")
        if not bucket_name:
            raise ValueError(
                "GENOMICS_S3_BUCKET environment variable is required. "
                "Set it to your S3 bucket name for caching."
            )
            
        return cls(
            bucket_name=bucket_name,
            region=os.getenv("AWS_REGION", "us-east-1"),
            prefix=os.getenv("GENOMICS_S3_PREFIX", "genomics-cache")
        )
    
    def get_sequence_key(self, accession_id: str) -> str:
        """Generate S3 key for a protein sequence."""
        key = self.sequence_key_pattern.format(accession_id=accession_id)
        return f"{self.prefix}/{key}" if self.prefix else key
    
    def get_variant_key(self, accession_id: str) -> str:
        """Generate S3 key for variant data."""
        key = self.variant_key_pattern.format(accession_id=accession_id)
        return f"{self.prefix}/{key}" if self.prefix else key


@dataclass 
class AWSConfig:
    """AWS configuration for the genomics pipeline."""
    
    s3: S3Config
    profile: Optional[str] = None
    
    @classmethod
    def from_env(cls) -> "AWSConfig":
        """Load AWS config from environment variables."""
        return cls(
            s3=S3Config.from_env(),
            profile=os.getenv("AWS_PROFILE")
        )


# Default configuration - can be overridden by environment variables
def get_aws_config() -> AWSConfig:
    """Get AWS configuration, preferring environment variables."""
    try:
        return AWSConfig.from_env()
    except ValueError as e:
        # Provide helpful error message if not configured
        raise ValueError(
            f"AWS configuration error: {e}\n\n"
            "Required environment variables:\n"
            "  GENOMICS_S3_BUCKET=your-bucket-name\n"
            "\n"
            "Optional environment variables:\n"
            "  AWS_REGION=us-east-1 (default)\n"
            "  AWS_PROFILE=your-profile (optional)\n"
            "  GENOMICS_S3_PREFIX=genomics-cache (default)\n"
        )