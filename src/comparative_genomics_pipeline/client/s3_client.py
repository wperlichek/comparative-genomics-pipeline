import boto3
import logging
from typing import Optional
import json
from botocore.exceptions import ClientError, NoCredentialsError

from ..config.aws_config import AWSConfig

logger = logging.getLogger(__name__)


class S3Client:
    """
    Simple S3 client for caching UniProt/NCBI API responses.
    
    Usage:
        from ..config.aws_config import get_aws_config
        
        config = get_aws_config()
        client = S3Client(config)
        
        # Cache a FASTA sequence
        await client.put_sequence("P35498", fasta_content)
        
        # Retrieve cached sequence
        cached = await client.get_sequence("P35498")
    """
    
    def __init__(self, aws_config: AWSConfig):
        self.config = aws_config
        self.bucket_name = aws_config.s3.bucket_name
        self.enabled = True
        
        try:
            # Initialize boto3 client with optional profile
            session = boto3.Session(profile_name=aws_config.profile)
            self.s3_client = session.client('s3', region_name=aws_config.s3.region)
        except (ClientError, NoCredentialsError) as e:
            logger.warning(f"S3 caching disabled - AWS not configured: {e}")
            self.enabled = False
            self.s3_client = None
        
    async def get_sequence(self, accession_id: str) -> Optional[str]:
        """
        Retrieve cached FASTA sequence by UniProt accession ID.
        
        Args:
            accession_id: UniProt accession ID (e.g., "P35498")
            
        Returns:
            FASTA sequence string if found, None if not cached or S3 disabled
        """
        if not self.enabled:
            return None
            
        key = self.config.s3.get_sequence_key(accession_id)
        
        try:
            response = self.s3_client.get_object(Bucket=self.bucket_name, Key=key)
            content = response['Body'].read().decode('utf-8')
            logger.info(f"Retrieved cached sequence for {accession_id}")
            return content
            
        except ClientError as e:
            error_code = e.response['Error']['Code']
            if error_code == 'NoSuchKey':
                logger.debug(f"No cached sequence found for {accession_id}")
                return None
            else:
                logger.error(f"S3 error retrieving {accession_id}: {e}")
                return None
            
    async def put_sequence(self, accession_id: str, fasta_content: str) -> bool:
        """
        Cache a FASTA sequence in S3.
        
        Args:
            accession_id: UniProt accession ID
            fasta_content: FASTA sequence content
            
        Returns:
            True if successfully cached, False if failed or S3 disabled
        """
        if not self.enabled:
            return False
            
        key = self.config.s3.get_sequence_key(accession_id)
        
        try:
            self.s3_client.put_object(
                Bucket=self.bucket_name,
                Key=key,
                Body=fasta_content,
                ContentType='text/plain'
            )
            logger.info(f"Cached sequence for {accession_id}")
            return True
            
        except ClientError as e:
            logger.error(f"S3 error caching {accession_id}: {e}")
            return False
            
    async def get_variants(self, accession_id: str) -> Optional[dict]:
        """
        Retrieve cached variant data by UniProt accession ID.
        
        Args:
            accession_id: UniProt accession ID
            
        Returns:
            Variant data dict if found, None if not cached or S3 disabled
        """
        if not self.enabled:
            return None
            
        key = self.config.s3.get_variant_key(accession_id)
        
        try:
            response = self.s3_client.get_object(Bucket=self.bucket_name, Key=key)
            content = response['Body'].read().decode('utf-8')
            data = json.loads(content)
            logger.info(f"Retrieved cached variants for {accession_id}")
            return data
            
        except ClientError as e:
            error_code = e.response['Error']['Code']
            if error_code == 'NoSuchKey':
                logger.debug(f"No cached variants found for {accession_id}")
                return None
            else:
                logger.error(f"S3 error retrieving variants for {accession_id}: {e}")
                return None
                
        except json.JSONDecodeError as e:
            logger.error(f"Error parsing cached variants for {accession_id}: {e}")
            return None
            
    async def put_variants(self, accession_id: str, variant_data: dict) -> bool:
        """
        Cache variant data in S3.
        
        Args:
            accession_id: UniProt accession ID
            variant_data: Variant data dictionary
            
        Returns:
            True if successfully cached, False if failed or S3 disabled
        """
        if not self.enabled:
            return False
            
        key = self.config.s3.get_variant_key(accession_id)
        
        try:
            content = json.dumps(variant_data, indent=2)
            self.s3_client.put_object(
                Bucket=self.bucket_name,
                Key=key,
                Body=content,
                ContentType='application/json'
            )
            logger.info(f"Cached variants for {accession_id}")
            return True
            
        except ClientError as e:
            logger.error(f"S3 error caching variants for {accession_id}: {e}")
            return False
            
        except json.JSONEncodeError as e:
            logger.error(f"Error encoding variants for {accession_id}: {e}")
            return False
            
    async def list_cached_sequences(self) -> list[str]:
        """
        List all cached sequence accession IDs.
        
        Returns:
            List of accession IDs that have cached sequences, empty if S3 disabled
        """
        if not self.enabled:
            return []
            
        try:
            response = self.s3_client.list_objects_v2(
                Bucket=self.bucket_name,
                Prefix=self.config.s3.prefix + '/sequences/' if self.config.s3.prefix else 'sequences/'
            )
            
            if 'Contents' not in response:
                return []
                
            accession_ids = []
            for obj in response['Contents']:
                key = obj['Key']
                # Extract accession ID from "prefix/sequences/P35498.fasta"
                if key.endswith('.fasta'):
                    accession_id = key.split('/')[-1].replace('.fasta', '')
                    accession_ids.append(accession_id)
                    
            logger.info(f"Found {len(accession_ids)} cached sequences")
            return accession_ids
            
        except ClientError as e:
            logger.error(f"S3 error listing sequences: {e}")
            return []