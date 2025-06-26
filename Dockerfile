# syntax=docker/dockerfile:1
FROM python:3.10-slim

# Set work directory
WORKDIR /app

# Install system dependencies (if needed)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY pyproject.toml ./
COPY . .
RUN pip install --upgrade pip && pip install -e .

# Create data directories (if not present)
RUN mkdir -p /app/data/input /app/data/output /app/data/results

# Set environment variables (optional, for reproducibility)
ENV PYTHONUNBUFFERED=1

# Default command (can be overridden)
CMD ["comparative-genomics-pipeline"]
