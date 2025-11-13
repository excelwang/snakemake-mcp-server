import pytest
import os
import tempfile
import gzip
from snakemake_mcp_server.wrapper_runner import run_wrapper

@pytest.mark.asyncio
async def test_megahit_wrapper(wrappers_path):
    """Test the megahit wrapper with container, benchmark, resources, and shadow."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create dummy input files
        r1_path = os.path.join(temp_dir, "sample1_R1.fastq.gz")
        r2_path = os.path.join(temp_dir, "sample1_R2.fastq.gz")
        
        with gzip.open(r1_path, "wt") as f:
            f.write("@read1\nAGCT\n+\nIIII\n")
        with gzip.open(r2_path, "wt") as f:
            f.write("@read1\nAGCT\n+\nIIII\n")

        output_dir = os.path.join(temp_dir, "assembly")
        os.makedirs(output_dir)
        output_file = os.path.join("assembly", "final.contigs.fasta")
        benchmark_file = "benchmark.txt"

        result = await run_wrapper(
            wrapper_name="bio/megahit",
            wrappers_path=wrappers_path,
            inputs={
                "reads": ["sample1_R1.fastq.gz", "sample1_R2.fastq.gz"]
            },
            outputs={"contigs": output_file},
            params={"extra": "--min-count 10 --k-list 21,29,39,59,79,99,119,141"},
            container_img="docker://continuumio/miniconda3:4.4.10",
            benchmark=os.path.join(temp_dir, benchmark_file),
            resources={"mem_mb": 250000},
            workdir=temp_dir,
        )

        assert result["status"] == "success"
        assert result["exit_code"] == 0
        assert os.path.exists(os.path.join(temp_dir, output_file))
        assert os.path.exists(os.path.join(temp_dir, benchmark_file))
