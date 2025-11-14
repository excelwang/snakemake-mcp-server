import click
import logging
import os
import sys
from ..mcp_factory import create_mcp_from_fastapi

logger = logging.getLogger(__name__)

@click.command(
    help="Start the Snakemake server with MCP protocol support. "
         "This provides MCP protocol endpoints derived from FastAPI definitions."
)
@click.option("--host", default="127.0.0.1", help="Host to bind to. Default: 127.0.0.1")
@click.option("--port", default=8083, type=int, help="Port to bind to. Default: 8083")
@click.option("--log-level", default="INFO", type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              help="Logging level. Default: INFO")
@click.pass_context
def mcp(ctx, host, port, log_level):
    """Start the Snakemake server with MCP protocol support."""
    # Reconfigure logging to respect the user's choice
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        force=True  # This is crucial to override the initial config
    )
    
    # Get paths from context (already strings now)
    wrappers_path = ctx.obj['WRAPPERS_PATH']
    workflows_dir = ctx.obj['WORKFLOWS_DIR']
    
    logger.setLevel(log_level)
    
    logger.info(f"Starting Snakemake Server with MCP protocol support...")
    logger.info(f"Server will be available at http://{host}:{port}")
    logger.info(f"MCP endpoints will be available at http://{host}:{port}/mcp")
    logger.info(f"Using snakebase from: {ctx.obj['SNAKEBASE_DIR']}")
    
    if not os.path.isdir(wrappers_path):
        logger.error(f"Wrappers directory not found at: {wrappers_path}")
        sys.exit(1)
    
    if not os.path.isdir(workflows_dir):
        logger.error(f"Workflows directory not found at: {workflows_dir}")
        sys.exit(1)

    mcp_server = create_mcp_from_fastapi(wrappers_path, workflows_dir)

    try:
        mcp_server.run(
            transport="http",
            host=host,
            port=port,
            log_level=log_level
        )
    except KeyboardInterrupt:
        logger.info("Server stopped by user")
    except Exception as e:
        logger.error(f"Server failed to start: {e}")
        sys.exit(1)
